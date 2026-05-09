"""
fast_ld_endpoint.py — FastAPI endpoint for the windows-driven fast_ld engine.

Endpoint: POST /api/ld/split_heatmap

Request body (JSON):
    chrom            str
    window_range     [int, int]                   inclusive window indices
    groups           {name: [sample_id, ...], ...}    1..4 groups
    triangle_assign  {"lower": name, "upper": name}   optional, for split rendering
    shelf_bp         [int, int]                   optional, for shelf-ratio
    snp_cap          int = 5000                   reject if more unique SNPs
    thin_to          int                          optional: thin to N SNPs
    threads          int = 4

Response (JSON):
    chrom
    window_range
    n_snps                    matrix dimension N
    n_pairs                   = N*(N-1)/2
    sites                     {idx, pos, maf_<g>, var_<g>, n_complete_<g>, ...}
    matrices                  {<group_name>: {pairs_b64, n_samples,
                                              median_r2_overall, shelf_ratio,
                                              pct_pairs_above_0_8,
                                              decay_deciles, ...}}
    summary                   global summary block
    triangle_assign           pass-through
    timing                    {compute_seconds, total_wallclock_seconds}
    cache_state               'miss' | 'hit'

Cache
-----
Composite keyed: each (group, window_range, chrom, snp_cap, thin_to,
fast_ld_engine_hash) entry caches its uint8 pair bytes + per-group sites
slice. Two-group requests share the sites computation. Swapping the
upper-triangle group (HOM1 → HOM2) reuses the lower-triangle pair bytes.

The cache key encoding sorts sample IDs so order doesn't matter.
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import logging
from pathlib import Path
from typing import Any, Awaitable, Callable, Dict, List, Optional, Tuple

from fastapi import HTTPException
from pydantic import BaseModel, Field, field_validator

log = logging.getLogger(__name__)


# =============================================================================
# Request models
# =============================================================================

class _Range2(BaseModel):
    """A 2-tuple validated as a list."""
    @classmethod
    def coerce(cls, v: Any) -> Tuple[int, int]:
        if v is None:
            return None  # type: ignore
        if not isinstance(v, (list, tuple)) or len(v) != 2:
            raise ValueError("must be [lo, hi]")
        a, b = int(v[0]), int(v[1])
        return (a, b)


class FastLDReq(BaseModel):
    chrom: str
    window_range: List[int]                      # [lo, hi] inclusive
    groups: Dict[str, List[str]]                 # 1..4 named groups
    triangle_assign: Optional[Dict[str, str]] = None
    shelf_bp: Optional[List[int]] = None         # [start, end]
    snp_cap: int = 5000
    thin_to: Optional[int] = None
    threads: int = 4

    @field_validator("window_range")
    @classmethod
    def _check_window_range(cls, v: List[int]) -> List[int]:
        if len(v) != 2:
            raise ValueError("window_range must be [lo, hi]")
        if v[0] < 0 or v[1] < v[0]:
            raise ValueError(f"invalid window_range: {v}")
        return v

    @field_validator("shelf_bp")
    @classmethod
    def _check_shelf(cls, v: Optional[List[int]]) -> Optional[List[int]]:
        if v is None:
            return None
        if len(v) != 2 or v[0] >= v[1]:
            raise ValueError("shelf_bp must be [start, end] with start < end")
        return v

    @field_validator("groups")
    @classmethod
    def _check_groups(cls, v: Dict[str, List[str]]) -> Dict[str, List[str]]:
        if not v:
            raise ValueError("groups must not be empty")
        if len(v) > 4:
            raise ValueError("at most 4 groups")
        for name, ids in v.items():
            if not name or not name.replace("_", "").isalnum():
                raise ValueError(f"group name '{name}' must be alphanumeric/underscore")
            if not ids or len(ids) < 5:
                raise ValueError(f"group '{name}' has only {len(ids)} samples (min 5)")
        return v

    @field_validator("snp_cap")
    @classmethod
    def _check_cap(cls, v: int) -> int:
        if v < 10:
            raise ValueError("snp_cap too small (min 10)")
        return v


# =============================================================================
# Cache key
# =============================================================================

def _canonicalize_groups(groups: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """Sort sample IDs within each group so cache key is order-independent."""
    return {name: sorted(ids) for name, ids in groups.items()}


def _cache_key(req: FastLDReq, engine_hash: str) -> str:
    body = {
        "chrom": req.chrom,
        "window_range": list(req.window_range),
        "groups": _canonicalize_groups(req.groups),
        "shelf_bp": req.shelf_bp,
        "snp_cap": req.snp_cap,
        "thin_to": req.thin_to,
        "engine": engine_hash,
    }
    encoded = json.dumps(body, sort_keys=True, separators=(",", ":"))
    h = hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:24]
    return f"fast_ld.split.{h}"


# =============================================================================
# Public handler
# =============================================================================

async def handle_split_heatmap(
    req: FastLDReq,
    *,
    fast_ld_bin: Path,
    fast_ld_engine_hash: str,
    dosage_dir: Path,
    windows_cache,                  # WindowsJsonCache instance
    sample_filter: Optional[Callable[[List[str]], Tuple[List[str], List[str]]]] = None,
    cache_get: Optional[Callable[[str], Optional[Dict[str, Any]]]] = None,
    cache_put: Optional[Callable[[str, Dict[str, Any]], None]] = None,
    runner: Optional[Callable[..., Awaitable[Dict[str, Any]]]] = None,
) -> Dict[str, Any]:
    """Top-level handler called by FastAPI.

    Parameters
    ----------
    req : FastLDReq
        Validated request.
    fast_ld_bin : Path
        Path to the fast_ld C binary.
    fast_ld_engine_hash : str
        Content-hash of the binary, used in cache keys.
    dosage_dir : Path
        Directory containing `<chrom>.dosage.tsv.gz` files.
    windows_cache : WindowsJsonCache
        Lazy windows JSON resolver.
    sample_filter : callable
        Optional `(ids) -> (kept, dropped)` to drop unknown sample IDs
        upfront (mirrors popstats_groupwise idiom).
    cache_get / cache_put : optional
        Cache layer. If None, no caching.
    runner : optional async callable
        Override for `_run_fast_ld()`. Used by tests to inject a fake.
    """
    import time

    t0 = time.monotonic()

    # 1. Cache lookup
    key = _cache_key(req, fast_ld_engine_hash)
    if cache_get is not None:
        cached = cache_get(key)
        if cached is not None:
            cached["cache_state"] = "hit"
            cached["cache_key"] = key
            cached["timing"] = {
                **cached.get("timing", {}),
                "total_wallclock_seconds": round(time.monotonic() - t0, 4),
            }
            return cached

    # 2. Resolve dosage path & windows JSON
    dosage_path = Path(dosage_dir) / f"{req.chrom}.dosage.tsv.gz"
    if not dosage_path.exists():
        raise HTTPException(404,
            f"dosage file not found for chrom '{req.chrom}': {dosage_path}")
    try:
        windows_json = windows_cache.get_or_build(req.chrom)
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))
    except RuntimeError as e:
        raise HTTPException(500, str(e))

    # 3. Filter samples through canonical master list (drops unknowns
    #    with a warning rather than failing the whole request).
    groups_filtered: Dict[str, List[str]] = {}
    drops_summary: Dict[str, List[str]] = {}
    for name, ids in req.groups.items():
        if sample_filter is not None:
            kept, dropped = sample_filter(ids)
            if dropped:
                drops_summary[name] = dropped[:10]   # cap reported drops
        else:
            kept = ids
        if len(kept) < 5:
            raise HTTPException(400,
                f"group '{name}' has only {len(kept)} known samples "
                f"after filtering (min 5)")
        groups_filtered[name] = kept

    # 4. Build the request for fast_ld_wrapper
    inner_req: Dict[str, Any] = {
        "dosage_path": str(dosage_path),
        "windows_json": str(windows_json),
        "window_range": list(req.window_range),
        "groups": groups_filtered,
        "snp_cap": req.snp_cap,
        "threads": req.threads,
    }
    if req.thin_to is not None:
        inner_req["thin_to"] = req.thin_to
    if req.shelf_bp is not None:
        inner_req["shelf_bp"] = list(req.shelf_bp)
    if req.triangle_assign is not None:
        inner_req["triangle_assign"] = req.triangle_assign

    # 5. Run fast_ld
    if runner is None:
        runner = _run_fast_ld_default
    try:
        result = await runner(inner_req, fast_ld_bin=fast_ld_bin)
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))
    except ValueError as e:
        raise HTTPException(400, str(e))
    except RuntimeError as e:
        raise HTTPException(500, str(e))

    # 6. Shape the response
    wallclock = round(time.monotonic() - t0, 4)
    payload = _shape_response(req, result, wallclock,
                              drops_summary, groups_filtered)
    payload["cache_state"] = "miss"
    payload["cache_key"] = key

    if cache_put is not None:
        # Store everything except cache_state (set fresh on each access).
        to_store = {k: v for k, v in payload.items()
                    if k not in ("cache_state",)}
        try:
            cache_put(key, to_store)
        except Exception as e:
            log.warning("cache_put failed for %s: %s", key, e)

    return payload


def _shape_response(req: FastLDReq, raw: Dict[str, Any],
                    wallclock_s: float,
                    drops_summary: Dict[str, List[str]],
                    groups_filtered: Dict[str, List[str]]) -> Dict[str, Any]:
    """Convert fast_ld_wrapper output to API response shape."""
    summary = raw["summary"]
    sites = raw["sites"]
    matrices_in = raw["matrices"]
    # Strip the dense float matrices — over the wire we only need pairs_b64
    # plus per-group summary numbers. The atlas decodes pairs_b64 client-side.
    matrices_out: Dict[str, Any] = {}
    for name, m in matrices_in.items():
        matrices_out[name] = {
            "n_samples": m["n_samples"],
            "n_pairs": m["n_pairs"],
            "pairs_b64": m["pairs_b64"],
            "median_r2_overall": m["median_r2_overall"],
            "median_r2_shelf": m["median_r2_shelf"],
            "median_r2_flank": m["median_r2_flank"],
            "shelf_ratio": m["shelf_ratio"],
            "pct_pairs_above_0_8": m["pct_pairs_above_0_8"],
            "decay_deciles": m["decay_deciles"],
        }

    triangle_assign = req.triangle_assign or _default_triangle_assign(matrices_out)

    payload = {
        "chrom": req.chrom,
        "window_range": list(req.window_range),
        "n_snps": summary["n_snps_used"],
        "n_pairs": summary["n_pairs"],
        "shelf_bp": req.shelf_bp,
        "thinning_applied": summary["thinning_applied"],
        "n_snps_unique_in_range": summary["n_snps_unique_in_range"],
        "snp_cap": summary["snp_cap"],
        "thin_to": summary["thin_to"],
        "sites": sites,
        "matrices": matrices_out,
        "triangle_assign": triangle_assign,
        "summary": summary,
        "timing": {
            "compute_seconds": summary.get("compute_seconds_total"),
            "total_wallclock_seconds": wallclock_s,
        },
        "engine_version": summary.get("engine_version"),
    }
    if drops_summary:
        payload["dropped_samples"] = drops_summary
    return payload


def _default_triangle_assign(matrices: Dict[str, Any]) -> Optional[Dict[str, str]]:
    """If exactly two groups, default to lower=largest, upper=other.
    Otherwise None — caller picks manually."""
    if len(matrices) != 2:
        return None
    names = list(matrices.keys())
    a, b = names[0], names[1]
    na = matrices[a]["n_samples"]
    nb = matrices[b]["n_samples"]
    if na >= nb:
        return {"lower": a, "upper": b}
    return {"lower": b, "upper": a}


# =============================================================================
# Default runner — invokes fast_ld_wrapper.compute_ld
# =============================================================================

async def _run_fast_ld_default(req: Dict[str, Any], *,
                               fast_ld_bin: Path) -> Dict[str, Any]:
    """Call the synchronous wrapper in a thread so we don't block the
    event loop."""
    # Lazy import so the module is importable for tests that inject a runner.
    from fast_ld_wrapper import compute_ld
    loop = asyncio.get_running_loop()
    def _go() -> Dict[str, Any]:
        return compute_ld(req, bin_path=fast_ld_bin, return_matrices=False)
    return await loop.run_in_executor(None, _go)


# =============================================================================
# Engine hash helper (used by /api/health and cache keys)
# =============================================================================

def fast_ld_engine_hash(bin_path: Path) -> str:
    """sha256 of the binary, first 16 hex chars."""
    h = hashlib.sha256()
    with open(bin_path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]
