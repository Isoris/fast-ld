"""
fast_ld_wrapper.py — Python adapter between a JSON LD request and the
fast_ld C binary (window-mode).

Use:
    from fast_ld_wrapper import compute_ld

    result = compute_ld({
        "dosage_path": "/data/dosage_per_chrom/C_gar_LG28.dosage.tsv.gz",
        "windows_json": "/data/windows_per_chrom/C_gar_LG28.windows.json",
        "window_range": [2150, 2400],
        "groups": {
            "HOM_INV": ["S001", "S017", ...],
            "HOM_REF": ["S003", "S011", ...],
        },
        "shelf_bp": [15_200_000, 18_100_000],   # optional
        "snp_cap": 5000,                         # optional, default 5000
        "thin_to": None,                         # optional, override snp_cap
        "triangle_assign": {"lower": "HOM_REF", "upper": "HOM_INV"},
        "threads": 4,
    })

    # result["matrices"][groupname]["pairs_b64"]  → uint8 upper triangle, b64
    # result["matrices"][groupname]["r2"]         → dense float32 matrix
    # result["sites"]                              → per-SNP table dict
    # result["summary"]                            → numbers for caption
    # result["triangle_assign"]                    → passed through

The wrapper writes a control file, calls the binary, reads back binaries,
and returns a JSON-shaped result. It is the only entry point the Python
server should use.
"""

from __future__ import annotations

import base64
import math
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Module-level binary path — set this from the server config, or it will
# fall back to expecting `fast_ld` on $PATH.
FAST_LD_BIN: Optional[Path] = None

MAX_GROUPS = 4   # must match MAX_GROUPS in fast_ld.c
SITES_RECORD_DTYPE = np.dtype([
    ("idx", "<i4"),
    ("pos", "<i4"),
    ("maf", "<f4", (MAX_GROUPS,)),
    ("var", "<f4", (MAX_GROUPS,)),
    ("n_complete", "<i4", (MAX_GROUPS,)),
])
assert SITES_RECORD_DTYPE.itemsize == 8 + 12 * MAX_GROUPS


# ----- helpers ---------------------------------------------------------------

def _validate_request(req: Dict[str, Any]) -> None:
    """Raise ValueError on malformed requests before invoking C."""
    for k in ("dosage_path", "windows_json", "window_range", "groups"):
        if k not in req:
            raise ValueError(f"request missing required key: {k}")
    wr = req["window_range"]
    if not (isinstance(wr, (list, tuple)) and len(wr) == 2):
        raise ValueError("window_range must be [lo, hi]")
    lo, hi = int(wr[0]), int(wr[1])
    if lo < 0 or hi < lo:
        raise ValueError(f"invalid window_range: [{lo}, {hi}]")
    groups = req["groups"]
    if not groups or len(groups) > MAX_GROUPS:
        raise ValueError(
            f"groups must have 1..{MAX_GROUPS} entries (got {len(groups)})")
    for name, ids in groups.items():
        if not name.replace("_", "").isalnum():
            raise ValueError(
                f"group name '{name}' must be alphanumeric/underscore")
        if len(ids) < 5:
            raise ValueError(f"group '{name}' has only {len(ids)} samples (min 5)")
    if "shelf_bp" in req and req["shelf_bp"] is not None:
        sb = req["shelf_bp"]
        if not (isinstance(sb, (list, tuple)) and len(sb) == 2):
            raise ValueError("shelf_bp must be [start, end]")


def _write_control(path: Path, *, dosage_path: Path, windows_json: Path,
                   window_range: Tuple[int, int],
                   shelf_bp: Optional[Tuple[int, int]],
                   snp_cap: int, thin_to: Optional[int], threads: int,
                   out_dir: Path, groups: Dict[str, List[str]]) -> None:
    with open(path, "w") as fh:
        fh.write(f"dosage_path={dosage_path}\n")
        fh.write(f"windows_json={windows_json}\n")
        fh.write(f"window_range={window_range[0]}-{window_range[1]}\n")
        if shelf_bp is not None:
            fh.write(f"shelf_start={int(shelf_bp[0])}\n")
            fh.write(f"shelf_end={int(shelf_bp[1])}\n")
        fh.write(f"snp_cap={int(snp_cap)}\n")
        if thin_to is not None:
            fh.write(f"thin_to={int(thin_to)}\n")
        fh.write(f"out_dir={out_dir}\n")
        fh.write(f"threads={int(threads)}\n")
        for name, ids in groups.items():
            fh.write(f"group={name}:{','.join(ids)}\n")


def _read_summary(path: Path, group_names: List[str]) -> Dict[str, Any]:
    raw: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            k, _, v = line.partition("\t")
            raw[k] = v

    def num(key: str) -> Optional[float]:
        v = raw.get(key)
        if v is None:
            return None
        try:
            f = float(v)
            return None if (math.isnan(f) or math.isinf(f)) else f
        except ValueError:
            return None

    out: Dict[str, Any] = {
        "engine_version": raw.get("engine_version"),
        "chrom": raw.get("chrom"),
        "window_lo": int(raw.get("window_lo", -1)),
        "window_hi": int(raw.get("window_hi", -1)),
        "n_snps_unique_in_range": int(raw.get("n_snps_unique_in_range", 0)),
        "n_snps_used": int(raw.get("n_snps_used", 0)),
        "n_pairs": int(raw.get("n_pairs", 0)),
        "thinning_applied": int(raw.get("thinning_applied", 0)) == 1,
        "snp_cap": int(raw.get("snp_cap", 0)),
        "thin_to": int(raw.get("thin_to", 0)) or None,
        "compute_seconds_total": num("compute_seconds_total"),
        "shelf": None,
        "groups": {},
    }
    if int(raw.get("shelf_start", -1)) >= 0:
        out["shelf"] = {
            "start_bp": int(raw["shelf_start"]),
            "end_bp": int(raw["shelf_end"]),
        }
    for gn in group_names:
        out["groups"][gn] = {
            "n_samples": int(raw.get(f"{gn}.n_samples", 0)),
            "compute_seconds": num(f"{gn}.compute_seconds"),
            "median_r2_overall": num(f"{gn}.median_r2_overall"),
            "pct_pairs_above_0_8": num(f"{gn}.pct_pairs_above_0_8"),
            "median_r2_shelf": num(f"{gn}.median_r2_shelf"),
            "median_r2_flank": num(f"{gn}.median_r2_flank"),
            "shelf_ratio": num(f"{gn}.shelf_ratio"),
            "decay_deciles": [num(f"{gn}.decay_decile_{b}") for b in range(10)],
        }
    return out


def _read_sites(path: Path, n_snps: int,
                group_names: List[str]) -> Dict[str, Any]:
    arr = np.fromfile(path, dtype=SITES_RECORD_DTYPE)
    if len(arr) != n_snps:
        raise RuntimeError(f"sites.bin has {len(arr)} records, expected {n_snps}")
    out: Dict[str, Any] = {
        "idx": arr["idx"].tolist(),
        "pos": arr["pos"].tolist(),
    }
    for g, name in enumerate(group_names):
        out[f"maf_{name}"] = arr["maf"][:, g].tolist()
        out[f"var_{name}"] = arr["var"][:, g].tolist()
        out[f"n_complete_{name}"] = arr["n_complete"][:, g].tolist()
    return out


def _pairs_to_matrix(raw: bytes, n: int) -> np.ndarray:
    n_pairs = n * (n - 1) // 2
    if len(raw) != n_pairs:
        raise RuntimeError(
            f"pairs file has {len(raw)} bytes, expected {n_pairs}")
    arr = np.frombuffer(raw, dtype=np.uint8)
    M = np.full((n, n), np.nan, dtype=np.float32)
    iu = np.triu_indices(n, k=1)
    vals = arr.astype(np.float32) / 255.0
    M[iu] = vals
    M[(iu[1], iu[0])] = vals
    return M


# ----- public entry point ----------------------------------------------------

def compute_ld(req: Dict[str, Any], *,
               bin_path: Optional[Path] = None,
               scratch_dir: Optional[Path] = None,
               return_matrices: bool = True) -> Dict[str, Any]:
    """Run fast_ld on a request and return a JSON-shaped result.

    Parameters
    ----------
    req : dict
        Request as documented at module top.
    bin_path : Path, optional
        Path to the fast_ld binary. Defaults to module-level FAST_LD_BIN
        or `fast_ld` on $PATH.
    scratch_dir : Path, optional
        Where to put per-request scratch dirs (default: system temp).
    return_matrices : bool
        If True, dense float32 r² matrices are included in the result
        alongside the base64 raw bytes. If False, only base64 (smaller
        memory footprint when shipping straight to a browser).
    """
    _validate_request(req)
    binary = bin_path or FAST_LD_BIN or Path("fast_ld")
    binary = Path(binary)

    dosage_path = Path(req["dosage_path"])
    windows_json = Path(req["windows_json"])
    if not dosage_path.exists():
        raise FileNotFoundError(f"dosage file not found: {dosage_path}")
    if not windows_json.exists():
        raise FileNotFoundError(f"windows JSON not found: {windows_json}")

    window_range = (int(req["window_range"][0]), int(req["window_range"][1]))
    shelf_bp = None
    if "shelf_bp" in req and req["shelf_bp"] is not None:
        shelf_bp = (int(req["shelf_bp"][0]), int(req["shelf_bp"][1]))
    snp_cap = int(req.get("snp_cap", 5000))
    thin_to = req.get("thin_to")
    if thin_to is not None:
        thin_to = int(thin_to)
    threads = int(req.get("threads", 4))
    groups = req["groups"]
    group_names = list(groups.keys())

    with tempfile.TemporaryDirectory(prefix="fast_ld_",
                                     dir=str(scratch_dir) if scratch_dir
                                     else None) as td:
        td_p = Path(td)
        out_dir = td_p / "out"
        out_dir.mkdir()
        ctrl = td_p / "control.tsv"
        _write_control(ctrl,
                       dosage_path=dosage_path,
                       windows_json=windows_json,
                       window_range=window_range,
                       shelf_bp=shelf_bp,
                       snp_cap=snp_cap, thin_to=thin_to, threads=threads,
                       out_dir=out_dir, groups=groups)
        proc = subprocess.run([str(binary), str(ctrl)],
                              capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError(
                f"fast_ld failed (rc={proc.returncode}):\n{proc.stderr}")

        summary = _read_summary(out_dir / "summary.tsv", group_names)
        n_snps = summary["n_snps_used"]
        sites = _read_sites(out_dir / "sites.bin", n_snps, group_names)

        matrices: Dict[str, Any] = {}
        for name in group_names:
            raw_bytes = (out_dir / f"pairs.{name}.bin").read_bytes()
            entry: Dict[str, Any] = {
                "n_samples": summary["groups"][name]["n_samples"],
                "n_pairs": n_snps * (n_snps - 1) // 2,
                "pairs_b64": base64.b64encode(raw_bytes).decode("ascii"),
                "median_r2_overall": summary["groups"][name]["median_r2_overall"],
                "median_r2_shelf": summary["groups"][name]["median_r2_shelf"],
                "median_r2_flank": summary["groups"][name]["median_r2_flank"],
                "shelf_ratio": summary["groups"][name]["shelf_ratio"],
                "pct_pairs_above_0_8": summary["groups"][name]["pct_pairs_above_0_8"],
                "decay_deciles": summary["groups"][name]["decay_deciles"],
            }
            if return_matrices:
                entry["r2"] = _pairs_to_matrix(raw_bytes, n_snps).tolist()
            matrices[name] = entry

        return {
            "summary": summary,
            "sites": sites,
            "matrices": matrices,
            "triangle_assign": req.get("triangle_assign"),
            "stderr": proc.stderr,
        }
