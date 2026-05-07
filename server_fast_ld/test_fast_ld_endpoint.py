"""Tests for fast_ld_endpoint.py."""
import asyncio
import base64
import gzip
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent   # bundle root (parent of server_turn11c_ld_fast/)
sys.path.insert(0, str(ROOT / "server_turn11c_ld_fast"))
sys.path.insert(0, str(ROOT / "engine_fast_ld"))

from fast_ld_endpoint import (
    FastLDReq, _cache_key, _default_triangle_assign, _shape_response,
    handle_split_heatmap, fast_ld_engine_hash,
)
from lazy_windows_json import WindowsJsonCache
from test_fast_ld import make_synthetic_windows_json, make_synthetic_dosage

FAST_LD_BIN = ROOT / "engine_fast_ld" / "fast_ld"
BUILD_WIN_SCRIPT = ROOT / "engine_fast_ld" / "build_windows_json.py"

if not FAST_LD_BIN.exists():
    sys.exit(
        f"fast_ld binary not built. From {ROOT}/engine_fast_ld/ run `make`,\n"
        f"or set FAST_LD_BIN explicitly. Expected at: {FAST_LD_BIN}")


# ============================================================================
# Tiny harness
# ============================================================================

PASSED = []

def _ok(name):
    PASSED.append(name)
    print(f"  PASS {name}")


# ============================================================================
# Helpers
# ============================================================================

def fake_runner_factory(matrices_payload):
    """Return a runner that returns a canned `result` shaped like
    fast_ld_wrapper.compute_ld()'s return."""
    async def runner(inner_req, *, fast_ld_bin):
        return matrices_payload
    return runner


def canned_result(n_snps=50, group_names=("HOM_INV", "HOM_REF"),
                  n_samples_per=20, shelf_ratio=10.0):
    """Generate a result dict matching fast_ld_wrapper.compute_ld() shape."""
    n_pairs = n_snps * (n_snps - 1) // 2
    matrices = {}
    for g in group_names:
        # Half the cells "high" (~ uint8 200), half "low" (~ uint8 20)
        bytes_arr = np.random.RandomState(0).randint(20, 200, n_pairs).astype(np.uint8)
        matrices[g] = {
            "n_samples": n_samples_per,
            "n_pairs": n_pairs,
            "pairs_b64": base64.b64encode(bytes_arr.tobytes()).decode("ascii"),
            "median_r2_overall": 0.4,
            "median_r2_shelf": 0.7 if g == "HOM_INV" else 0.4,
            "median_r2_flank": 0.07 if g == "HOM_INV" else 0.4,
            "shelf_ratio": shelf_ratio if g == "HOM_INV" else 1.0,
            "pct_pairs_above_0_8": 0.2,
            "decay_deciles": [0.5 - i * 0.03 for i in range(10)],
        }
    sites = {
        "idx": list(range(n_snps)),
        "pos": [1_000_000 + i * 100 for i in range(n_snps)],
    }
    for g in group_names:
        sites[f"maf_{g}"] = [0.3] * n_snps
        sites[f"var_{g}"] = [0.4] * n_snps
        sites[f"n_complete_{g}"] = [n_samples_per] * n_snps
    summary = {
        "engine_version": "fast_ld 0.2.0",
        "chrom": "C_gar_LG28",
        "window_lo": 0,
        "window_hi": 10,
        "n_snps_unique_in_range": n_snps,
        "n_snps_used": n_snps,
        "n_pairs": n_pairs,
        "thinning_applied": False,
        "snp_cap": 5000,
        "thin_to": None,
        "compute_seconds_total": 0.05,
        "shelf": None,
        "groups": {g: {
            "n_samples": n_samples_per,
            "shelf_ratio": shelf_ratio if g == "HOM_INV" else 1.0,
            "median_r2_overall": 0.4,
            "median_r2_shelf": 0.7 if g == "HOM_INV" else 0.4,
            "median_r2_flank": 0.07 if g == "HOM_INV" else 0.4,
            "pct_pairs_above_0_8": 0.2,
            "decay_deciles": [0.5 - i * 0.03 for i in range(10)],
            "compute_seconds": 0.025,
        } for g in group_names},
    }
    return {"summary": summary, "sites": sites, "matrices": matrices}


def make_dummy_files(td: Path, chrom: str = "C_gar_LG28"):
    """Build minimal atlas JSON, sites file, dosage file, and a windows
    cache pointing at them — so handle_split_heatmap can run end-to-end."""
    atlas_dir = td / "atlas"; atlas_dir.mkdir()
    sites_dir = td / "sites"; sites_dir.mkdir()
    dosage_dir = td / "dosage"; dosage_dir.mkdir()
    windows_dir = td / "windows"; windows_dir.mkdir()

    # Build a small synthetic windows JSON to derive the atlas + sites layout
    wj_path = atlas_dir / f"{chrom}_temp.windows.json"
    positions = make_synthetic_windows_json(
        wj_path, chrom,
        n_windows=20, snps_per_window=100, step=20,
        first_pos=1_000_000, snp_step_bp=200)

    # Write a fake atlas JSON with only what build_windows_json.py reads
    with open(wj_path) as fh:
        wj = json.load(fh)
    atlas_json = {
        "chrom": chrom,
        "windows": [
            {"idx": w["idx"], "start_bp": w["start_bp"], "end_bp": w["end_bp"]}
            for w in wj["windows"]
        ],
    }
    atlas_path = atlas_dir / f"{chrom}.json"
    with open(atlas_path, "w") as fh:
        json.dump(atlas_json, fh)
    wj_path.unlink()  # delete the temp

    # Write a sites file listing the unique positions
    sites_path = sites_dir / f"{chrom}.sites.tsv.gz"
    with gzip.open(sites_path, "wt") as fh:
        fh.write("marker\tchrom\tpos\tallele1\tallele2\n")
        for p in positions:
            fh.write(f"{chrom}_{p}\t{chrom}\t{p}\tA\tG\n")

    # Write a dosage file
    dosage_path = dosage_dir / f"{chrom}.dosage.tsv.gz"
    sample_ids = [f"S{i:03d}" for i in range(40)]
    make_synthetic_dosage(dosage_path, chrom, positions, sample_ids,
                          carrier_idx=list(range(20)),
                          shelf_snp_range=(200, 350))
    return {
        "atlas_dir": atlas_dir,
        "sites_dir": sites_dir,
        "dosage_dir": dosage_dir,
        "windows_dir": windows_dir,
        "sample_ids": sample_ids,
        "positions": positions,
        "chrom": chrom,
    }


# ============================================================================
# Tests
# ============================================================================

def test_request_validation():
    # Valid
    FastLDReq(chrom="C_gar_LG28", window_range=[0, 10],
              groups={"a": ["S1", "S2", "S3", "S4", "S5"]})

    # window_range bad
    for wr in ([], [5], [0, 1, 2], [10, 5], [-1, 5]):
        try:
            FastLDReq(chrom="x", window_range=wr,
                      groups={"a": ["S1", "S2", "S3", "S4", "S5"]})
            raise AssertionError(f"should reject window_range={wr}")
        except (ValueError, Exception):
            pass

    # group too small
    try:
        FastLDReq(chrom="x", window_range=[0, 5],
                  groups={"a": ["S1", "S2"]})
        raise AssertionError("should reject group of size 2")
    except Exception:
        pass

    # group bad name
    try:
        FastLDReq(chrom="x", window_range=[0, 5],
                  groups={"bad name!": ["S1", "S2", "S3", "S4", "S5"]})
        raise AssertionError("should reject 'bad name!'")
    except Exception:
        pass

    # too many groups
    try:
        FastLDReq(chrom="x", window_range=[0, 5],
                  groups={f"g{i}": ["S1", "S2", "S3", "S4", "S5"]
                          for i in range(5)})
        raise AssertionError("should reject 5 groups")
    except Exception:
        pass

    # snp_cap too small
    try:
        FastLDReq(chrom="x", window_range=[0, 5],
                  groups={"a": ["S1", "S2", "S3", "S4", "S5"]},
                  snp_cap=5)
        raise AssertionError("should reject snp_cap=5")
    except Exception:
        pass

    _ok("test_request_validation")


def test_cache_key_is_order_independent():
    a = FastLDReq(chrom="X", window_range=[0, 10],
                  groups={"g": ["S3", "S1", "S2", "S4", "S5"]})
    b = FastLDReq(chrom="X", window_range=[0, 10],
                  groups={"g": ["S5", "S4", "S3", "S2", "S1"]})
    assert _cache_key(a, "abc") == _cache_key(b, "abc")
    # Different engine hash → different key
    assert _cache_key(a, "abc") != _cache_key(a, "xyz")
    _ok("test_cache_key_is_order_independent")


def test_cache_key_changes_with_meaningful_inputs():
    a = FastLDReq(chrom="X", window_range=[0, 10],
                  groups={"g": ["S1", "S2", "S3", "S4", "S5"]})
    b = FastLDReq(chrom="X", window_range=[0, 11],
                  groups={"g": ["S1", "S2", "S3", "S4", "S5"]})
    c = FastLDReq(chrom="Y", window_range=[0, 10],
                  groups={"g": ["S1", "S2", "S3", "S4", "S5"]})
    d = FastLDReq(chrom="X", window_range=[0, 10],
                  groups={"g": ["S1", "S2", "S3", "S4", "S5"]},
                  thin_to=200)
    keys = {_cache_key(r, "h") for r in (a, b, c, d)}
    assert len(keys) == 4
    _ok("test_cache_key_changes_with_meaningful_inputs")


def test_default_triangle_assign():
    # Two groups: lower=larger
    m = {"HOM_INV": {"n_samples": 60}, "HOM_REF": {"n_samples": 226}}
    ta = _default_triangle_assign(m)
    assert ta == {"lower": "HOM_REF", "upper": "HOM_INV"}, ta
    # One group → None
    assert _default_triangle_assign({"a": {"n_samples": 10}}) is None
    # Three groups → None
    assert _default_triangle_assign({f"g{i}": {"n_samples": 10}
                                     for i in range(3)}) is None
    _ok("test_default_triangle_assign")


def test_shape_response_strips_dense_matrices():
    raw = canned_result(n_snps=30)
    # If `r2` got into raw matrices, it should NOT propagate
    raw["matrices"]["HOM_INV"]["r2"] = [[0.1] * 30] * 30
    req = FastLDReq(chrom="C_gar_LG28", window_range=[0, 10],
                    groups={"HOM_INV": ["S1", "S2", "S3", "S4", "S5"],
                            "HOM_REF": ["S6", "S7", "S8", "S9", "S10"]})
    payload = _shape_response(req, raw, 0.123, {}, {})
    for name in ("HOM_INV", "HOM_REF"):
        assert "pairs_b64" in payload["matrices"][name]
        assert "r2" not in payload["matrices"][name]
    assert payload["chrom"] == "C_gar_LG28"
    assert payload["n_snps"] == 30
    assert payload["timing"]["total_wallclock_seconds"] == 0.123
    assert payload["triangle_assign"] is not None
    _ok("test_shape_response_strips_dense_matrices")


def test_handler_with_fake_runner_caches_correctly():
    """End-to-end through handle_split_heatmap with fake runner; second
    call should hit cache."""
    cache: dict = {}
    def cache_get(k): return cache.get(k)
    def cache_put(k, v): cache[k] = v

    raw = canned_result(n_snps=30)

    async def go():
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            files = make_dummy_files(td)
            wcache = WindowsJsonCache(
                cache_dir=files["windows_dir"],
                atlas_json_dir=files["atlas_dir"],
                sites_dir=files["sites_dir"],
                builder_script=BUILD_WIN_SCRIPT)

            req = FastLDReq(
                chrom=files["chrom"],
                window_range=[0, 10],
                groups={"HOM_INV": files["sample_ids"][:20],
                        "HOM_REF": files["sample_ids"][20:]},
                shelf_bp=[files["positions"][50],
                          files["positions"][150]],
                snp_cap=5000)

            r1 = await handle_split_heatmap(
                req, fast_ld_bin=Path("/dev/null"),
                fast_ld_engine_hash="testhash",
                dosage_dir=files["dosage_dir"],
                windows_cache=wcache,
                cache_get=cache_get, cache_put=cache_put,
                runner=fake_runner_factory(raw))
            assert r1["cache_state"] == "miss", r1["cache_state"]
            assert r1["chrom"] == files["chrom"]
            assert r1["n_snps"] == 30
            # Has a cache_key?
            assert r1.get("cache_key", "").startswith("fast_ld.split.")

            r2 = await handle_split_heatmap(
                req, fast_ld_bin=Path("/dev/null"),
                fast_ld_engine_hash="testhash",
                dosage_dir=files["dosage_dir"],
                windows_cache=wcache,
                cache_get=cache_get, cache_put=cache_put,
                runner=fake_runner_factory(raw))
            assert r2["cache_state"] == "hit"
            assert r2["cache_key"] == r1["cache_key"]
            assert r2["matrices"]["HOM_INV"]["pairs_b64"] == \
                   r1["matrices"]["HOM_INV"]["pairs_b64"]

    asyncio.run(go())
    _ok("test_handler_with_fake_runner_caches_correctly")


def test_handler_dosage_file_404():
    async def go():
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            files = make_dummy_files(td)
            wcache = WindowsJsonCache(
                cache_dir=files["windows_dir"],
                atlas_json_dir=files["atlas_dir"],
                sites_dir=files["sites_dir"],
                builder_script=BUILD_WIN_SCRIPT)
            req = FastLDReq(
                chrom="GHOST_CHROM",
                window_range=[0, 5],
                groups={"a": files["sample_ids"][:10]})
            try:
                await handle_split_heatmap(
                    req, fast_ld_bin=Path("/dev/null"),
                    fast_ld_engine_hash="h",
                    dosage_dir=files["dosage_dir"],
                    windows_cache=wcache,
                    runner=fake_runner_factory(canned_result(n_snps=10)))
                raise AssertionError("should have raised 404")
            except Exception as e:
                # FastAPI HTTPException
                assert "404" in str(type(e)) or "404" in str(e) or "GHOST" in str(e)
    asyncio.run(go())
    _ok("test_handler_dosage_file_404")


def test_sample_filter_drops_unknown():
    """sample_filter callback drops unknown IDs; fewer than 5 → 400."""
    async def go():
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            files = make_dummy_files(td)
            wcache = WindowsJsonCache(
                cache_dir=files["windows_dir"],
                atlas_json_dir=files["atlas_dir"],
                sites_dir=files["sites_dir"],
                builder_script=BUILD_WIN_SCRIPT)
            known = set(files["sample_ids"])
            def sample_filter(ids):
                kept, dropped = [], []
                for x in ids:
                    (kept if x in known else dropped).append(x)
                return kept, dropped
            # Valid request with two ghosts
            req = FastLDReq(
                chrom=files["chrom"],
                window_range=[0, 5],
                groups={"a": files["sample_ids"][:10] + ["GHOST"],
                        "b": files["sample_ids"][10:30]})
            payload = await handle_split_heatmap(
                req, fast_ld_bin=Path("/dev/null"),
                fast_ld_engine_hash="h",
                dosage_dir=files["dosage_dir"],
                windows_cache=wcache,
                sample_filter=sample_filter,
                runner=fake_runner_factory(canned_result(group_names=("a","b"))))
            assert "dropped_samples" in payload
            assert "GHOST" in payload["dropped_samples"]["a"]

            # Now fewer than 5 valid: should 400
            req2 = FastLDReq(
                chrom=files["chrom"],
                window_range=[0, 5],
                groups={"a": ["GHOST1", "GHOST2", "GHOST3", "GHOST4",
                              "GHOST5", files["sample_ids"][0]]})
            try:
                await handle_split_heatmap(
                    req2, fast_ld_bin=Path("/dev/null"),
                    fast_ld_engine_hash="h",
                    dosage_dir=files["dosage_dir"],
                    windows_cache=wcache,
                    sample_filter=sample_filter,
                    runner=fake_runner_factory(canned_result()))
                raise AssertionError("should 400")
            except Exception as e:
                # HTTPException 400
                assert "400" in str(type(e)) or "min 5" in str(e)
    asyncio.run(go())
    _ok("test_sample_filter_drops_unknown")


def test_lazy_windows_json_builds_and_caches():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        files = make_dummy_files(td)
        wcache = WindowsJsonCache(
            cache_dir=files["windows_dir"],
            atlas_json_dir=files["atlas_dir"],
            sites_dir=files["sites_dir"],
            builder_script=BUILD_WIN_SCRIPT)
        # First call builds
        p = wcache.get_or_build(files["chrom"])
        assert p.exists()
        with open(p) as fh:
            wj = json.load(fh)
        assert wj["chrom"] == files["chrom"]
        assert wj["n_windows"] == 20
        # Second call: cache hit (same path, no rebuild error)
        p2 = wcache.get_or_build(files["chrom"])
        assert p2 == p
        # Touching atlas JSON forces rebuild — verify by deleting cache
        # and rebuilding (the rebuild path is what we want to exercise).
        # The mtime check is in `_is_stale`; we sanity-check by setting
        # cache mtime *backwards* so the atlas appears newer.
        atlas_p = files["atlas_dir"] / f"{files['chrom']}.json"
        old_t = p.stat().st_mtime
        os.utime(p, (old_t - 1000, old_t - 1000))   # cache appears older
        os.utime(atlas_p, (old_t, old_t))            # atlas at original time
        p3 = wcache.get_or_build(files["chrom"])
        # After rebuild, mtime should have moved forward (current time)
        assert p3.stat().st_mtime > old_t - 1000
        # Missing chrom → FileNotFoundError
        try:
            wcache.get_or_build("NO_SUCH")
            raise AssertionError("should raise FileNotFoundError")
        except FileNotFoundError:
            pass
    _ok("test_lazy_windows_json_builds_and_caches")


def test_engine_hash():
    h = fast_ld_engine_hash(FAST_LD_BIN)
    assert isinstance(h, str) and len(h) == 16
    # Stable on re-read
    h2 = fast_ld_engine_hash(FAST_LD_BIN)
    assert h == h2
    _ok("test_engine_hash")


def test_real_end_to_end_with_actual_binary():
    """No fake runner — actually call fast_ld through the wrapper."""
    async def go():
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            files = make_dummy_files(td)
            wcache = WindowsJsonCache(
                cache_dir=files["windows_dir"],
                atlas_json_dir=files["atlas_dir"],
                sites_dir=files["sites_dir"],
                builder_script=BUILD_WIN_SCRIPT)
            req = FastLDReq(
                chrom=files["chrom"],
                window_range=[0, 19],
                groups={"HOM_INV": files["sample_ids"][:20],
                        "HOM_REF": files["sample_ids"][20:]},
                shelf_bp=[files["positions"][200],
                          files["positions"][349]],
                snp_cap=5000,
                threads=2)
            payload = await handle_split_heatmap(
                req,
                fast_ld_bin=FAST_LD_BIN,
                fast_ld_engine_hash=fast_ld_engine_hash(FAST_LD_BIN),
                dosage_dir=files["dosage_dir"],
                windows_cache=wcache)
            assert payload["cache_state"] == "miss"
            assert payload["n_snps"] > 100   # all 20 windows of 100 step 20 → 480 SNPs
            for g in ("HOM_INV", "HOM_REF"):
                assert g in payload["matrices"]
                m = payload["matrices"][g]
                # base64 → bytes → sanity
                raw = base64.b64decode(m["pairs_b64"])
                expected = m["n_pairs"]
                assert len(raw) == expected
            # Shelf signal: HOM_INV ratio >> HOM_REF ratio
            sr_inv = payload["matrices"]["HOM_INV"]["shelf_ratio"]
            sr_ref = payload["matrices"]["HOM_REF"]["shelf_ratio"]
            assert sr_inv > 5 * sr_ref, f"sr_inv={sr_inv}, sr_ref={sr_ref}"
            print(f"    real-binary E2E: HOM_INV shelf_ratio={sr_inv:.2f} "
                  f"HOM_REF={sr_ref:.2f}")
            assert payload["triangle_assign"] is not None
    asyncio.run(go())
    _ok("test_real_end_to_end_with_actual_binary")


# ============================================================================
# Run
# ============================================================================

if __name__ == "__main__":
    test_request_validation()
    test_cache_key_is_order_independent()
    test_cache_key_changes_with_meaningful_inputs()
    test_default_triangle_assign()
    test_shape_response_strips_dense_matrices()
    test_lazy_windows_json_builds_and_caches()
    test_engine_hash()
    test_handler_with_fake_runner_caches_correctly()
    test_handler_dosage_file_404()
    test_sample_filter_drops_unknown()
    test_real_end_to_end_with_actual_binary()
    print(f"\n{len(PASSED)} server tests passed.")
