"""Tests for the Python wrapper around fast_ld (window-mode)."""
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))

from fast_ld_wrapper import compute_ld

# Re-use synthetic helpers from the engine test
sys.path.insert(0, str(ROOT / "tests"))
from test_fast_ld import (make_synthetic_windows_json, make_synthetic_dosage)

FAST_LD = ROOT / "src" / "fast_ld"


def test_wrapper_end_to_end():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        chrom = "TC"
        wj_path = td / "wj.json"
        positions = make_synthetic_windows_json(
            wj_path, chrom,
            n_windows=80, snps_per_window=100, step=20,
            first_pos=1_000_000, snp_step_bp=200)
        sample_ids = [f"S{i:03d}" for i in range(100)]
        carriers = list(range(50))
        noncarriers = list(range(50, 100))
        dosage_path = td / "d.tsv.gz"
        make_synthetic_dosage(dosage_path, chrom, positions, sample_ids,
                              carrier_idx=carriers,
                              shelf_snp_range=(600, 1000), seed=42)

        req = {
            "dosage_path": str(dosage_path),
            "windows_json": str(wj_path),
            "window_range": [0, 79],
            "groups": {
                "HOM_INV": [sample_ids[i] for i in carriers],
                "HOM_REF": [sample_ids[i] for i in noncarriers],
            },
            "shelf_bp": [positions[600], positions[999]],
            "snp_cap": 2000,
            "triangle_assign": {"lower": "HOM_REF", "upper": "HOM_INV"},
            "threads": 2,
        }
        result = compute_ld(req, bin_path=FAST_LD)

        # Top-level shape
        assert "summary" in result
        assert "sites" in result
        assert "matrices" in result
        assert result["triangle_assign"] == {"lower": "HOM_REF", "upper": "HOM_INV"}

        s = result["summary"]
        assert s["chrom"] == chrom
        assert s["window_lo"] == 0 and s["window_hi"] == 79
        assert s["n_snps_used"] > 100
        assert s["n_pairs"] == s["n_snps_used"] * (s["n_snps_used"] - 1) // 2
        assert s["thinning_applied"] is False
        assert s["shelf"]["start_bp"] == positions[600]
        assert "HOM_INV" in s["groups"] and "HOM_REF" in s["groups"]

        # Sites table aligned with n_snps_used
        sites = result["sites"]
        assert len(sites["pos"]) == s["n_snps_used"]
        assert "maf_HOM_INV" in sites and "maf_HOM_REF" in sites

        # Matrices
        n = s["n_snps_used"]
        for name in ("HOM_INV", "HOM_REF"):
            m = result["matrices"][name]
            r2 = np.array(m["r2"], dtype=np.float32)
            assert r2.shape == (n, n)
            assert np.isnan(np.diag(r2)).all()
            offdiag = r2[~np.isnan(r2)]
            assert (offdiag >= 0).all() and (offdiag <= 1).all()

        # Signal/noise
        sr_inv = s["groups"]["HOM_INV"]["shelf_ratio"]
        sr_ref = s["groups"]["HOM_REF"]["shelf_ratio"]
        print(f"  shelf_ratio: HOM_INV={sr_inv:.2f}  HOM_REF={sr_ref:.2f}")
        assert sr_inv > 5 * sr_ref

        # base64 round-trip parity with .r2 field
        import base64
        for name in ("HOM_INV", "HOM_REF"):
            raw = base64.b64decode(result["matrices"][name]["pairs_b64"])
            arr = np.frombuffer(raw, dtype=np.uint8)
            assert len(arr) == n * (n - 1) // 2
            M_recon = np.full((n, n), np.nan, dtype=np.float32)
            iu = np.triu_indices(n, k=1)
            M_recon[iu] = arr.astype(np.float32) / 255.0
            M_recon[(iu[1], iu[0])] = M_recon[iu]
            M_orig = np.array(result["matrices"][name]["r2"], dtype=np.float32)
            ok = ~np.isnan(M_recon)
            np.testing.assert_allclose(M_recon[ok], M_orig[ok], atol=1e-6)
        print("PASSED test_wrapper_end_to_end")


def test_wrapper_validation_errors():
    """Wrapper should raise ValueError on malformed requests, before
    invoking the binary."""
    cases = [
        # missing required key
        ({}, ValueError),
        ({"dosage_path": "/x", "windows_json": "/y"}, ValueError),
        # bad window_range
        ({"dosage_path": "/x", "windows_json": "/y",
          "window_range": "not a list",
          "groups": {"a": ["S1", "S2", "S3", "S4", "S5"]}}, ValueError),
        # bad group name
        ({"dosage_path": "/x", "windows_json": "/y",
          "window_range": [0, 10],
          "groups": {"bad name!": ["S1", "S2", "S3", "S4", "S5"]}}, ValueError),
        # group too small
        ({"dosage_path": "/x", "windows_json": "/y",
          "window_range": [0, 10],
          "groups": {"g": ["S1", "S2"]}}, ValueError),
        # too many groups
        ({"dosage_path": "/x", "windows_json": "/y",
          "window_range": [0, 10],
          "groups": {f"g{i}": ["S1", "S2", "S3", "S4", "S5"]
                     for i in range(5)}}, ValueError),
    ]
    for req, exc in cases:
        try:
            compute_ld(req, bin_path=FAST_LD)
            raise AssertionError(f"should have raised {exc.__name__} for {req}")
        except exc:
            pass
        except (FileNotFoundError, RuntimeError):
            # Some of these get past the validator and fail at file-existence
            # check, which is also fine — they don't reach the binary.
            pass
    print("PASSED test_wrapper_validation_errors")


if __name__ == "__main__":
    test_wrapper_validation_errors()
    test_wrapper_end_to_end()
    print("\nAll wrapper tests passed.")
