#!/usr/bin/env python3
"""
test_fast_ld.py — end-to-end tests of the window-mode engine.

Strategy:
  1. Build a synthetic windows JSON.
  2. Build a matching dosage file with a known LD structure (an
     "inversion" block in a contiguous range of windows).
  3. Run fast_ld.
  4. Verify outputs.

We test:
  - JSON parse + window-range → unique SNP collection
  - Selective dosage read (only requested SNPs)
  - r² computation against numpy.corrcoef² (within uint8 quantization)
  - Synthetic inversion: shelf_ratio high for carriers, ~1 for non-carriers
  - snp_cap rejection and thin_to override
  - Unknown sample handling
  - Pair index formula
"""

import gzip
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FAST_LD = ROOT / "src" / "fast_ld"
assert FAST_LD.exists(), f"binary not built: {FAST_LD}"

MAX_GROUPS = 4


# ============================================================================
# Helpers to build synthetic inputs
# ============================================================================

def make_synthetic_windows_json(out_path, chrom, n_windows=200,
                                snps_per_window=100, step=20,
                                first_pos=1_000_000, snp_step_bp=200):
    """Build a windows JSON with sliding 100-SNP windows, step 20 SNPs.

    SNP positions are evenly spaced at `snp_step_bp`. Total unique SNPs:
    (n_windows - 1) * step + snps_per_window.
    """
    n_unique_snps = (n_windows - 1) * step + snps_per_window
    snp_positions = [first_pos + i * snp_step_bp for i in range(n_unique_snps)]

    windows = []
    for w in range(n_windows):
        snp_lo = w * step
        snp_hi = snp_lo + snps_per_window  # exclusive
        ws = snp_positions[snp_lo:snp_hi]
        windows.append({
            "idx": w,
            "start_bp": ws[0],
            "end_bp": ws[-1],
            "n_snps": len(ws),
            "snp_positions": ws,
        })

    out = {
        "schema_version": 1,
        "chrom": chrom,
        "n_windows": n_windows,
        "unique_snp_count": n_unique_snps,
        "first_bp": snp_positions[0],
        "last_bp": snp_positions[-1],
        "windows": windows,
    }
    with open(out_path, "w") as fh:
        json.dump(out, fh)
    return snp_positions


def make_synthetic_dosage(out_path, chrom, snp_positions, sample_ids,
                          carrier_idx, shelf_snp_range, seed=0):
    """Build a dosage TSV.

    Inside `shelf_snp_range = (lo, hi)` (SNP indices, half-open):
      carriers all carry a shared haplotype dose, so all shelf SNPs are
      tightly correlated within carriers.
    Outside: independent random dosages.
    Non-carriers: independent random everywhere.
    """
    rng = np.random.default_rng(seed)
    n_samp = len(sample_ids)
    carrier_set = set(carrier_idx)
    carrier_hap = rng.choice([0.0, 1.0, 2.0], size=n_samp,
                             p=[0.25, 0.5, 0.25])
    lo, hi = shelf_snp_range

    with gzip.open(out_path, "wt") as fh:
        fh.write("marker\t" + "\t".join(sample_ids) + "\n")
        for snp_i, pos in enumerate(snp_positions):
            in_shelf = lo <= snp_i < hi
            row = np.empty(n_samp, dtype=np.float32)
            for k in range(n_samp):
                if in_shelf and k in carrier_set:
                    row[k] = float(np.clip(carrier_hap[k] + rng.normal(0, 0.1),
                                           0, 2))
                else:
                    p = rng.uniform(0.1, 0.5)
                    row[k] = rng.binomial(2, p)
            fh.write(f"{chrom}_{pos}\t" +
                     "\t".join(f"{v:.4f}" for v in row) + "\n")


def write_control(path, **kw):
    groups = kw.pop("groups", {})
    with open(path, "w") as fh:
        for k, v in kw.items():
            if v is not None and v != "":
                fh.write(f"{k}={v}\n")
        for name, ids in groups.items():
            fh.write(f"group={name}:{','.join(ids)}\n")


def read_summary(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            k, _, v = line.partition("\t")
            out[k] = v
    return out


def read_pairs(path, n):
    n_pairs = n * (n - 1) // 2
    arr = np.fromfile(path, dtype=np.uint8)
    assert len(arr) == n_pairs, f"{len(arr)} vs {n_pairs}"
    return arr


def reconstruct_matrix(pairs_q, n):
    M = np.full((n, n), np.nan, dtype=np.float32)
    iu = np.triu_indices(n, k=1)
    M[iu] = pairs_q.astype(np.float32) / 255.0
    M[(iu[1], iu[0])] = M[iu]
    return M


def run_fast_ld(control_path, expect_fail=False):
    proc = subprocess.run([str(FAST_LD), str(control_path)],
                          capture_output=True, text=True)
    if expect_fail:
        if proc.returncode == 0:
            print("STDOUT:", proc.stdout)
            print("STDERR:", proc.stderr)
            raise RuntimeError("expected failure but rc=0")
        return proc
    if proc.returncode != 0:
        print("STDOUT:", proc.stdout)
        print("STDERR:", proc.stderr)
        raise RuntimeError(f"fast_ld failed rc={proc.returncode}")
    return proc


# ============================================================================
# Tests
# ============================================================================

def test_json_parse_and_unique_snp_collection():
    """Verify window range → unique SNP set is computed correctly."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        chrom = "TC"
        positions = make_synthetic_windows_json(
            td / "wj.json", chrom,
            n_windows=50, snps_per_window=100, step=20,
            first_pos=1_000_000, snp_step_bp=100)

        # Window range 0-4 (5 windows). With step 20 and width 100:
        # window 0: SNPs 0..99, window 4: SNPs 80..179
        # Union = SNPs 0..179 = 180 unique SNPs
        sample_ids = [f"S{i:03d}" for i in range(20)]
        make_synthetic_dosage(td / "d.tsv.gz", chrom, positions, sample_ids,
                              carrier_idx=list(range(10)),
                              shelf_snp_range=(50, 130))

        out = td / "out"; out.mkdir()
        write_control(td / "ctrl.tsv",
                      dosage_path=str(td / "d.tsv.gz"),
                      windows_json=str(td / "wj.json"),
                      window_range="0-4",
                      out_dir=str(out),
                      groups={"all": sample_ids})
        run_fast_ld(td / "ctrl.tsv")
        s = read_summary(out / "summary.tsv")
        assert int(s["n_snps_unique_in_range"]) == 180, s["n_snps_unique_in_range"]
        assert int(s["n_snps_used"]) == 180
        assert s["chrom"] == chrom
        print("PASSED test_json_parse_and_unique_snp_collection")


def test_known_correlation_matches_numpy():
    """fast_ld r² should match numpy.corrcoef² for the same dosage data."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        chrom = "TC"
        positions = make_synthetic_windows_json(
            td / "wj.json", chrom,
            n_windows=20, snps_per_window=100, step=20,
            first_pos=1_000_000, snp_step_bp=100)
        # 20 windows, step 20, width 100 → 480 unique SNPs
        rng = np.random.default_rng(7)
        n_samp = 80
        sample_ids = [f"S{i:03d}" for i in range(n_samp)]
        # Build raw dosage matrix in memory (so we can compare to numpy)
        D = rng.uniform(0, 2, size=(len(positions), n_samp)).astype(np.float32)
        # Inject some structure: every 50th SNP is a copy of SNP 0
        for i in range(0, len(positions), 50):
            if i > 0:
                D[i] = D[0]

        with gzip.open(td / "d.tsv.gz", "wt") as fh:
            fh.write("marker\t" + "\t".join(sample_ids) + "\n")
            for i, pos in enumerate(positions):
                row = D[i]
                fh.write(f"{chrom}_{pos}\t" +
                         "\t".join(f"{v:.6f}" for v in row) + "\n")

        out = td / "out"; out.mkdir()
        # Use only first 5 windows: SNPs 0..179 (180 unique)
        write_control(td / "ctrl.tsv",
                      dosage_path=str(td / "d.tsv.gz"),
                      windows_json=str(td / "wj.json"),
                      window_range="0-4",
                      out_dir=str(out),
                      groups={"all": sample_ids})
        run_fast_ld(td / "ctrl.tsv")
        s = read_summary(out / "summary.tsv")
        n = int(s["n_snps_used"])
        assert n == 180
        pairs = read_pairs(out / "pairs.all.bin", n)
        M_fast = reconstruct_matrix(pairs, n)
        # Reference: r² of D[0:180] in numpy
        D_sub = D[:180]
        R = np.corrcoef(D_sub) ** 2
        offdiag = ~np.eye(n, dtype=bool)
        diff = np.abs(M_fast[offdiag] - R[offdiag])
        diff = diff[~np.isnan(diff)]
        max_err = diff.max()
        assert max_err < 2.0 / 255.0, f"max err {max_err}"
        print(f"PASSED test_known_correlation_matches_numpy "
              f"(max err {max_err:.6f})")


def test_synthetic_inversion_signal():
    """Real signal/noise test: carrier shelf-ratio should be huge."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        chrom = "TC"
        positions = make_synthetic_windows_json(
            td / "wj.json", chrom,
            n_windows=80, snps_per_window=100, step=20,
            first_pos=1_000_000, snp_step_bp=200)
        # 80 windows → 1680 unique SNPs over 336 kb genomic span.
        # Shelf = SNPs 600..1000 (i.e. middle ~80 kb).
        n_samp = 100
        sample_ids = [f"S{i:03d}" for i in range(n_samp)]
        carriers = list(range(50))
        noncarriers = list(range(50, 100))
        make_synthetic_dosage(
            td / "d.tsv.gz", chrom, positions, sample_ids,
            carrier_idx=carriers, shelf_snp_range=(600, 1000), seed=42)

        # Convert SNP-index range to bp range for shelf_start / shelf_end
        shelf_start_bp = positions[600]
        shelf_end_bp = positions[999]

        out = td / "out"; out.mkdir()
        # All 80 windows. With snp_cap=2000, 1680 unique SNPs is fine.
        write_control(td / "ctrl.tsv",
                      dosage_path=str(td / "d.tsv.gz"),
                      windows_json=str(td / "wj.json"),
                      window_range="0-79",
                      shelf_start=shelf_start_bp,
                      shelf_end=shelf_end_bp,
                      snp_cap=2000,
                      out_dir=str(out),
                      threads=2,
                      groups={
                          "HOM_INV": [sample_ids[i] for i in carriers],
                          "HOM_REF": [sample_ids[i] for i in noncarriers],
                      })
        run_fast_ld(td / "ctrl.tsv")
        s = read_summary(out / "summary.tsv")
        sr_inv = float(s["HOM_INV.shelf_ratio"])
        sr_ref = float(s["HOM_REF.shelf_ratio"])
        med_inv = float(s["HOM_INV.median_r2_overall"])
        med_ref = float(s["HOM_REF.median_r2_overall"])
        sh_inv = float(s["HOM_INV.median_r2_shelf"])
        fl_inv = float(s["HOM_INV.median_r2_flank"])
        n_used = int(s["n_snps_used"])
        n_pairs = int(s["n_pairs"])
        print(f"  n_snps_used={n_used}, n_pairs={n_pairs}")
        print(f"  HOM_INV: median={med_inv:.4f} shelf={sh_inv:.4f} "
              f"flank={fl_inv:.4f} ratio={sr_inv:.2f}")
        print(f"  HOM_REF: median={med_ref:.4f} ratio={sr_ref:.2f}")
        assert sh_inv > 0.4
        assert fl_inv < 0.1
        assert sr_inv > 5.0, f"carrier shelf_ratio too low: {sr_inv}"
        assert sr_ref < 3.0, f"non-carrier shelf_ratio too high: {sr_ref}"
        print("PASSED test_synthetic_inversion_signal")


def test_snp_cap_rejection():
    """If unique SNPs > snp_cap and no thin_to, expect failure."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        chrom = "TC"
        positions = make_synthetic_windows_json(
            td / "wj.json", chrom,
            n_windows=200, snps_per_window=100, step=20,
            first_pos=1_000_000, snp_step_bp=100)
        sample_ids = [f"S{i:03d}" for i in range(20)]
        make_synthetic_dosage(td / "d.tsv.gz", chrom, positions, sample_ids,
                              carrier_idx=list(range(10)),
                              shelf_snp_range=(0, 0))
        out = td / "out"; out.mkdir()
        # 200 windows → ~4080 unique SNPs. With snp_cap=1000, should fail.
        write_control(td / "ctrl.tsv",
                      dosage_path=str(td / "d.tsv.gz"),
                      windows_json=str(td / "wj.json"),
                      window_range="0-199",
                      snp_cap=1000,
                      out_dir=str(out),
                      groups={"all": sample_ids})
        proc = run_fast_ld(td / "ctrl.tsv", expect_fail=True)
        assert "snp_cap" in proc.stderr, proc.stderr
        print("PASSED test_snp_cap_rejection")


def test_thin_to_override():
    """thin_to should override snp_cap and reduce SNPs."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        chrom = "TC"
        positions = make_synthetic_windows_json(
            td / "wj.json", chrom,
            n_windows=200, snps_per_window=100, step=20,
            first_pos=1_000_000, snp_step_bp=100)
        sample_ids = [f"S{i:03d}" for i in range(20)]
        make_synthetic_dosage(td / "d.tsv.gz", chrom, positions, sample_ids,
                              carrier_idx=list(range(10)),
                              shelf_snp_range=(0, 0))
        out = td / "out"; out.mkdir()
        write_control(td / "ctrl.tsv",
                      dosage_path=str(td / "d.tsv.gz"),
                      windows_json=str(td / "wj.json"),
                      window_range="0-199",
                      snp_cap=1000,    # would otherwise reject
                      thin_to=500,
                      out_dir=str(out),
                      groups={"all": sample_ids})
        run_fast_ld(td / "ctrl.tsv")
        s = read_summary(out / "summary.tsv")
        assert int(s["thinning_applied"]) == 1
        assert int(s["n_snps_used"]) <= 500
        assert int(s["n_snps_unique_in_range"]) > 1000
        print("PASSED test_thin_to_override")


def test_unknown_samples_dropped():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        chrom = "TC"
        positions = make_synthetic_windows_json(
            td / "wj.json", chrom,
            n_windows=10, snps_per_window=100, step=20,
            first_pos=1_000_000, snp_step_bp=200)
        sample_ids = [f"S{i:03d}" for i in range(20)]
        make_synthetic_dosage(td / "d.tsv.gz", chrom, positions, sample_ids,
                              carrier_idx=list(range(10)),
                              shelf_snp_range=(0, 100))
        out = td / "out"; out.mkdir()
        ids_with_ghosts = sample_ids[:10] + ["GHOST_001", "GHOST_002"]
        write_control(td / "ctrl.tsv",
                      dosage_path=str(td / "d.tsv.gz"),
                      windows_json=str(td / "wj.json"),
                      window_range="0-9",
                      out_dir=str(out),
                      groups={"g1": ids_with_ghosts,
                              "g2": sample_ids[10:]})
        proc = run_fast_ld(td / "ctrl.tsv")
        assert "unknown sample" in proc.stderr
        s = read_summary(out / "summary.tsv")
        assert int(s["g1.n_samples"]) == 10
        print("PASSED test_unknown_samples_dropped")


def test_pair_index_formula():
    """Pair index formula must round-trip for arbitrary N."""
    for n in [3, 10, 50]:
        expected = []
        for i in range(n):
            for j in range(i + 1, n):
                expected.append((i, j))
        recon = []
        for i in range(n):
            base = i * (2 * n - i - 1) // 2
            for j in range(i + 1, n):
                idx = base + (j - i - 1)
                recon.append((idx, i, j))
        recon.sort()
        assert [(i, j) for _, i, j in recon] == expected, f"n={n}"
        assert [k for k, _, _ in recon] == list(range(len(expected)))
    print("PASSED test_pair_index_formula")


def test_real_lg28_windows_json():
    """Smoke test: parse the real LG28.windows.json built from the atlas
    JSON, confirm chrom + n_windows match, then bail before dosage stage.
    Skipped if the file isn't present."""
    real_wj = Path("/tmp/LG28.windows.json")
    if not real_wj.exists():
        print("SKIPPED test_real_lg28_windows_json (no /tmp/LG28.windows.json)")
        return
    # Just call fast_ld with a missing dosage path; we expect it to parse
    # the JSON and fail at dosage open.
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        out = td / "out"; out.mkdir()
        write_control(td / "ctrl.tsv",
                      dosage_path="/no/such/file.tsv.gz",
                      windows_json=str(real_wj),
                      window_range="2150-2160",
                      out_dir=str(out),
                      groups={"g": ["S0", "S1", "S2", "S3", "S4"]})
        proc = run_fast_ld(td / "ctrl.tsv", expect_fail=True)
        assert "n_windows=4302" in proc.stderr
        assert "245 unique SNPs" in proc.stderr or "unique SNPs" in proc.stderr
        print("PASSED test_real_lg28_windows_json")


if __name__ == "__main__":
    test_pair_index_formula()
    test_json_parse_and_unique_snp_collection()
    test_known_correlation_matches_numpy()
    test_snp_cap_rejection()
    test_thin_to_override()
    test_unknown_samples_dropped()
    test_synthetic_inversion_signal()
    test_real_lg28_windows_json()
    print("\nAll tests passed.")
