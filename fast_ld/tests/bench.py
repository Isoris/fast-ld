"""Benchmark fast_ld at scales matching real catfish atlas use cases."""
import gzip
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FAST_LD = ROOT / "src" / "fast_ld"


def make_realistic_dataset(td, n_samples, n_windows, snps_per_window, step,
                           n_carriers, shelf_window_range):
    """Build a realistic dataset: sliding-window structure, inversion shelf."""
    chrom = "C_gar_LG28"
    snp_step_bp = 200  # ~5 SNPs/kb
    n_unique_snps = (n_windows - 1) * step + snps_per_window
    first_pos = 1_000_000
    snp_positions = [first_pos + i * snp_step_bp for i in range(n_unique_snps)]

    # Windows JSON
    windows = []
    for w in range(n_windows):
        snp_lo = w * step
        snp_hi = snp_lo + snps_per_window
        ws = snp_positions[snp_lo:snp_hi]
        windows.append({
            "idx": w, "start_bp": ws[0], "end_bp": ws[-1],
            "n_snps": len(ws), "snp_positions": ws,
        })
    wj_path = td / "wj.json"
    with open(wj_path, "w") as fh:
        json.dump({
            "schema_version": 1, "chrom": chrom, "n_windows": n_windows,
            "unique_snp_count": n_unique_snps,
            "first_bp": snp_positions[0], "last_bp": snp_positions[-1],
            "windows": windows,
        }, fh)

    # Dosage
    rng = np.random.default_rng(0)
    sample_ids = [f"S{i:04d}" for i in range(n_samples)]
    carriers = list(range(n_carriers))
    carrier_hap = rng.choice([0.0, 1.0, 2.0], size=n_carriers,
                             p=[0.25, 0.5, 0.25])
    shelf_lo, shelf_hi = shelf_window_range
    shelf_snp_lo = shelf_lo * step
    shelf_snp_hi = shelf_hi * step + snps_per_window
    dosage_path = td / "d.tsv.gz"

    t = time.time()
    with gzip.open(dosage_path, "wt", compresslevel=1) as fh:
        fh.write("marker\t" + "\t".join(sample_ids) + "\n")
        for snp_i, pos in enumerate(snp_positions):
            row = np.empty(n_samples, dtype=np.float32)
            in_shelf = shelf_snp_lo <= snp_i < shelf_snp_hi
            if in_shelf:
                row[:n_carriers] = np.clip(
                    carrier_hap + rng.normal(0, 0.1, n_carriers), 0, 2)
            else:
                p = rng.uniform(0.1, 0.5)
                row[:n_carriers] = rng.binomial(2, p, n_carriers)
            p2 = rng.uniform(0.1, 0.5)
            row[n_carriers:] = rng.binomial(2, p2, n_samples - n_carriers)
            fh.write(f"{chrom}_{pos}\t" +
                     "\t".join(f"{v:.4f}" for v in row) + "\n")
    print(f"  built dataset in {time.time()-t:.1f}s "
          f"(dosage {dosage_path.stat().st_size/1e6:.1f} MB)")

    return wj_path, dosage_path, sample_ids, carriers, snp_positions


def bench(label, n_samples=226, n_windows=4302,
          snps_per_window=100, step=20, n_carriers=113,
          window_range=(0, 4301), shelf_window_range=(2000, 2200),
          snp_cap=5000, thin_to=None, threads=4):
    print(f"\n=== {label} ===")
    print(f"  n_samples={n_samples}, n_windows={n_windows} ({snps_per_window}-SNP/step{step})")
    print(f"  window_range={window_range}, snp_cap={snp_cap}, thin_to={thin_to}, threads={threads}")
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        wj_path, dosage_path, sample_ids, carriers, positions = make_realistic_dataset(
            td, n_samples, n_windows, snps_per_window, step, n_carriers,
            shelf_window_range)
        out = td / "out"; out.mkdir()
        ctrl = td / "ctrl.tsv"
        with open(ctrl, "w") as fh:
            fh.write(f"dosage_path={dosage_path}\n")
            fh.write(f"windows_json={wj_path}\n")
            fh.write(f"window_range={window_range[0]}-{window_range[1]}\n")
            fh.write(f"snp_cap={snp_cap}\n")
            if thin_to:
                fh.write(f"thin_to={thin_to}\n")
            shelf_start_bp = positions[shelf_window_range[0] * step]
            shelf_end_bp = positions[(shelf_window_range[1] + 1) * step + snps_per_window - 1
                                     if (shelf_window_range[1] + 1) * step + snps_per_window <= len(positions)
                                     else -1]
            fh.write(f"shelf_start={shelf_start_bp}\n")
            fh.write(f"shelf_end={shelf_end_bp}\n")
            fh.write(f"out_dir={out}\n")
            fh.write(f"threads={threads}\n")
            fh.write(f"group=HOM_INV:{','.join(sample_ids[:n_carriers])}\n")
            fh.write(f"group=HOM_REF:{','.join(sample_ids[n_carriers:])}\n")

        t0 = time.time()
        proc = subprocess.run([str(FAST_LD), str(ctrl)],
                              capture_output=True, text=True)
        wall = time.time() - t0
        if proc.returncode != 0:
            print("  FAILED:")
            print(proc.stderr)
            return
        for line in proc.stderr.split("\n"):
            if any(s in line for s in ["unique SNPs", "computed in", "thin_to",
                                       "selective read", "total"]):
                print(f"  {line}")
        print(f"  wallclock: {wall:.2f}s")


if __name__ == "__main__":
    # 1. Tiny inversion: 4 windows, ~150 unique SNPs (the pathological case
    #    the windows-bin design failed at)
    bench("4-window inversion (the small case)",
          n_windows=300, snps_per_window=100, step=20,
          window_range=(100, 103),
          shelf_window_range=(101, 102),
          snp_cap=5000)

    # 2. Typical inversion candidate: 50 windows, ~1080 unique SNPs
    bench("50-window inversion candidate",
          n_windows=300, snps_per_window=100, step=20,
          window_range=(100, 149),
          shelf_window_range=(115, 135),
          snp_cap=5000)

    # 3. Large candidate: 250 windows, ~5080 unique SNPs (right at the cap)
    bench("250-window large region",
          n_windows=400, snps_per_window=100, step=20,
          window_range=(50, 299),
          shelf_window_range=(150, 200),
          snp_cap=8000)

    # 4. Whole chromosome with thin_to override
    bench("whole-chromosome with thin_to=2000",
          n_windows=4302, snps_per_window=100, step=20,
          window_range=(0, 4301),
          shelf_window_range=(2000, 2200),
          snp_cap=5000, thin_to=2000, threads=8)
