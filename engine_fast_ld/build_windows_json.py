#!/usr/bin/env python3
"""
build_windows_json.py — preprocessing step that builds a slim per-chrom
windows JSON from the atlas's chromosome JSON + the canonical sites file.

Produces a lookup table:
    window_idx → {start_bp, end_bp, snp_positions}

Where snp_positions is the list of SNP bp positions (from sites.tsv.gz,
which is paired with the dosage file from STEP_A01) that fall inside the
window's [start_bp, end_bp].

Usage
-----
    python3 build_windows_json.py \\
        --atlas-json /path/to/LG28.json \\
        --sites /path/to/LG28.sites.tsv.gz \\
        --out /path/to/LG28.windows.json

Notes
-----
* Reads windows from atlas JSON's `windows` array (each entry has
  `idx`, `start_bp`, `end_bp`).
* Reads SNP positions from a sites file. Two formats supported:
    1. STEP_A01 output: gzipped TSV with header
       `marker, chrom, pos, allele1, allele2`
    2. Upstream ANGSD sites: 4-col TSV (no header) `chrom, pos, major, minor`
* Each window's snp_positions list is built by binary-searching the sorted
  SNP position array for `[start_bp, end_bp]` (inclusive on both ends).
* Windows overlap by construction (sliding 100-SNP windows with step 20),
  so the same SNP appears in multiple windows. That is expected.
* Output JSON also carries `unique_snp_count` (deduped across all windows)
  for downstream sanity checks.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import sys
import time
from pathlib import Path
from typing import List, Tuple


def read_sites_file(sites_path: Path, target_chrom: str) -> List[int]:
    """Read SNP positions for `target_chrom` from a sites file.

    Auto-detects format:
      - STEP_A01: gzipped TSV with header starting with 'marker'
      - ANGSD: plain or gzipped TSV with no header, 4 columns
    """
    is_gz = sites_path.suffix == ".gz"
    opener = gzip.open if is_gz else open
    positions: List[int] = []
    n_total = 0
    with opener(sites_path, "rt") as fh:
        first = fh.readline()
        if not first:
            raise RuntimeError(f"empty sites file: {sites_path}")
        # Detect header by checking if first field looks like a column name
        is_step_a01 = first.startswith("marker\t") or first.startswith("marker ")
        if not is_step_a01:
            # First line is data — process it
            _process_line(first, target_chrom, positions)
            n_total += 1
        for line in fh:
            _process_line(line, target_chrom, positions)
            n_total += 1
    if not positions:
        raise RuntimeError(
            f"no SNPs matching chrom='{target_chrom}' in {sites_path} "
            f"(read {n_total} lines)")
    # Sort & dedupe — sites files SHOULD be sorted but be defensive
    positions.sort()
    deduped: List[int] = []
    for p in positions:
        if not deduped or deduped[-1] != p:
            deduped.append(p)
    return deduped


def _process_line(line: str, target_chrom: str, out: List[int]) -> None:
    line = line.rstrip("\n")
    if not line:
        return
    fields = line.split("\t")
    if len(fields) < 2:
        return
    # Determine columns:
    # STEP_A01 fmt: marker chrom pos allele1 allele2  (5 cols)
    # ANGSD fmt:    chrom pos major minor              (4 cols)
    if len(fields) >= 5:
        chrom = fields[1]
        pos_field = fields[2]
    else:
        chrom = fields[0]
        pos_field = fields[1]
    if chrom != target_chrom:
        return
    try:
        out.append(int(pos_field))
    except ValueError:
        return  # skip malformed


def build_windows(atlas_json_path: Path, sites_path: Path,
                  out_path: Path) -> None:
    t0 = time.time()
    with open(atlas_json_path) as fh:
        atlas = json.load(fh)
    chrom = atlas["chrom"]
    atlas_windows = atlas["windows"]
    n_windows = len(atlas_windows)
    print(f"[build_windows_json] chrom={chrom}, n_windows={n_windows}",
          file=sys.stderr)

    print(f"[build_windows_json] reading sites from {sites_path}",
          file=sys.stderr)
    snp_positions = read_sites_file(sites_path, chrom)
    print(f"[build_windows_json] loaded {len(snp_positions):,} SNP positions "
          f"on {chrom}", file=sys.stderr)

    # For each window, binary-search the SNP positions for [start_bp, end_bp].
    out_windows = []
    total_membership = 0
    n_empty = 0
    snps_per_win_counts = []
    for w in atlas_windows:
        start_bp = int(w["start_bp"])
        end_bp = int(w["end_bp"])
        lo = bisect.bisect_left(snp_positions, start_bp)
        hi = bisect.bisect_right(snp_positions, end_bp)
        snps_in_win = snp_positions[lo:hi]
        if not snps_in_win:
            n_empty += 1
        snps_per_win_counts.append(len(snps_in_win))
        total_membership += len(snps_in_win)
        out_windows.append({
            "idx": int(w["idx"]),
            "start_bp": start_bp,
            "end_bp": end_bp,
            "n_snps": len(snps_in_win),
            "snp_positions": snps_in_win,
        })

    if n_empty > 0:
        print(f"[build_windows_json] WARNING: {n_empty} windows contain "
              f"zero SNPs (chromosomally unrepresented in sites file)",
              file=sys.stderr)

    # Sanity statistics for the user
    if snps_per_win_counts:
        avg = sum(snps_per_win_counts) / len(snps_per_win_counts)
        sorted_counts = sorted(snps_per_win_counts)
        med = sorted_counts[len(sorted_counts) // 2]
        mn = sorted_counts[0]
        mx = sorted_counts[-1]
        print(f"[build_windows_json] SNPs per window: "
              f"min={mn}, median={med}, mean={avg:.1f}, max={mx}",
              file=sys.stderr)
        print(f"[build_windows_json] total membership entries: "
              f"{total_membership:,} (sum across windows; SNPs counted "
              f"once per window they belong to)", file=sys.stderr)

    out = {
        "schema_version": 1,
        "chrom": chrom,
        "n_windows": n_windows,
        "unique_snp_count": len(snp_positions),
        "first_bp": snp_positions[0] if snp_positions else None,
        "last_bp": snp_positions[-1] if snp_positions else None,
        "atlas_json_source": str(atlas_json_path),
        "sites_file_source": str(sites_path),
        "windows": out_windows,
    }
    with open(out_path, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    sz = out_path.stat().st_size
    print(f"[build_windows_json] wrote {out_path} ({sz/1e6:.2f} MB) "
          f"in {time.time()-t0:.2f}s", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atlas-json", required=True, type=Path,
                    help="Atlas chromosome JSON (e.g. LG28.json)")
    ap.add_argument("--sites", required=True, type=Path,
                    help="Sites file: gzipped STEP_A01 sites.tsv.gz OR "
                         "upstream sites.thin_W.majmin.tsv")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output windows JSON path")
    args = ap.parse_args()

    if not args.atlas_json.exists():
        sys.exit(f"missing atlas JSON: {args.atlas_json}")
    if not args.sites.exists():
        sys.exit(f"missing sites file: {args.sites}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    build_windows(args.atlas_json, args.sites, args.out)


if __name__ == "__main__":
    main()
