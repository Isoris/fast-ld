#!/usr/bin/env python3
"""
STEP_LD_01_build_binary_cache.py  (v4)

Parse raw ngsLD pairs_r2.tsv ONCE, bin into symmetric numpy matrix,
save as pickle cache. Supports three cache types:
  - localdense: from full-marker LD (08_ld_pruned81), short-range
  - chrsparse:  from thinned-marker LD (08_ld_chrsparse), long-range
  - invzoom:    candidate-specific zoom LD

Usage:
    python3 STEP_LD_01_build_binary_cache.py \
        --infile <pairs_r2.tsv> --chrom <name> --chr-len <bp> \
        --bin-bp <size> --cache-type invzoom --outdir <dir> \
        [--outfile <exact_cache_path>] \
        [--recreate-binary-matrix]
"""
import argparse
import os
import sys
import pickle
import time
import numpy as np


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--chr-len", type=int, required=True)
    ap.add_argument("--bin-bp", type=int, required=True)
    ap.add_argument("--cache-type", default="localdense", choices=["localdense", "chrsparse", "invzoom"])
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--outfile", default=None, help="Exact output cache path to write")
    ap.add_argument("--recreate-binary-matrix", action="store_true")
    ap.add_argument("--use-previous-binary", action="store_true")
    return ap.parse_args()


def nice_bin_label(bp):
    return f"{bp//1000}kb" if bp < 1_000_000 else f"{bp//1_000_000}Mb"


def parse_raw(infile):
    t0 = time.time()
    p1, p2, r2 = [], [], []
    nbad = 0

    with open(infile) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                nbad += 1
                continue
            try:
                a = int(parts[0].rsplit(":", 1)[1])
                b = int(parts[1].rsplit(":", 1)[1])
                v = float(parts[3])
            except (ValueError, IndexError):
                nbad += 1
                continue

            if v < 0 or v > 1 or not np.isfinite(v):
                continue

            p1.append(a)
            p2.append(b)
            r2.append(v)

    dt = time.time() - t0
    print(f"[cache] Parsed {len(p1):,} pairs in {dt:.1f}s ({nbad} bad)", file=sys.stderr)
    return np.array(p1, np.int64), np.array(p2, np.int64), np.array(r2, np.float32)


def build_matrix(p1, p2, r2, chr_len, bin_bp):
    n = (chr_len + bin_bp - 1) // bin_bp

    rsum = np.zeros((n, n), np.float64)
    cnt = np.zeros((n, n), np.int32)

    mask = (p1 >= 1) & (p2 >= 1) & (p1 <= chr_len) & (p2 <= chr_len)

    b1 = ((p1[mask] - 1) // bin_bp).astype(np.int32)
    b2 = ((p2[mask] - 1) // bin_bp).astype(np.int32)

    bx = np.minimum(b1, b2)
    by = np.maximum(b1, b2)

    np.add.at(rsum, (bx, by), r2[mask].astype(np.float64))
    np.add.at(cnt, (bx, by), 1)

    with np.errstate(divide="ignore", invalid="ignore"):
        mean_r2 = np.where(cnt > 0, rsum / cnt, np.nan).astype(np.float32)

    tri = np.triu_indices(n, k=1)
    mean_r2[tri[1], tri[0]] = mean_r2[tri[0], tri[1]]
    cnt[tri[1], tri[0]] = cnt[tri[0], tri[1]]

    print(f"[cache] Matrix {n}×{n}, {np.count_nonzero(cnt):,} filled tiles", file=sys.stderr)
    return mean_r2, cnt, n


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    label = nice_bin_label(args.bin_bp)
    default_outpath = os.path.join(args.outdir, f"{args.chrom}.{args.cache_type}.bin{label}.ld_cache.pkl")
    outpath = args.outfile if args.outfile else default_outpath
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    if args.use_previous_binary:
        if os.path.exists(outpath):
            print(f"[cache] Reusing: {outpath}", file=sys.stderr)
            return
        print(f"[cache] ERROR: --use-previous-binary but missing: {outpath}", file=sys.stderr)
        sys.exit(1)

    if os.path.exists(outpath) and not args.recreate_binary_matrix:
        print(f"[cache] Exists, reusing: {outpath}", file=sys.stderr)
        return

    p1, p2, r2 = parse_raw(args.infile)
    mat, cnt, nb = build_matrix(p1, p2, r2, args.chr_len, args.bin_bp)

    cache = dict(
        chrom=args.chrom,
        chr_len=args.chr_len,
        bin_bp=args.bin_bp,
        cache_type=args.cache_type,
        n_bins=nb,
        r2_matrix=mat,
        count_matrix=cnt,
        source_file=os.path.abspath(args.infile),
        total_pairs=len(p1),
        retained_pairs=int(np.sum(cnt) // 2),
        created=time.strftime("%Y-%m-%d %H:%M:%S"),
        outfile=os.path.abspath(outpath),
    )

    with open(outpath, "wb") as f:
        pickle.dump(cache, f, protocol=pickle.HIGHEST_PROTOCOL)

    sz = os.path.getsize(outpath) / 1e6
    print(f"[cache] Saved: {outpath} ({sz:.1f} MB)", file=sys.stderr)


if __name__ == "__main__":
    main()
