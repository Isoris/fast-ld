#!/usr/bin/env python3
# plot_ld_contrast_score_vs_pos.py

import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True, help="marker_ld_contrast.tsv")
    ap.add_argument("--outfile", required=True, help="output PNG")
    ap.add_argument("--min-positive-frac", type=float, default=None,
                    help="Optional filter on score_positive_fraction")
    args = ap.parse_args()

    df = pd.read_csv(args.infile, sep="\t")
    if df.empty:
        raise SystemExit("Input table is empty")

    df["pos"] = df["marker"].astype(str).str.rsplit(":", n=1).str[-1].astype(int)

    if args.min_positive_frac is not None and "score_positive_fraction" in df.columns:
        df = df[df["score_positive_fraction"] >= args.min_positive_frac].copy()

    if df.empty:
        raise SystemExit("No rows left after filtering")

    plt.figure(figsize=(12, 4))
    plt.scatter(df["pos"], df["score_mean"], s=8)
    plt.axhline(0, linewidth=1)
    plt.xlabel("Position (bp)")
    plt.ylabel("LD contrast score")
    plt.tight_layout()
    plt.savefig(args.outfile, dpi=350)


if __name__ == "__main__":
    main()
