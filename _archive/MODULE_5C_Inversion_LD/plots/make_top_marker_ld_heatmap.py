#!/usr/bin/env python3
# make_top_marker_ld_heatmap.py

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-table", required=True, help="marker_ld_contrast.top.tsv")
    ap.add_argument("--pairs", required=True, help="candidate invzoom pairs_r2.tsv")
    ap.add_argument("--outfile", required=True, help="output PNG")
    ap.add_argument("--top-n", type=int, default=100, help="number of top markers to keep")
    ap.add_argument("--title", default="Focused LD heatmap")
    args = ap.parse_args()

    top = pd.read_csv(args.top_table, sep="\t")
    if top.empty:
        raise SystemExit("Top-marker table is empty")

    top_markers = set(top["marker"].head(args.top_n).astype(str))

    df = pd.read_csv(args.pairs, sep=r"\s+", header=None, names=["m1", "m2", "dist", "r2"])
    df = df[df["m1"].isin(top_markers) & df["m2"].isin(top_markers)].copy()

    if df.empty:
        raise SystemExit("No LD pairs found among selected top markers")

    markers = sorted(
        set(df["m1"]) | set(df["m2"]),
        key=lambda x: int(str(x).rsplit(":", 1)[1])
    )
    idx = {m: i for i, m in enumerate(markers)}

    mat = np.full((len(markers), len(markers)), np.nan, dtype=float)
    for _, r in df.iterrows():
        i = idx[r["m1"]]
        j = idx[r["m2"]]
        mat[i, j] = r["r2"]
        mat[j, i] = r["r2"]

    np.fill_diagonal(mat, 1.0)

    plt.figure(figsize=(8, 7))
    plt.imshow(mat, origin="lower", aspect="auto")
    plt.colorbar(label="r²")
    plt.title(args.title)
    plt.tight_layout()
    plt.savefig(args.outfile, dpi=350)


if __name__ == "__main__":
    main()
