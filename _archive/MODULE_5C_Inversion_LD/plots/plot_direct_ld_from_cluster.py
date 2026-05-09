#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--cluster-markers", required=True)
    ap.add_argument("--outfile", required=True)
    ap.add_argument("--title", default="Direct LD heatmap of largest profile cluster")
    args = ap.parse_args()

    cluster = set(x.strip() for x in open(args.cluster_markers) if x.strip())
    if len(cluster) == 0:
        raise SystemExit("cluster marker file is empty")

    df = pd.read_csv(args.pairs, sep=r"\s+", header=None, names=["m1", "m2", "dist", "r2"])
    df = df[df["m1"].isin(cluster) & df["m2"].isin(cluster)].copy()
    if df.empty:
        raise SystemExit("no LD pairs found among cluster markers")

    markers = sorted(set(df["m1"]) | set(df["m2"]), key=lambda x: int(str(x).rsplit(":", 1)[1]))
    idx = {m: i for i, m in enumerate(markers)}

    mat = np.full((len(markers), len(markers)), np.nan, dtype=float)
    for _, r in df.iterrows():
        i = idx[r["m1"]]
        j = idx[r["m2"]]
        mat[i, j] = r["r2"]
        mat[j, i] = r["r2"]

    np.fill_diagonal(mat, 1.0)

    plt.figure(figsize=(8, 7))
    plt.imshow(mat, origin="lower", aspect="auto", vmin=0, vmax=1)
    plt.colorbar(label="r²")
    plt.title(args.title)
    plt.tight_layout()
    plt.savefig(args.outfile, dpi=350)


if __name__ == "__main__":
    main()
