#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile-sim-tsv", required=True)
    ap.add_argument("--cluster-markers", default=None)
    ap.add_argument("--outfile", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.profile_sim_tsv, sep="\t")
    df = df.dropna(subset=["profile_corr"]).copy()
    if df.empty:
        raise SystemExit("profile similarity table is empty")

    markers = sorted(set(df["marker1"]) | set(df["marker2"]), key=lambda x: int(str(x).rsplit(":", 1)[1]))

    if args.cluster_markers:
        wanted = set(x.strip() for x in open(args.cluster_markers) if x.strip())
        markers = [m for m in markers if m in wanted]
        df = df[df["marker1"].isin(wanted) & df["marker2"].isin(wanted)].copy()

    if len(markers) == 0:
        raise SystemExit("no markers left after filtering")

    idx = {m: i for i, m in enumerate(markers)}
    mat = np.full((len(markers), len(markers)), np.nan, dtype=float)

    for _, r in df.iterrows():
        i = idx[r["marker1"]]
        j = idx[r["marker2"]]
        mat[i, j] = r["profile_corr"]
        mat[j, i] = r["profile_corr"]

    np.fill_diagonal(mat, 1.0)

    plt.figure(figsize=(8, 7))
    plt.imshow(mat, origin="lower", aspect="auto", vmin=-1, vmax=1)
    plt.colorbar(label="profile correlation")
    plt.title("Marker profile-similarity heatmap")
    plt.tight_layout()
    plt.savefig(args.outfile, dpi=350)


if __name__ == "__main__":
    main()
