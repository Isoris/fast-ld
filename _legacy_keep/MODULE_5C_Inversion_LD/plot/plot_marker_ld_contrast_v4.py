#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contrast-tsv", required=True)
    ap.add_argument("--outprefix", required=True)
    ap.add_argument("--cluster-markers", default=None)
    args = ap.parse_args()

    df = pd.read_csv(args.contrast_tsv, sep="\t")
    if df.empty:
        raise SystemExit("contrast table is empty")

    df["pos"] = df["marker"].astype(str).str.rsplit(":", n=1).str[-1].astype(int)
    df["in_cluster"] = False

    if args.cluster_markers:
        cluster = set(x.strip() for x in open(args.cluster_markers) if x.strip())
        df["in_cluster"] = df["marker"].isin(cluster)

    plt.figure(figsize=(12, 4))
    plt.scatter(df.loc[~df["in_cluster"], "pos"], df.loc[~df["in_cluster"], "score_mean"], s=8, alpha=0.6, label="other")
    if df["in_cluster"].any():
        plt.scatter(df.loc[df["in_cluster"], "pos"], df.loc[df["in_cluster"], "score_mean"], s=14, alpha=0.9, label="largest profile cluster")
    plt.axhline(0, linewidth=1)
    plt.xlabel("Position (bp)")
    plt.ylabel("LD contrast score")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.outprefix + ".score_vs_pos.png", dpi=200)

    plt.figure(figsize=(12, 4))
    plt.scatter(df["pos"], df["score_positive_fraction"], s=8)
    plt.xlabel("Position (bp)")
    plt.ylabel("score_positive_fraction")
    plt.tight_layout()
    plt.savefig(args.outprefix + ".positive_fraction_vs_pos.png", dpi=200)


if __name__ == "__main__":
    main()
