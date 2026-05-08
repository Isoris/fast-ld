#!/usr/bin/env python3
"""
run_ld_candidate_suite_v4.py

Three-product LD heatmap orchestrator:
  - Whole-chromosome: chrsparse (2000 markers, unlimited distance)
  - Candidate zoom:   invzoom (2000 markers in zoom region, unlimited distance)
  - Fallback zoom:    localdense (all markers, 500kb distance — thin diagonal only)

For each candidate produces:
  - whole_chr: triangle_down + full
  - inv_zoom:  triangle_down + full
  - neg_ctrl:  full
"""
import argparse
import os
import sys
import subprocess
import csv


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate-tsv", required=True)
    ap.add_argument("--chrsparse-dir", required=True)
    ap.add_argument("--invzoom-dir", required=True, help="Dir with <chrom>.cand<id>.invzoom.pairs_r2.tsv")
    ap.add_argument("--chr-len-tsv", required=True)
    ap.add_argument("--cache-dir", default="ld_cache_v4")
    ap.add_argument("--outdir", default="ld_plots_v4")
    ap.add_argument("--target-bins", type=int, default=250)
    ap.add_argument("--pad-frac", type=float, default=0.05)
    ap.add_argument("--transform", default="power")
    ap.add_argument("--power", type=float, default=0.5)
    ap.add_argument("--cmap", default="apex")
    ap.add_argument("--dpi", type=int, default=400)
    ap.add_argument("--script-dir", default=".")
    ap.add_argument("--candidate-id", default=None)
    return ap.parse_args()


def read_chr_lens(p):
    d = {}
    with open(p) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[r["chrom"]] = int(r["chr_len"])
    return d


def round_bin(raw):
    if raw < 5000:
        return 5000
    if raw < 25000:
        return max(5000, (raw // 5000) * 5000)
    if raw < 100000:
        return max(10000, (raw // 10000) * 10000)
    return max(25000, (raw // 25000) * 25000)


def bl(bp):
    return f"{bp//1000}kb" if bp < 1_000_000 else f"{bp//1_000_000}Mb"


def neg_ctrl(cl, inv_s, inv_e, zs, ze):
    zl = ze - zs + 1
    ic = (inv_s + inv_e) / 2
    cc = []
    if inv_s - 1 >= zl:
        cc += [(1, zl), (inv_s - zl, inv_s - 1)]
    if cl - inv_e >= zl:
        cc += [(inv_e + 1, inv_e + zl), (cl - zl + 1, cl)]
    if not cc:
        return None, None
    return max(cc, key=lambda c: abs((c[0] + c[1]) / 2 - ic))


def run_cmd(cmd):
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Command failed:\n" + " ".join(cmd))


def main():
    args = parse_args()
    os.makedirs(args.cache_dir, exist_ok=True)
    for sub in ["whole_chr", "inv_zoom", "neg_ctrl"]:
        os.makedirs(os.path.join(args.outdir, sub), exist_ok=True)

    chr_lens = read_chr_lens(args.chr_len_tsv)
    cands = []
    with open(args.candidate_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if args.candidate_id and str(r["candidate_id"]) != str(args.candidate_id):
                continue
            cands.append(r)

    if not cands:
        print("[LD-v4] No candidates", file=sys.stderr)
        return

    build = os.path.join(args.script_dir, "STEP_LD_01_build_binary_cache.py")
    plot = os.path.join(args.script_dir, "STEP_LD_02_plot_heatmaps.py")

    chr_cands = {}
    for c in cands:
        chr_cands.setdefault(c["chrom"], []).append(c)

    summary = []

    for chrom, cc in chr_cands.items():
        cl = chr_lens.get(chrom)
        if not cl:
            continue

        # ── Build chromosome-sparse cache ─────────────────────────
        sparse_file = os.path.join(args.chrsparse_dir, f"{chrom}.chrsparse.pairs_r2.tsv")
        bin_whole = round_bin(cl // args.target_bins)
        cache_whole = os.path.join(
            args.cache_dir,
            f"{chrom}.chrsparse.bin{bl(bin_whole)}.ld_cache.pkl"
        )

        if os.path.exists(sparse_file) and not os.path.exists(cache_whole):
            run_cmd([
                sys.executable, build,
                "--infile", sparse_file,
                "--chrom", chrom,
                "--chr-len", str(cl),
                "--bin-bp", str(bin_whole),
                "--cache-type", "chrsparse",
                "--outdir", args.cache_dir,
                "--outfile", cache_whole
            ])

        has_whole = os.path.exists(cache_whole)

        for c in cc:
            cid = c["candidate_id"]
            inv_s = int(c["start_bp"])
            inv_e = int(c["end_bp"])
            inv_len = inv_e - inv_s

            pad = max(int(cl * args.pad_frac), inv_len // 2)
            zs = max(1, inv_s - pad)
            ze = min(cl, inv_e + pad)
            zoom_span = ze - zs

            # ── Build invzoom cache ───────────────────────────────
            zoom_file = os.path.join(args.invzoom_dir, f"{chrom}.cand{cid}.invzoom.pairs_r2.tsv")
            bin_zoom = round_bin(zoom_span // args.target_bins)
            cache_zoom = os.path.join(
                args.cache_dir,
                f"{chrom}.cand{cid}.{zs}_{ze}.invzoom.bin{bl(bin_zoom)}.ld_cache.pkl"
            )

            if os.path.exists(zoom_file) and not os.path.exists(cache_zoom):
                run_cmd([
                    sys.executable, build,
                    "--infile", zoom_file,
                    "--chrom", chrom,
                    "--chr-len", str(cl),
                    "--bin-bp", str(bin_zoom),
                    "--cache-type", "invzoom",
                    "--outdir", args.cache_dir,
                    "--outfile", cache_zoom
                ])

            has_zoom = os.path.exists(cache_zoom)

            common = [
                "--transform", args.transform,
                "--power", str(args.power),
                "--cmap", args.cmap,
                "--na-color", "white",
                "--candidate-id", str(cid),
                "--dpi", str(args.dpi)
            ]

            # ── Whole chromosome plots ────────────────────────────
            if has_whole:
                for mode in ["triangle_down", "full"]:
                    pfx = os.path.join(args.outdir, "whole_chr", f"{chrom}.cand{cid}.wholechr")
                    run_cmd([
                        sys.executable, plot,
                        "--cache", cache_whole,
                        "--outprefix", pfx,
                        "--region-start", "1",
                        "--region-end", str(cl),
                        "--highlight-start", str(inv_s),
                        "--highlight-end", str(inv_e),
                        "--mode", mode
                    ] + common)

            # ── Zoom plots ────────────────────────────────────────
            if has_zoom:
                for mode in ["triangle_down", "full"]:
                    pfx = os.path.join(args.outdir, "inv_zoom", f"{chrom}.cand{cid}.invzoom")
                    run_cmd([
                        sys.executable, plot,
                        "--cache", cache_zoom,
                        "--outprefix", pfx,
                        "--region-start", str(zs),
                        "--region-end", str(ze),
                        "--highlight-start", str(inv_s),
                        "--highlight-end", str(inv_e),
                        "--mode", mode
                    ] + common)

                ns, ne = neg_ctrl(cl, inv_s, inv_e, zs, ze)
                if ns and ne:
                    pfx = os.path.join(args.outdir, "neg_ctrl", f"{chrom}.cand{cid}.negctrl")
                    run_cmd([
                        sys.executable, plot,
                        "--cache", cache_zoom,
                        "--outprefix", pfx,
                        "--region-start", str(ns),
                        "--region-end", str(ne),
                        "--mode", "full"
                    ] + common)
            else:
                print(f"[LD-v4] WARN: no invzoom for cand {cid}, skipping zoom plots", file=sys.stderr)

            summary.append(dict(
                candidate_id=cid,
                chrom=chrom,
                inv_s=inv_s,
                inv_e=inv_e,
                zoom_s=zs,
                zoom_e=ze,
                bin_whole=bin_whole,
                bin_zoom=bin_zoom,
                has_whole=has_whole,
                has_zoom=has_zoom
            ))

        print(f"[LD-v4] Done: {chrom} ({len(cc)} cands)", file=sys.stderr)

    if summary:
        sf = os.path.join(args.outdir, "ld_plot_summary.tsv")
        with open(sf, "w") as f:
            w = csv.DictWriter(f, fieldnames=summary[0].keys(), delimiter="\t")
            w.writeheader()
            w.writerows(summary)

    print("[LD-v4] All done.", file=sys.stderr)


if __name__ == "__main__":
    main()
