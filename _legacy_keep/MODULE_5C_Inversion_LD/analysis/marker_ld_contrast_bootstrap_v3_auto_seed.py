#!/usr/bin/env python3
"""
marker_ld_contrast_bootstrap_v3_auto_seed.py

Automatic seed-marker version with Ind0/Ind1/... -> real sample-name remapping.

Inputs:
  --candidate-tsv              Candidate master table with candidate_id, chrom, start_bp, end_bp
  --candidate-id               Candidate to analyze
  --sites-tsv                  Chromosome sites TSV/GZ with at least: marker, chrom, pos
  --dosage-tsv                 Chromosome dosage TSV/GZ with first column marker and sample columns
  --sample-order-file          Optional TXT file mapping Ind0..IndN dosage columns to real sample IDs
  --full-a-samples             TXT file, one sample per line
  --full-b-samples             TXT file, one sample per line
  --half-samples               TXT file, one sample per line
  --candidate-pairs            Candidate-region pairs_r2.tsv
  --chrom-pairs                Chromosome-wide/background pairs_r2.tsv
  --outprefix                  Output prefix

Optional:
  --neutral-mode flank|outside_all
  --candidate-padding-bp INT
  --min-neutral-markers INT
  --n-bootstrap INT
  --min-pairs INT
  --min-sep FLOAT
  --require-half-between
  --top-n-seed INT
  --random-seed INT

Outputs:
  <outprefix>.seed_markers.tsv
  <outprefix>.seed_markers.txt
  <outprefix>.neutral_markers.txt
  <outprefix>.marker_ld_contrast.tsv
  <outprefix>.marker_ld_contrast.top.tsv
  <outprefix>.summary.txt
"""

import argparse
import gzip
import random
import re
from collections import defaultdict
from statistics import mean
import math


def openth(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate-tsv", required=True)
    ap.add_argument("--candidate-id", required=True)

    ap.add_argument("--sites-tsv", required=True)
    ap.add_argument("--dosage-tsv", required=True)
    ap.add_argument("--sample-order-file", default=None)

    ap.add_argument("--full-a-samples", required=True)
    ap.add_argument("--full-b-samples", required=True)
    ap.add_argument("--half-samples", required=True)

    ap.add_argument("--candidate-pairs", required=True)
    ap.add_argument("--chrom-pairs", required=True)

    ap.add_argument("--outprefix", required=True)

    ap.add_argument("--neutral-mode", choices=["flank", "outside_all"], default="flank")
    ap.add_argument("--candidate-padding-bp", type=int, default=0)
    ap.add_argument("--min-neutral-markers", type=int, default=100)

    ap.add_argument("--n-bootstrap", type=int, default=200)
    ap.add_argument("--min-pairs", type=int, default=5)

    ap.add_argument("--min-sep", type=float, default=0.6)
    ap.add_argument("--require-half-between", action="store_true")
    ap.add_argument("--top-n-seed", type=int, default=300)

    ap.add_argument("--random-seed", type=int, default=1)
    return ap.parse_args()


def read_candidate(candidate_tsv, candidate_id):
    with openth(candidate_tsv, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {x: i for i, x in enumerate(header)}
        for req in ["candidate_id", "chrom", "start_bp", "end_bp"]:
            if req not in idx:
                raise ValueError(f"Missing required column in candidate table: {req}")

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if str(parts[idx["candidate_id"]]) == str(candidate_id):
                return {
                    "candidate_id": str(parts[idx["candidate_id"]]),
                    "chrom": parts[idx["chrom"]],
                    "start_bp": int(parts[idx["start_bp"]]),
                    "end_bp": int(parts[idx["end_bp"]]),
                }

    raise ValueError(f"Candidate {candidate_id} not found in {candidate_tsv}")


def normalize_marker_id(marker):
    if ":" in marker:
        return marker
    if "_" in marker:
        a, b = marker.rsplit("_", 1)
        if b.isdigit():
            return f"{a}:{b}"
    return marker

def read_sample_list(path):
    vals = []
    with open(path) as f:
        for line in f:
            x = line.strip()
            if x:
                vals.append(x)
    return vals


def read_sample_order(path):
    vals = []
    with open(path) as f:
        for line in f:
            x = line.strip()
            if x:
                vals.append(x)
    return vals


def parse_marker_pos(marker, default_chrom=None):
    if ":" in marker:
        a, b = marker.rsplit(":", 1)
        try:
            return a, int(b)
        except ValueError:
            return None, None
    if "_" in marker:
        a, b = marker.rsplit("_", 1)
        try:
            return a, int(b)
        except ValueError:
            return None, None
    try:
        return default_chrom, int(marker)
    except ValueError:
        return None, None

def read_sites(path):
    with openth(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {x: i for i, x in enumerate(header)}
        for req in ["marker", "chrom", "pos"]:
            if req not in idx:
                raise ValueError(f"Missing required column in sites file: {req}")

        rows = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rows.append({
                "marker": normalize_marker_id(parts[idx["marker"]]),
                #"marker": parts[idx["marker"]],
                "chrom": parts[idx["chrom"]],
                "pos": int(parts[idx["pos"]]),
            })
    return rows


def read_dosage(path, sample_order_file=None):
    """
    Returns:
      samples: list of sample column names
      dosage: dict marker -> dict sample -> float
    """
    dosage = {}
    with openth(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[0] != "marker":
            raise ValueError("First column of dosage file must be 'marker'")
        samples = header[1:]

        # Remap Ind0/Ind1/... to real sample IDs if needed
        if all(re.match(r"^Ind[0-9]+$", s) for s in samples):
            if sample_order_file is None:
                raise ValueError(
                    "Dosage columns are Ind-style but --sample-order-file was not provided"
                )
            sample_order = read_sample_order(sample_order_file)
            if len(sample_order) != len(samples):
                raise ValueError(
                    f"sample_order_file length ({len(sample_order)}) does not match dosage sample count ({len(samples)})"
                )
            samples = sample_order

        for line in f:
            parts = line.rstrip("\n").split("\t")
            marker = normalize_marker_id(parts[0])
            #marker = parts[0]
            row = {}
            for s, v in zip(samples, parts[1:]):
                if v in ("NA", ".", ""):
                    row[s] = None
                else:
                    try:
                        row[s] = float(v)
                    except ValueError:
                        row[s] = None
            dosage[marker] = row

    return samples, dosage


def safe_mean(vals):
    vals = [x for x in vals if x is not None and math.isfinite(x)]
    return None if not vals else mean(vals)


def half_between(muA, muH, muB):
    if muA is None or muH is None or muB is None:
        return False
    lo = min(muA, muB)
    hi = max(muA, muB)
    return lo <= muH <= hi


def compute_seed_markers(candidate, sites_rows, dosage, full_a, full_b, half, min_sep=0.6, require_half_between=False, top_n_seed=300):
    chrom = candidate["chrom"]
    start_bp = candidate["start_bp"]
    end_bp = candidate["end_bp"]

    rows = []
    for rec in sites_rows:
        if rec["chrom"] != chrom:
            continue
        if not (start_bp <= rec["pos"] <= end_bp):
            continue

        marker = rec["marker"]
        if marker not in dosage:
            continue

        row = dosage[marker]
        muA = safe_mean([row.get(s) for s in full_a])
        muB = safe_mean([row.get(s) for s in full_b])
        muH = safe_mean([row.get(s) for s in half])

        if muA is None or muB is None:
            continue

        sepAB = abs(muA - muB)
        hb = half_between(muA, muH, muB)

        if muA <= muB:
            orientation = "A_low_B_high"
        else:
            orientation = "A_high_B_low"

        score = sepAB
        if hb:
            score += 0.25

        selected = sepAB >= min_sep
        if require_half_between and not hb:
            selected = False

        rows.append({
            "marker": marker,
            "chrom": rec["chrom"],
            "pos": rec["pos"],
            "muA": muA,
            "muH": muH,
            "muB": muB,
            "sepAB": sepAB,
            "half_between": hb,
            "orientation": orientation,
            "seed_score": score,
            "selected_initial": selected,
        })

    rows.sort(key=lambda x: (x["selected_initial"], x["seed_score"]), reverse=True)

    selected = [r for r in rows if r["selected_initial"]]
    if top_n_seed > 0:
        selected = selected[:top_n_seed]

    selected_markers = {r["marker"] for r in selected}

    for r in rows:
        r["selected_as_seed"] = r["marker"] in selected_markers

    return rows, [r["marker"] for r in selected]


def read_pairs(path):
    adj = defaultdict(dict)
    markers = set()

    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            a = normalize_marker_id(parts[0])
            b = normalize_marker_id(parts[1])
           # a = parts[0]
           # b = parts[1]
            try:
                r2 = float(parts[3])
            except ValueError:
                continue
            adj[a][b] = r2
            adj[b][a] = r2
            markers.add(a)
            markers.add(b)

    return adj, markers


def choose_neutral_markers(all_chrom_markers, chrom, start_bp, end_bp, mode="flank", padding_bp=0):
    cand_start = start_bp - padding_bp
    cand_end = end_bp + padding_bp
    cand_len = end_bp - start_bp + 1

    left, right, outside_all = [], [], []

    for m in all_chrom_markers:
        mch, pos = parse_marker_pos(m, default_chrom=chrom)
        if mch != chrom or pos is None:
            continue
        if pos < cand_start:
            left.append((m, pos))
            outside_all.append((m, pos))
        elif pos > cand_end:
            right.append((m, pos))
            outside_all.append((m, pos))

    if mode == "outside_all":
        return [m for m, _ in sorted(outside_all, key=lambda x: x[1])]

    left_flank = []
    right_flank = []

    for m, pos in sorted(left, key=lambda x: x[1], reverse=True):
        if pos >= cand_start - cand_len:
            left_flank.append((m, pos))

    for m, pos in sorted(right, key=lambda x: x[1]):
        if pos <= cand_end + cand_len:
            right_flank.append((m, pos))

    flank = sorted(left_flank + right_flank, key=lambda x: x[1])
    return [m for m, _ in flank]


def mean_r2_to_set(marker, refset, adj):
    vals = []
    nbrs = adj.get(marker, {})
    for x in refset:
        if x == marker:
            continue
        if x in nbrs:
            vals.append(nbrs[x])
    if not vals:
        return None, 0
    return mean(vals), len(vals)


def bootstrap_score(marker, seed_markers, neutral_markers, cand_adj, chrom_adj, n_boot=200, min_pairs=5):
    boot_scores = []

    for _ in range(n_boot):
        seed_draw = [random.choice(seed_markers) for _ in range(len(seed_markers))]
        neutral_draw = [random.choice(neutral_markers) for _ in range(len(neutral_markers))]

        s_mean, s_n = mean_r2_to_set(marker, seed_draw, cand_adj)
        n_mean, n_n = mean_r2_to_set(marker, neutral_draw, chrom_adj)

        if s_mean is None or n_mean is None:
            continue
        if s_n < min_pairs or n_n < min_pairs:
            continue

        boot_scores.append(s_mean - n_mean)

    if not boot_scores:
        return None

    boot_scores.sort()
    n = len(boot_scores)
    lo = boot_scores[max(0, int(0.025 * n) - 1)]
    hi = boot_scores[min(n - 1, int(0.975 * n))]
    pos_frac = sum(x > 0 for x in boot_scores) / n

    return {
        "score_mean": mean(boot_scores),
        "score_min": boot_scores[0],
        "score_max": boot_scores[-1],
        "score_ci_low": lo,
        "score_ci_high": hi,
        "score_positive_fraction": pos_frac,
        "n_boot_used": n,
    }


def main():
    args = parse_args()
    random.seed(args.random_seed)

    cand = read_candidate(args.candidate_tsv, args.candidate_id)

    full_a = read_sample_list(args.full_a_samples)
    full_b = read_sample_list(args.full_b_samples)
    half = read_sample_list(args.half_samples)

    sites_rows = read_sites(args.sites_tsv)
    dosage_samples, dosage = read_dosage(args.dosage_tsv, sample_order_file=args.sample_order_file)

    # Keep only groups that exist in dosage after remapping
    full_a = [s for s in full_a if s in dosage_samples]
    full_b = [s for s in full_b if s in dosage_samples]
    half = [s for s in half if s in dosage_samples]

    if len(full_a) == 0 or len(full_b) == 0:
        raise RuntimeError(
            f"After sample-name reconciliation, group sizes are FULL_A={len(full_a)}, FULL_B={len(full_b)}, HALF={len(half)}"
        )

    seed_rows, seed_markers = compute_seed_markers(
        candidate=cand,
        sites_rows=sites_rows,
        dosage=dosage,
        full_a=full_a,
        full_b=full_b,
        half=half,
        min_sep=args.min_sep,
        require_half_between=args.require_half_between,
        top_n_seed=args.top_n_seed,
    )

    if not seed_markers:
        raise RuntimeError("No seed markers selected. Relax thresholds.")

    cand_adj, cand_markers = read_pairs(args.candidate_pairs)
    chrom_adj, chrom_markers = read_pairs(args.chrom_pairs)

    neutral_markers = choose_neutral_markers(
        all_chrom_markers=chrom_markers,
        chrom=cand["chrom"],
        start_bp=cand["start_bp"],
        end_bp=cand["end_bp"],
        mode=args.neutral_mode,
        padding_bp=args.candidate_padding_bp,
    )
    neutral_markers = [m for m in neutral_markers if m not in cand_markers]

    if len(neutral_markers) < args.min_neutral_markers:
        raise RuntimeError(
            f"Too few neutral markers ({len(neutral_markers)}). "
            f"Try --neutral-mode outside_all or lower --min-neutral-markers."
        )

    out_seed_tsv = args.outprefix + ".seed_markers.tsv"
    out_seed_txt = args.outprefix + ".seed_markers.txt"
    out_neutral = args.outprefix + ".neutral_markers.txt"
    out_tsv = args.outprefix + ".marker_ld_contrast.tsv"
    out_top = args.outprefix + ".marker_ld_contrast.top.tsv"
    out_summary = args.outprefix + ".summary.txt"

    with open(out_seed_tsv, "w") as f:
        f.write("\t".join([
            "marker", "chrom", "pos", "muA", "muH", "muB",
            "sepAB", "half_between", "orientation",
            "seed_score", "selected_initial", "selected_as_seed"
        ]) + "\n")
        for r in seed_rows:
            f.write("\t".join([
                r["marker"],
                r["chrom"],
                str(r["pos"]),
                "NA" if r["muA"] is None else f"{r['muA']:.6f}",
                "NA" if r["muH"] is None else f"{r['muH']:.6f}",
                "NA" if r["muB"] is None else f"{r['muB']:.6f}",
                f"{r['sepAB']:.6f}",
                str(r["half_between"]),
                r["orientation"],
                f"{r['seed_score']:.6f}",
                str(r["selected_initial"]),
                str(r["selected_as_seed"]),
            ]) + "\n")

    with open(out_seed_txt, "w") as f:
        for m in seed_markers:
            f.write(m + "\n")

    with open(out_neutral, "w") as f:
        for m in neutral_markers:
            f.write(m + "\n")

    rows = []
    for i, marker in enumerate(sorted(cand_markers), start=1):
        res = bootstrap_score(
            marker=marker,
            seed_markers=seed_markers,
            neutral_markers=neutral_markers,
            cand_adj=cand_adj,
            chrom_adj=chrom_adj,
            n_boot=args.n_bootstrap,
            min_pairs=args.min_pairs,
        )
        if res is None:
            continue
        rows.append({"marker": marker, **res})
        if i % 1000 == 0:
            print(f"[INFO] scored {i} candidate markers")

    rows.sort(key=lambda x: x["score_mean"], reverse=True)

    with open(out_tsv, "w") as f:
        f.write("\t".join([
            "marker", "score_mean", "score_min", "score_max",
            "score_ci_low", "score_ci_high",
            "score_positive_fraction", "n_boot_used"
        ]) + "\n")
        for r in rows:
            f.write("\t".join([
                r["marker"],
                f"{r['score_mean']:.6f}",
                f"{r['score_min']:.6f}",
                f"{r['score_max']:.6f}",
                f"{r['score_ci_low']:.6f}",
                f"{r['score_ci_high']:.6f}",
                f"{r['score_positive_fraction']:.6f}",
                str(r["n_boot_used"]),
            ]) + "\n")

    top_n = min(200, len(rows))
    with open(out_top, "w") as f:
        f.write("\t".join([
            "marker", "score_mean", "score_ci_low",
            "score_ci_high", "score_positive_fraction", "n_boot_used"
        ]) + "\n")
        for r in rows[:top_n]:
            f.write("\t".join([
                r["marker"],
                f"{r['score_mean']:.6f}",
                f"{r['score_ci_low']:.6f}",
                f"{r['score_ci_high']:.6f}",
                f"{r['score_positive_fraction']:.6f}",
                str(r["n_boot_used"]),
            ]) + "\n")

    with open(out_summary, "w") as f:
        f.write(f"candidate_id\t{cand['candidate_id']}\n")
        f.write(f"chrom\t{cand['chrom']}\n")
        f.write(f"start_bp\t{cand['start_bp']}\n")
        f.write(f"end_bp\t{cand['end_bp']}\n")
        f.write(f"full_a_samples_used\t{len(full_a)}\n")
        f.write(f"full_b_samples_used\t{len(full_b)}\n")
        f.write(f"half_samples_used\t{len(half)}\n")
        f.write(f"seed_markers\t{len(seed_markers)}\n")
        f.write(f"candidate_markers_in_pairs\t{len(cand_markers)}\n")
        f.write(f"neutral_markers\t{len(neutral_markers)}\n")
        f.write(f"neutral_mode\t{args.neutral_mode}\n")
        f.write(f"n_bootstrap\t{args.n_bootstrap}\n")
        f.write(f"min_pairs\t{args.min_pairs}\n")
        f.write(f"min_sep\t{args.min_sep}\n")
        f.write(f"require_half_between\t{args.require_half_between}\n")
        f.write(f"sample_order_file\t{args.sample_order_file}\n")

    print(f"[DONE] wrote {out_seed_tsv}")
    print(f"[DONE] wrote {out_seed_txt}")
    print(f"[DONE] wrote {out_neutral}")
    print(f"[DONE] wrote {out_tsv}")
    print(f"[DONE] wrote {out_top}")
    print(f"[DONE] wrote {out_summary}")


if __name__ == "__main__":
    main()
