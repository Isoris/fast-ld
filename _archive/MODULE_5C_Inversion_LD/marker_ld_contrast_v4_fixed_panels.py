#!/usr/bin/env python3
"""
marker_ld_contrast_v4_fixed_panels.py

Upgrade over v3:
- automatic seed markers from dosage groups
- Ind0/Ind1/... remapping to real sample IDs
- FIXED bootstrap seed panels and FIXED bootstrap neutral panels shared across all markers
- per-marker contrast matrix
- marker profile similarity matrix
- coherent cluster extraction

Outputs:
  <outprefix>.seed_markers.tsv
  <outprefix>.seed_markers.txt
  <outprefix>.neutral_markers.txt
  <outprefix>.bootstrap_panels.summary.tsv
  <outprefix>.marker_ld_contrast.tsv
  <outprefix>.marker_ld_contrast.top.tsv
  <outprefix>.marker_bootstrap_matrix.tsv.gz
  <outprefix>.marker_profile_similarity.tsv.gz
  <outprefix>.marker_profile_cluster_assignments.tsv
  <outprefix>.largest_profile_cluster.markers.txt
  <outprefix>.summary.txt
"""

import argparse
import gzip
import math
import random
import re
from collections import Counter, defaultdict
from statistics import mean


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

    ap.add_argument("--neutral-mode", choices=["flank", "outside_all"], default="outside_all")
    ap.add_argument("--candidate-padding-bp", type=int, default=0)
    ap.add_argument("--min-neutral-markers", type=int, default=100)

    ap.add_argument("--n-bootstrap", type=int, default=200)
    ap.add_argument("--seed-panel-size", type=int, default=200)
    ap.add_argument("--neutral-panel-size", type=int, default=200)
    ap.add_argument("--min-pairs", type=int, default=5)

    ap.add_argument("--min-sep", type=float, default=0.6)
    ap.add_argument("--require-half-between", action="store_true")
    ap.add_argument("--top-n-seed", type=int, default=300)

    ap.add_argument("--profile-min-corr", type=float, default=0.7)
    ap.add_argument("--profile-min-shared-panels", type=int, default=20)

    ap.add_argument("--random-seed", type=int, default=1)
    return ap.parse_args()


def normalize_marker_id(marker):
    if ":" in marker:
        return marker
    if "_" in marker:
        a, b = marker.rsplit("_", 1)
        if b.isdigit():
            return f"{a}:{b}"
    return marker


def parse_marker_pos(marker, default_chrom=None):
    marker = normalize_marker_id(marker)
    if ":" in marker:
        a, b = marker.rsplit(":", 1)
        try:
            return a, int(b)
        except ValueError:
            return None, None
    try:
        return default_chrom, int(marker)
    except ValueError:
        return None, None


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
                "chrom": parts[idx["chrom"]],
                "pos": int(parts[idx["pos"]]),
            })
    return rows


def read_dosage(path, sample_order_file=None):
    dosage = {}
    with openth(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[0] != "marker":
            raise ValueError("First column of dosage file must be 'marker'")
        samples = header[1:]

        if all(re.match(r"^Ind[0-9]+$", s) for s in samples):
            if sample_order_file is None:
                raise ValueError("Dosage columns are Ind-style but --sample-order-file was not provided")
            sample_order = read_sample_order(sample_order_file)
            if len(sample_order) != len(samples):
                raise ValueError(
                    f"sample_order_file length ({len(sample_order)}) does not match dosage sample count ({len(samples)})"
                )
            samples = sample_order

        for line in f:
            parts = line.rstrip("\n").split("\t")
            marker = normalize_marker_id(parts[0])
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
            try:
                r2 = float(parts[3])
            except ValueError:
                continue
            adj[a][b] = r2
            adj[b][a] = r2
            markers.add(a)
            markers.add(b)

    return adj, markers


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

        orientation = "A_low_B_high" if muA <= muB else "A_high_B_low"

        score = sepAB + (0.25 if hb else 0.0)

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


def choose_neutral_markers(all_chrom_markers, chrom, start_bp, end_bp, mode="outside_all", padding_bp=0):
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


def build_fixed_panels(marker_pool, n_boot, panel_size, rng):
    marker_pool = list(marker_pool)
    if len(marker_pool) == 0:
        return []
    panels = []
    for _ in range(n_boot):
        panels.append([rng.choice(marker_pool) for _ in range(panel_size)])
    return panels


def score_all_markers_against_fixed_panels(markers_to_score, seed_panels, neutral_panels, cand_adj, chrom_adj, min_pairs):
    results = []
    matrix_rows = []

    for idx_m, marker in enumerate(markers_to_score, start=1):
        contrasts = []
        used_seed_counts = []
        used_neutral_counts = []

        for bidx, (seed_panel, neutral_panel) in enumerate(zip(seed_panels, neutral_panels), start=1):
            s_mean, s_n = mean_r2_to_set(marker, seed_panel, cand_adj)
            n_mean, n_n = mean_r2_to_set(marker, neutral_panel, chrom_adj)

            if s_mean is None or n_mean is None:
                contrast = None
            elif s_n < min_pairs or n_n < min_pairs:
                contrast = None
            else:
                contrast = s_mean - n_mean

            matrix_rows.append({
                "marker": marker,
                "bootstrap_id": bidx,
                "seed_n_pairs": s_n,
                "neutral_n_pairs": n_n,
                "contrast": contrast,
            })

            if contrast is not None:
                contrasts.append(contrast)
                used_seed_counts.append(s_n)
                used_neutral_counts.append(n_n)

        if contrasts:
            contrasts_sorted = sorted(contrasts)
            n = len(contrasts_sorted)
            lo = contrasts_sorted[max(0, int(0.025 * n) - 1)]
            hi = contrasts_sorted[min(n - 1, int(0.975 * n))]
            pos_frac = sum(x > 0 for x in contrasts_sorted) / n

            results.append({
                "marker": marker,
                "score_mean": mean(contrasts_sorted),
                "score_min": contrasts_sorted[0],
                "score_max": contrasts_sorted[-1],
                "score_ci_low": lo,
                "score_ci_high": hi,
                "score_positive_fraction": pos_frac,
                "n_boot_used": n,
                "mean_seed_n_pairs": mean(used_seed_counts),
                "mean_neutral_n_pairs": mean(used_neutral_counts),
            })

        if idx_m % 1000 == 0:
            print(f"[INFO] scored {idx_m} markers")

    results.sort(key=lambda x: x["score_mean"], reverse=True)
    return results, matrix_rows


def pearson_corr(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(pairs) < 3:
        return None, len(pairs)
    xs = [x for x, _ in pairs]
    ys = [y for _, y in pairs]
    mx = mean(xs)
    my = mean(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx <= 0 or vy <= 0:
        return None, len(pairs)
    cov = sum((x - mx) * (y - my) for x, y in pairs)
    return cov / math.sqrt(vx * vy), len(pairs)


def build_profile_similarity(marker_matrix, min_shared_panels=20):
    by_marker = defaultdict(list)
    for r in marker_matrix:
        by_marker[r["marker"]].append((r["bootstrap_id"], r["contrast"]))

    vectors = {}
    for marker, vals in by_marker.items():
        vals = sorted(vals, key=lambda x: x[0])
        vectors[marker] = [v for _, v in vals]

    markers = sorted(vectors.keys())
    edges = []
    sim_rows = []

    for i in range(len(markers)):
        m1 = markers[i]
        for j in range(i + 1, len(markers)):
            m2 = markers[j]
            corr, n_shared = pearson_corr(vectors[m1], vectors[m2])
            sim_rows.append({
                "marker1": m1,
                "marker2": m2,
                "profile_corr": corr,
                "n_shared_panels": n_shared,
            })
            if corr is not None and n_shared >= min_shared_panels:
                edges.append((m1, m2, corr))

    return sim_rows, edges


def extract_largest_cluster(edges, min_corr=0.7):
    graph = defaultdict(set)
    for a, b, corr in edges:
        if corr is not None and corr >= min_corr:
            graph[a].add(b)
            graph[b].add(a)

    seen = set()
    components = []
    for node in graph:
        if node in seen:
            continue
        stack = [node]
        comp = []
        seen.add(node)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in graph[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        components.append(sorted(comp))

    if not components:
        return []

    components.sort(key=len, reverse=True)
    return components[0]


def write_tsv(path, header, rows):
    with openth(path, "wt") if path.endswith(".gz") else open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            vals = []
            for h in header:
                v = r.get(h, "")
                if v is None:
                    vals.append("NA")
                else:
                    vals.append(str(v))
            f.write("\t".join(vals) + "\n")


def main():
    args = parse_args()
    rng = random.Random(args.random_seed)

    cand = read_candidate(args.candidate_tsv, args.candidate_id)

    full_a = read_sample_list(args.full_a_samples)
    full_b = read_sample_list(args.full_b_samples)
    half = read_sample_list(args.half_samples)

    sites_rows = read_sites(args.sites_tsv)
    dosage_samples, dosage = read_dosage(args.dosage_tsv, sample_order_file=args.sample_order_file)

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

    seed_panel_size = min(args.seed_panel_size, len(seed_markers))
    neutral_panel_size = min(args.neutral_panel_size, len(neutral_markers))

    seed_panels = build_fixed_panels(seed_markers, args.n_bootstrap, seed_panel_size, rng)
    neutral_panels = build_fixed_panels(neutral_markers, args.n_bootstrap, neutral_panel_size, rng)

    results, matrix_rows = score_all_markers_against_fixed_panels(
        markers_to_score=sorted(cand_markers),
        seed_panels=seed_panels,
        neutral_panels=neutral_panels,
        cand_adj=cand_adj,
        chrom_adj=chrom_adj,
        min_pairs=args.min_pairs,
    )

    sim_rows, edges = build_profile_similarity(
        marker_matrix=matrix_rows,
        min_shared_panels=args.profile_min_shared_panels,
    )
    largest_cluster = extract_largest_cluster(edges, min_corr=args.profile_min_corr)

    out_seed_tsv = args.outprefix + ".seed_markers.tsv"
    out_seed_txt = args.outprefix + ".seed_markers.txt"
    out_neutral = args.outprefix + ".neutral_markers.txt"
    out_panel_summary = args.outprefix + ".bootstrap_panels.summary.tsv"
    out_tsv = args.outprefix + ".marker_ld_contrast.tsv"
    out_top = args.outprefix + ".marker_ld_contrast.top.tsv"
    out_matrix = args.outprefix + ".marker_bootstrap_matrix.tsv.gz"
    out_sim = args.outprefix + ".marker_profile_similarity.tsv.gz"
    out_cluster_tsv = args.outprefix + ".marker_profile_cluster_assignments.tsv"
    out_cluster_txt = args.outprefix + ".largest_profile_cluster.markers.txt"
    out_summary = args.outprefix + ".summary.txt"

    write_tsv(
        out_seed_tsv,
        ["marker", "chrom", "pos", "muA", "muH", "muB", "sepAB", "half_between", "orientation", "seed_score", "selected_initial", "selected_as_seed"],
        seed_rows,
    )

    with open(out_seed_txt, "w") as f:
        for m in seed_markers:
            f.write(m + "\n")

    with open(out_neutral, "w") as f:
        for m in neutral_markers:
            f.write(m + "\n")

    panel_rows = []
    for i in range(args.n_bootstrap):
        panel_rows.append({
            "bootstrap_id": i + 1,
            "seed_panel_size": len(seed_panels[i]),
            "neutral_panel_size": len(neutral_panels[i]),
            "seed_unique_markers": len(set(seed_panels[i])),
            "neutral_unique_markers": len(set(neutral_panels[i])),
        })
    write_tsv(
        out_panel_summary,
        ["bootstrap_id", "seed_panel_size", "neutral_panel_size", "seed_unique_markers", "neutral_unique_markers"],
        panel_rows,
    )

    write_tsv(
        out_tsv,
        ["marker", "score_mean", "score_min", "score_max", "score_ci_low", "score_ci_high", "score_positive_fraction", "n_boot_used", "mean_seed_n_pairs", "mean_neutral_n_pairs"],
        results,
    )

    write_tsv(
        out_top,
        ["marker", "score_mean", "score_ci_low", "score_ci_high", "score_positive_fraction", "n_boot_used"],
        results[: min(200, len(results))],
    )

    write_tsv(
        out_matrix,
        ["marker", "bootstrap_id", "seed_n_pairs", "neutral_n_pairs", "contrast"],
        matrix_rows,
    )

    write_tsv(
        out_sim,
        ["marker1", "marker2", "profile_corr", "n_shared_panels"],
        sim_rows,
    )

    cluster_set = set(largest_cluster)
    cluster_rows = [{"marker": m, "in_largest_profile_cluster": m in cluster_set} for m in sorted({r["marker"] for r in results})]
    write_tsv(
        out_cluster_tsv,
        ["marker", "in_largest_profile_cluster"],
        cluster_rows,
    )

    with open(out_cluster_txt, "w") as f:
        for m in largest_cluster:
            f.write(m + "\n")

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
        f.write(f"seed_panel_size\t{seed_panel_size}\n")
        f.write(f"neutral_panel_size\t{neutral_panel_size}\n")
        f.write(f"min_pairs\t{args.min_pairs}\n")
        f.write(f"min_sep\t{args.min_sep}\n")
        f.write(f"require_half_between\t{args.require_half_between}\n")
        f.write(f"profile_min_corr\t{args.profile_min_corr}\n")
        f.write(f"profile_min_shared_panels\t{args.profile_min_shared_panels}\n")
        f.write(f"largest_profile_cluster_size\t{len(largest_cluster)}\n")
        f.write(f"sample_order_file\t{args.sample_order_file}\n")

    print(f"[DONE] wrote {out_seed_tsv}")
    print(f"[DONE] wrote {out_seed_txt}")
    print(f"[DONE] wrote {out_neutral}")
    print(f"[DONE] wrote {out_panel_summary}")
    print(f"[DONE] wrote {out_tsv}")
    print(f"[DONE] wrote {out_top}")
    print(f"[DONE] wrote {out_matrix}")
    print(f"[DONE] wrote {out_sim}")
    print(f"[DONE] wrote {out_cluster_tsv}")
    print(f"[DONE] wrote {out_cluster_txt}")
    print(f"[DONE] wrote {out_summary}")


if __name__ == "__main__":
    main()
