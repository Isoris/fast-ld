#!/usr/bin/env python3
"""
STEP_LD_02_plot_heatmaps.py  (v4)

LD heatmap plotter from binary cache.

Modes:
  full             — square symmetric heatmap
  upper_half       — upper triangle of square (masked lower)
  triangle_down    — rotated 45° pointing downward (standard Hi-C style)
  triangle_up      — rotated 45° pointing upward

Fixes in v4:
  - White background for NA/empty (not grey)
  - Both triangle orientations
  - Title spacing fix (no overlap with highlight bar)
  - Viridis palette option
  - Robust vmax with quantile mode
  - Proper colorbar scaling (vmax ≤ 1.0 for sqrt(r²))
"""
import argparse, os, sys, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Polygon, Rectangle
from matplotlib.collections import PatchCollection

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--outprefix", default="ld_heatmap")
    ap.add_argument("--region-start", type=int, default=None)
    ap.add_argument("--region-end", type=int, default=None)
    ap.add_argument("--highlight-start", type=int, default=None)
    ap.add_argument("--highlight-end", type=int, default=None)
    ap.add_argument("--mode", default="triangle_down",
                    choices=["full","upper_half","triangle_down","triangle_up"])
    ap.add_argument("--transform", default="power", choices=["raw","power","log1p"])
    ap.add_argument("--power", type=float, default=0.5)
    ap.add_argument("--log1p-k", type=float, default=10.0)
    ap.add_argument("--vmax-coef", type=float, default=2.0)
    ap.add_argument("--vmax-quantile", type=float, default=0.995,
                    help="Use this quantile of non-zero values for vmax")
    ap.add_argument("--manual-vmax", type=float, default=-1)
    ap.add_argument("--cmap", default="apex",
                    help="apex, viridis, magma, plasma, inferno, or comma-separated colors")
    ap.add_argument("--na-color", default="white")
    ap.add_argument("--bar-color", default="#8B1E3F")
    ap.add_argument("--candidate-id", default="")
    ap.add_argument("--group-counts", default="")
    ap.add_argument("--width", type=float, default=10.0)
    ap.add_argument("--height", type=float, default=None)
    ap.add_argument("--dpi", type=int, default=400)
    ap.add_argument("--formats", default="png")
    ap.add_argument("--title-size", type=int, default=12)
    ap.add_argument("--subtitle-size", type=int, default=8)
    ap.add_argument("--tick-size", type=int, default=9)
    ap.add_argument("--title-pad", type=float, default=20,
                    help="Padding above title in points")
    return ap.parse_args()

def load_cache(p):
    with open(p, "rb") as f: return pickle.load(f)

def get_cmap(name, na_color):
    if name == "apex":
        cols = ["#080808","#1a0a2e","#3d1261","#6a1b8a",
                "#9c2e7a","#cc4466","#e8734a","#f5a623"]
        cm = LinearSegmentedColormap.from_list("apex", cols, N=256)
    elif "," in name:
        cm = LinearSegmentedColormap.from_list("custom",
             [c.strip() for c in name.split(",")], N=256)
    else:
        cm = plt.get_cmap(name).copy()
    cm = cm.copy()
    cm.set_bad(color=na_color)
    return cm

def apply_transform(m, transform, power, k):
    if transform == "raw": return m
    elif transform == "power": return np.power(np.clip(m, 0, None), power)
    elif transform == "log1p": return np.log1p(k * np.clip(m, 0, None))
    return m

def compute_vmax(m, coef, quantile, manual):
    if manual > 0: return manual
    v = m[np.isfinite(m) & (m > 0)]
    if len(v) == 0: return 1.0
    # Use quantile-based vmax — more robust than median*coef
    return float(min(np.quantile(v, quantile) * coef, np.max(v)))

def extract_sub(cache, rs, re):
    bp = cache["bin_bp"]
    if rs is None: rs = 1
    if re is None: re = cache["chr_len"]
    bs = max(0, (rs - 1) // bp)
    be = min(cache["n_bins"], (re - 1) // bp + 1)
    return (cache["r2_matrix"][bs:be, bs:be].copy(),
            cache["count_matrix"][bs:be, bs:be].copy(), bs, be, bs * bp)

def nice_size(bp):
    if bp >= 1e6: return f"{bp/1e6:.2f} Mb"
    if bp >= 1e3: return f"{bp/1e3:.0f} kb"
    return f"{bp} bp"

def build_metadata_text(cache, rs, re, hs, he, bin_bp):
    parts = []
    if rs and re: parts.append(f"region: {rs:,}–{re:,}")
    if hs and he: parts.append(f"inv: {nice_size(he - hs)}")
    parts.append(f"bin: {nice_size(bin_bp)}")
    parts.append(f"type: {cache.get('cache_type', '?')}")
    rp = cache.get("retained_pairs", "?")
    parts.append(f"pairs: {rp:,}" if isinstance(rp, int) else f"pairs: {rp}")
    return " | ".join(parts)


# ═══════════════════════════════════════════════════════════════════════
# ROTATED TRIANGLE (using pcolormesh with rotated coordinates)
# ═══════════════════════════════════════════════════════════════════════

def plot_triangle(r2_sub, offset_bp, bin_bp, chrom, orientation,
                  transform, power, k, vmax_coef, vmax_q, manual_vmax,
                  cmap_name, na_color, bar_color,
                  rs, re, hs, he, cid, cache,
                  width, height, dpi, title_size, subtitle_size, tick_size,
                  title_pad, outpath):
    """Rotated triangle using coordinate transform for speed + quality."""
    n = r2_sub.shape[0]
    if n == 0: return

    # Only upper triangle (j >= i)
    mat = r2_sub.copy()
    mat[np.tril_indices(n, k=-1)] = np.nan

    display = np.ma.masked_invalid(mat)
    display = apply_transform(display, transform, power, k)
    vmax = compute_vmax(display.compressed(), vmax_coef, vmax_q, manual_vmax)
    cm = get_cmap(cmap_name, na_color)

    dx = bin_bp / 1e6
    x0 = offset_bp / 1e6

    if height is None:
        height = width * 0.5

    fig, ax = plt.subplots(figsize=(width, height))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Build rotated coordinate mesh
    # For bin (i,j): center_x = x0 + (i+j+1)/2 * dx, center_y = (j-i) * dx / 2
    # We build the mesh corners for pcolormesh
    xs = np.arange(n + 1) * dx + x0  # bin edges in genomic coords
    # Rotated coordinates: X = (col + row)/2, Y = (col - row)/2
    X = np.zeros((n + 1, n + 1))
    Y = np.zeros((n + 1, n + 1))
    for i in range(n + 1):
        for j in range(n + 1):
            X[i, j] = (xs[i] + xs[j]) / 2
            Y[i, j] = (xs[j] - xs[i]) / 2

    # For pcolormesh, we need the data aligned to the mesh
    # display[i,j] goes at mesh cell (i,j)
    pcm = ax.pcolormesh(X, Y, display, cmap=cm, vmin=0, vmax=vmax,
                        shading='flat', edgecolors='none', linewidth=0,
                        rasterized=True)

    x_min = x0
    x_max = x0 + n * dx
    y_max_val = n * dx / 2

    ax.set_xlim(x_min, x_max)

    if orientation == "triangle_down":
        ax.set_ylim(y_max_val * 1.02, -y_max_val * 0.02)  # inverted = pointing down
    else:  # triangle_up
        ax.set_ylim(-y_max_val * 0.02, y_max_val * 1.02)

    ax.set_aspect("equal")

    # Highlight lines
    if hs is not None and he is not None:
        h0, h1 = hs / 1e6, he / 1e6
        ax.axvline(h0, color=bar_color, lw=1.5, ls="--", alpha=0.8, zorder=5)
        ax.axvline(h1, color=bar_color, lw=1.5, ls="--", alpha=0.8, zorder=5)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cm, norm=Normalize(0, vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, shrink=0.6)
    tlab = {"raw": "r²", "power": "√r²" if power == 0.5 else f"r²^{power}",
            "log1p": f"log(1+{k}r²)"}
    cbar.set_label(tlab.get(transform, "r²"), fontsize=tick_size)

    ax.set_xlabel(f"{chrom} position (Mb)", fontsize=tick_size + 1)
    ax.set_ylabel("Pairwise distance (Mb)", fontsize=tick_size + 1)
    ax.tick_params(labelsize=tick_size)

    # Title with proper padding
    orient_label = "▼" if orientation == "triangle_down" else "▲"
    title = f"LD {orient_label} — {chrom}"
    if cid: title = f"Cand {cid} — {title}"
    ax.set_title(title, fontsize=title_size, fontweight="bold", pad=title_pad)

    # Metadata subtitle below title
    meta = build_metadata_text(cache, rs, re, hs, he, bin_bp)
    ax.text(0.5, 1.005, meta, transform=ax.transAxes,
            fontsize=subtitle_size, ha="center", color="grey")

    fig.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[LD-plot] {outpath}", file=sys.stderr)


# ═══════════════════════════════════════════════════════════════════════
# SQUARE HEATMAP (full or upper-half)
# ═══════════════════════════════════════════════════════════════════════

def plot_square(r2_sub, offset_bp, bin_bp, chrom, mode,
                transform, power, k, vmax_coef, vmax_q, manual_vmax,
                cmap_name, na_color, bar_color,
                rs, re, hs, he, cid, cache,
                width, height, dpi, title_size, subtitle_size, tick_size,
                title_pad, outpath):
    n = r2_sub.shape[0]
    if n == 0: return

    display = np.ma.masked_invalid(r2_sub.copy())
    if mode == "upper_half":
        display[np.tril_indices(n, k=-1)] = np.ma.masked

    display = apply_transform(display, transform, power, k)
    vmax = compute_vmax(display.compressed(), vmax_coef, vmax_q, manual_vmax)
    cm = get_cmap(cmap_name, na_color)

    if height is None: height = width

    x0 = offset_bp / 1e6
    x1 = x0 + n * bin_bp / 1e6

    fig, ax = plt.subplots(figsize=(width, height))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    im = ax.imshow(display, cmap=cm, norm=Normalize(0, vmax),
                   extent=[x0, x1, x1, x0], aspect="equal",
                   interpolation="nearest", origin="upper")

    if hs is not None and he is not None:
        h0, h1 = hs / 1e6, he / 1e6
        rect = Rectangle((h0, h0), h1 - h0, h1 - h0, lw=1.5,
                          edgecolor=bar_color, facecolor="none", ls="--", zorder=5)
        ax.add_patch(rect)
        bh = (x1 - x0) * 0.012
        ax.add_patch(Rectangle((max(h0, x0), x0 - bh * 2),
                                min(h1, x1) - max(h0, x0), bh,
                                fc=bar_color, ec='none', clip_on=False, zorder=10))
        ax.add_patch(Rectangle((x0 - bh * 2, max(h0, x0)),
                                bh, min(h1, x1) - max(h0, x0),
                                fc=bar_color, ec='none', clip_on=False, zorder=10))

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    tlab = {"raw": "r²", "power": "√r²" if power == 0.5 else f"r²^{power}",
            "log1p": f"log(1+{k}r²)"}
    cbar.set_label(tlab.get(transform, "r²"), fontsize=tick_size)

    ax.set_xlabel(f"{chrom} (Mb)", fontsize=tick_size + 1)
    ax.set_ylabel(f"{chrom} (Mb)", fontsize=tick_size + 1)
    ax.tick_params(labelsize=tick_size)

    mode_l = "Full" if mode == "full" else "Upper half"
    title = f"LD ({mode_l}) — {chrom}"
    if cid: title = f"Cand {cid} — {title}"
    ax.set_title(title, fontsize=title_size, fontweight="bold", pad=title_pad)

    meta = build_metadata_text(cache, rs, re, hs, he, bin_bp)
    ax.text(0.5, 1.005, meta, transform=ax.transAxes,
            fontsize=subtitle_size, ha="center", color="grey")

    fig.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[LD-plot] {outpath}", file=sys.stderr)


def main():
    args = parse_args()
    cache = load_cache(args.cache)
    r2_sub, cnt_sub, bs, be, offset = extract_sub(
        cache, args.region_start, args.region_end)
    rs = args.region_start or 1
    re = args.region_end or cache["chr_len"]
    fmts = [f.strip() for f in args.formats.split(",")]

    for fmt in fmts:
        outpath = f"{args.outprefix}.{args.mode}.{fmt}"
        os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)

        kw = dict(offset_bp=offset, bin_bp=cache["bin_bp"], chrom=cache["chrom"],
                  transform=args.transform, power=args.power, k=args.log1p_k,
                  vmax_coef=args.vmax_coef, vmax_q=args.vmax_quantile,
                  manual_vmax=args.manual_vmax,
                  cmap_name=args.cmap, na_color=args.na_color, bar_color=args.bar_color,
                  rs=rs, re=re, hs=args.highlight_start, he=args.highlight_end,
                  cid=args.candidate_id, cache=cache,
                  width=args.width, height=args.height, dpi=args.dpi,
                  title_size=args.title_size, subtitle_size=args.subtitle_size,
                  tick_size=args.tick_size, title_pad=args.title_pad,
                  outpath=outpath)

        if args.mode in ("triangle_down", "triangle_up"):
            plot_triangle(r2_sub, orientation=args.mode, **kw)
        else:
            plot_square(r2_sub, mode=args.mode, **kw)

if __name__ == "__main__":
    main()
