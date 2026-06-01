#!/usr/bin/env python3
# =============================================================================
# plot_pl_sensitivity.py
#
# Three publication-quality figures from the tidy frame:
#   FIG 1 (MAIN) : persistence (n_final / n_founding) -- rescue/extinction
#   FIG 2 (COMP) : fitness_lag -- residual maladaptation among survivors
#   FIG 3 (COMP) : fitness_end -- demographic-cost framing
#
# SHARED GRAMMAR (every figure):
#   * rows  = lineage (pdav-N, mixed, rot), in steepness order
#   * x     = migration rate (log-spaced: 1e-6, 1e-5, 1e-4)
#   * colour/offset = script (with vs without introgression)
#   * ALL TEN rep points shown (jittered), with mean +/- 95% CI overlaid.
#     Per the standing rule: never means-only; central tendency always carried
#     with its spread, and the raw points always visible beneath it.
#   * pl_mult handled by a small multiple within each row (3 panels), so the
#     full 3 (lineage) x 3 (pl) grid is visible with migration on x.
#
# SCAFFOLDING BANNER: if the data is synthetic, a clear watermark says so, so
# nobody mistakes the placeholder for a result.
# =============================================================================

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---- house style: clean, journal-ready, light theme -------------------------
mpl.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 300,
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "Georgia"],
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "axes.edgecolor": "#33312e",
    "xtick.color": "#33312e", "ytick.color": "#33312e",
    "axes.labelcolor": "#26241f", "text.color": "#26241f",
    "axes.grid": True, "grid.color": "#e7e2d9", "grid.linewidth": 0.7,
    "axes.axisbelow": True, "figure.facecolor": "white", "axes.facecolor": "#fcfbf7",
})

# ecotype-adjacent palette: introgression warm, no-introgression cool
C_WITH    = "#c2683d"   # terracotta  (with introgression)
C_WITHOUT = "#4f6f8f"   # slate blue  (without introgression)
LINEAGES  = ["pdav-N", "mixed", "rot"]
PL_ORDER  = ["0.5", "1", "2"]
MIG_ORDER = ["1e-06", "1e-05", "0.0001"]
MIG_LABEL = {"1e-06": r"$10^{-6}$", "1e-05": r"$10^{-5}$", "0.0001": r"$10^{-4}$"}


def ci95(vals):
    """Mean and half-width of a 95% CI (normal approx). Returns (mean, half)."""
    vals = np.asarray(vals, float)
    n = len(vals)
    if n < 2:
        return (float(vals.mean()) if n else np.nan, 0.0)
    m = vals.mean(); se = vals.std(ddof=1) / np.sqrt(n)
    return m, 1.96 * se


def _panel(ax, sub, ycol, ylim=None):
    """Draw one panel: x = migration, two scripts side by side, points + mean/CI."""
    xbase = np.arange(len(MIG_ORDER))
    rng = np.random.default_rng(7)
    for k, (scr, col, dx) in enumerate([("with", C_WITH, -0.13),
                                        ("without", C_WITHOUT, 0.13)]):
        means, halves, xs = [], [], []
        for j, mig in enumerate(MIG_ORDER):
            vals = sub[(sub["script"] == scr) & (sub["mig"] == mig)][ycol].values
            x = xbase[j] + dx
            xs.append(x)
            jit = rng.uniform(-0.045, 0.045, size=len(vals))
            ax.scatter(np.full(len(vals), x) + jit, vals, s=15, color=col,
                       alpha=0.5, edgecolor="none", zorder=2)
            m, h = ci95(vals); means.append(m); halves.append(h)
        ax.errorbar(xs, means, yerr=halves, fmt="o", color=col, ms=5.5,
                    capsize=3, lw=1.4, mec="white", mew=0.7, zorder=4)
        ax.plot(xs, means, color=col, lw=1.2, alpha=0.85, zorder=3)
    ax.set_xticks(xbase); ax.set_xticklabels([MIG_LABEL[m] for m in MIG_ORDER])
    if ylim:
        ax.set_ylim(*ylim)


def make_figure(df, ycol, title, ylabel, outpath, synthetic, ylim=None, ref_line=None):
    nrow, ncol = len(LINEAGES), len(PL_ORDER)
    fig, axes = plt.subplots(nrow, ncol, figsize=(9.8, 8.8),
                             sharex=True, sharey="row")
    for i, lin in enumerate(LINEAGES):
        for j, pl in enumerate(PL_ORDER):
            ax = axes[i, j]
            sub = df[(df["lineage"] == lin) & (df["pl"] == pl)]
            _panel(ax, sub, ycol, ylim=ylim)
            if ref_line is not None:
                ax.axhline(ref_line, ls=":", lw=0.8, color="#9a9488", zorder=1)
            if i == 0:
                ax.set_title(f"plasticity \u00d7{pl}", pad=8)
            if j == 0:
                ax.set_ylabel(ylabel, fontsize=9)
            if i == nrow - 1:
                ax.set_xlabel("migration rate")
        # one row label per lineage, placed to the LEFT of the row (no overlap)
        ymid = axes[i, 0].get_position().y0 + \
            (axes[i, 0].get_position().y1 - axes[i, 0].get_position().y0) / 2

    # row labels as rotated text in the far-left margin, after layout is set
    # (done post-tight_layout below to use final positions)

    # legend -- sits in its own band ABOVE the subtitle, no collision
    handles = [Line2D([0], [0], marker="o", color=C_WITH, lw=1.2, ms=6,
                      mec="white", label="with introgression"),
               Line2D([0], [0], marker="o", color=C_WITHOUT, lw=1.2, ms=6,
                      mec="white", label="without introgression")]

    fig.tight_layout(rect=[0.06, 0, 1, 0.88])

    # place header elements using figure coords AFTER layout
    fig.suptitle(title, y=0.985, fontsize=13, fontweight="bold")
    fig.text(0.53, 0.93,
             "rows: lineage (steepest climate shift at top)   |   "
             "points: 10 replicates, jittered   |   markers: mean \u00b1 95% CI",
             ha="center", fontsize=8.5, color="#6b6559")
    fig.legend(handles=handles, loc="center", ncol=2, frameon=False,
               bbox_to_anchor=(0.53, 0.905), fontsize=10)

    # rotated lineage row-labels in the left margin
    for i, lin in enumerate(LINEAGES):
        pos = axes[i, 0].get_position()
        ymid = pos.y0 + (pos.y1 - pos.y0) / 2
        fig.text(0.018, ymid, lin, rotation=90, va="center", ha="center",
                 fontsize=12, fontweight="bold", color="#3a372f")

    if synthetic:
        fig.text(0.5, 0.45, "SCAFFOLDING  \u00b7  SYNTHETIC PLACEHOLDER DATA",
                 ha="center", va="center", fontsize=26, color="#c9302c",
                 alpha=0.11, rotation=24, fontweight="bold", zorder=0)

    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    return outpath


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tidy", default="pl_sensitivity_tidy.csv")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--synthetic", action="store_true",
                    help="Stamp a SCAFFOLDING watermark (use for placeholder data).")
    args = ap.parse_args()

    df = pd.read_csv(args.tidy, dtype={"pl": str, "mig": str, "rep": str})
    df["pl"]  = pd.Categorical(df["pl"], PL_ORDER, ordered=True)
    df["mig"] = pd.Categorical(df["mig"], MIG_ORDER, ordered=True)
    df["lineage"] = pd.Categorical(df["lineage"], LINEAGES, ordered=True)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)   # savefig won't create dirs
    made = []
    made.append(make_figure(
        df, "persistence",
        "Fig 1 (MAIN). Does introgression rescue populations from climate collapse?",
        "persistence  (n_final / n_founding)",
        outdir / "fig1_persistence_MAIN.png", args.synthetic,
        ylim=(-0.03, 1.03)))
    made.append(make_figure(
        df, "fitness_lag",
        "Fig 2. Residual maladaptation among survivors (fitness lag)",
        "fitness lag  (optimum \u2212 phenotype)",
        outdir / "fig2_fitness_lag.png", args.synthetic,
        ref_line=0.0))
    made.append(make_figure(
        df, "fitness_end",
        "Fig 3. End-point mean fitness (demographic-cost framing)",
        "mean fitness at final generation",
        outdir / "fig3_fitness_end.png", args.synthetic,
        ylim=(-0.03, 1.03)))
    for p in made:
        print("wrote", p)


if __name__ == "__main__":
    main()
