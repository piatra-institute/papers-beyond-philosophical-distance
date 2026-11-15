"""Figures for the manuscript, drawn from the results dictionary of the recorded run."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FuncFormatter, NullFormatter  # noqa: E402

plt.rcParams.update({
    "font.size": 15, "axes.titlesize": 16, "axes.labelsize": 15, "legend.fontsize": 13,
    "xtick.labelsize": 14, "ytick.labelsize": 14, "pdf.fonttype": 42, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False,
})

COL = {"categorical": "#1b9e77", "excl": "#66a61e", "lin": "#d95f02", "pair": "#7570b3",
       "unres": "#e7298a", "impl": "#666666", "cont": "#1f78b4", "grid": "#e31a1c", "box": "#33a02c"}


def _power_label(ax, x, y, base: str, exponent: str, rest: str, color: str, size: float = 13,
                 exp_size: float = 12, second_line: str | None = None) -> None:
    """Direct label "base^exponent rest" with the exponent drawn as raised text.

    Mathtext sets exponents at 70% of the label size, which prints below 7 pt once the
    figure is scaled to the text width; a separate text keeps the exponent at exp_size.
    The annotations are not clipped, so their extents can be measured before the axis
    limits are final.
    """
    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    pt = 72.0 / fig.dpi
    t = ax.text(x, y, base, color=color, fontsize=size, va="center", ha="left")
    w = t.get_window_extent(renderer).width * pt
    e = ax.annotate(exponent, (x, y), xytext=(w + 0.4, 0.42 * size), textcoords="offset points",
                    color=color, fontsize=exp_size, fontstyle="italic", va="center", ha="left",
                    annotation_clip=False)
    we = e.get_window_extent(renderer).width * pt
    ax.annotate(rest, (x, y), xytext=(w + we + 0.4 * size, 0), textcoords="offset points",
                color=color, fontsize=size, va="center", ha="left", annotation_clip=False)
    if second_line:
        ax.annotate(second_line, (x, y), xytext=(1.6 * size, -1.25 * size), textcoords="offset points",
                    color=color, fontsize=size, va="center", ha="left", annotation_clip=False)


def figure_orders(res: dict, path: Path) -> None:
    th = res["theory"]
    fig = plt.figure(figsize=(12, 10.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], width_ratios=[1.5, 1], hspace=0.38, wspace=0.28)
    ax = fig.add_subplot(gs[0, 0])
    rows = th["boolean"]["rows"]
    specs = [
        ("unrestricted", None, COL["unres"], lambda n: 2.0 ** n, "D", 5),
        ("pairwise_clauses", "$2n$ two-literal clauses", COL["pair"], lambda n: 2.0 * n, "s", 6),
        ("horn", "$n+1$ answers + implications,\n   Horn rules, parity", COL["lin"], lambda n: n + 1.0, "o", 6),
        ("categorical_exclusions", "3  answers + exclusions", COL["excl"], lambda n: 3.0 + 0 * n, "^", 6),
        ("categorical", "2  categorical answers", COL["categorical"], lambda n: 2.0 + 0 * n, "v", 6),
        ("implications_only", "1  implications alone", COL["impl"], lambda n: 1.0 + 0 * n, "x", 6),
    ]
    for key, label, color, f, marker, nmax in specs:
        nn = np.arange(2, nmax + 1)
        ax.plot(nn, f(nn), color=color, lw=1.6, ls="--", alpha=0.75)
        pts = [(r["issues"], r["order_by_points"]) for r in rows if r["class"] == key and r["issues"] >= 2]
        ax.plot([q[0] for q in pts], [q[1] for q in pts], marker=marker, ls="none", color=color, ms=9)
        if label is None:
            _power_label(ax, nmax + 0.15, f(np.array([nmax]))[0], "2", "n", "unrestricted", color)
        else:
            ax.text(nmax + 0.15, f(np.array([nmax]))[0], label, color=color, fontsize=13, va="center")
    for key, marker in (("categorical_implications", "P"), ("parity", "X")):
        pts = [(r["issues"], r["order_by_points"]) for r in rows if r["class"] == key and r["issues"] >= 2]
        ax.plot([q[0] + 0.14 for q in pts], [q[1] for q in pts], marker=marker, ls="none", color=COL["lin"], ms=8)
    ax.set_yscale("log", base=2)
    ax.set_ylim(0.8, 40)
    ax.set_xlim(1.7, 9.9)
    ax.set_xlabel("binary issues $n$")
    ax.set_ylabel("obstruction order")
    ax.set_title("(a) Yes/no commitments", loc="left")
    ax.set_xticks(range(2, 7))
    ax.set_yticks([1, 2, 4, 8, 16, 32])
    ax.set_yticklabels(["1", "2", "4", "8", "16", "32"])

    ax = fig.add_subplot(gs[0, 1])
    d = np.arange(1, 6)
    ax.plot(d, 2.0 ** d, color=COL["grid"], lw=1.6, ls="--")
    ax.plot(d, d + 1.0, color=COL["cont"], lw=1.6, ls="--")
    ax.plot(d, 2.0 + 0 * d, color=COL["box"], lw=1.6, ls="--")
    _power_label(ax, 5.15, 32, "2", "d", "convex,", COL["grid"], second_line="on a grid")
    ax.text(5.15, 6, "$d+1$ convex,\n   continuum", color=COL["cont"], fontsize=13, va="center")
    ax.text(5.15, 2, "2  boxes", color=COL["box"], fontsize=13, va="center")
    gc = th["witnesses"]["grid_cube"]
    ax.plot([r["dimension"] for r in gc], [r["family_size"] for r in gc], "o", color=COL["grid"], ms=9,
            label="verified lower-bound family")
    sp = th["single_issue_and_products"]["rows"]
    ex = [r["order"] for r in sp if r["case"] == "grid_convex"]
    bx = [r["order"] for r in sp if r["case"] == "grid_boxes"]
    ax.plot([2.1] * len(ex), ex, "*", color="black", ms=13, label="exhaustive, $L\\times L$ grids")
    ax.plot([2.1] * len(bx), bx, "*", color="black", ms=13)
    ax.set_yscale("log", base=2)
    ax.set_ylim(1.6, 40)
    ax.set_xlim(0.7, 7.6)
    ax.set_yticks([2, 4, 8, 16, 32])
    ax.set_yticklabels(["2", "4", "8", "16", "32"])
    ax.set_xticks(range(1, 6))
    ax.set_xlabel("quality dimensions $d$")
    ax.set_ylabel("obstruction order")
    ax.set_title("(b) Regions of a conceptual space", loc="left")

    ax = fig.add_subplot(gs[1, :])
    pr = th["planar"]["rows"]
    m = np.array([r["positions"] for r in pr])
    for r in pr:
        for k, v in r["histogram"].items():
            ax.scatter([r["positions"]], [int(k)], s=12 + 6 * v, color="#999999", alpha=0.55, zorder=1)
    ax.scatter([], [], s=60, color="#999999", alpha=0.55, label="random sets (area = count)")
    ax.plot(m, [r["mean_order"] for r in pr], "o-", color="black", ms=7, lw=1.8, label="mean over random sets", zorder=3)
    cp = th["planar"]["convex_position"]
    ax.plot([r["positions"] for r in cp], [r["order"] for r in cp], "s--", color=COL["unres"], ms=7,
            label="convex position (order $= m$)")
    ax.axhline(3, color=COL["cont"], ls=":", lw=1.8, label="continuum plane (3)")
    ax.set_xlabel("occupied positions $m$ in the plane")
    ax.set_ylabel("obstruction order")
    ax.set_title("(c) Convex regions evaluated on the occupied positions of a planar map", loc="left")
    ax.set_xticks(range(int(m.min()), int(m.max()) + 1))
    ax.legend(loc="upper left", frameon=False, fontsize=13, ncol=2)
    ax.set_ylim(2.5, 13.5)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def figure_data(res: dict, path: Path) -> None:
    tr = res["trolley"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.0))
    ax = axes[0]
    labels = [c["pattern_switch_loop_footbridge"] for c in tr["cells"]]
    x = np.arange(len(labels))
    ratio = np.array([c["ratio"] for c in tr["cells"]])
    lo = np.array([c["ratio_ci95"][0] for c in tr["cells"]])
    hi = np.array([c["ratio_ci95"][1] for c in tr["cells"]])
    ax.bar(x, ratio - 1.0, bottom=1.0, color=["#d95f02" if r < 1 else "#1b9e77" for r in ratio], alpha=0.85)
    ax.errorbar(x, ratio, yerr=[ratio - lo, hi - ratio], fmt="none", ecolor="black", capsize=4, lw=1.4)
    ax.axhline(1.0, color="black", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=14)
    ax.set_xlabel("answers to Switch, Loop, Footbridge (1 = sacrifice)")
    ax.set_ylabel("observed / pairwise prediction")
    ax.set_title("(a) Trolley judgments, $N$ = {:,}".format(tr["participants"]), loc="left")

    ax = axes[1]
    pos = res["consciousness"]["positions"]
    pw = pos["power"]
    for r in pw["rows"]:
        xs = [c["third_order_bits"] for c in r["curve"]]
        ax.plot(xs, [c["power_alpha_005"] for c in r["curve"]], color="#1f78b4", alpha=0.45, lw=1.3)
    dp = pos["discovery_power"]
    lv = [v for v in dp["levels"] if v["datasets"] > 0]
    xb = np.array([v["target_bits"] for v in lv])
    yb = np.array([v["power"] for v in lv])
    lo_b = np.array([v["power_ci95"][0] for v in lv])
    hi_b = np.array([v["power_ci95"][1] for v in lv])
    ax.errorbar(xb, yb, yerr=[yb - lo_b, hi_b - yb], fmt="s-", color="#e31a1c", ms=7, lw=1.6, capsize=3)
    trolley_bits = tr["information"]["third_order_bits"]
    ax.axvline(trolley_bits, color="black", ls="--", lw=1.5)
    tt = pos["triplet_tests"]
    ax.axvline(tt["max_third_order_bits"], color="#666666", ls=":", lw=1.8)
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 0.3)
    ax.set_ylim(0, 1.02)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.text(trolley_bits * 0.9, 0.97, "trolley effect", rotation=90, ha="right", va="top", fontsize=14)
    ax.text(tt["max_third_order_bits"] * 0.94, 0.3, "largest triplet", rotation=90, ha="right",
            va="bottom", fontsize=14, color="#444444")
    ax.set_xlabel("planted third-order information (bits)")
    ax.set_ylabel("detection probability")
    ax.set_title("(b) Power at $N$ = %d" % pos["respondents_complete"], loc="left")
    handles = [Line2D([0], [0], color="#1f78b4", lw=1.6),
               Line2D([0], [0], color="#e31a1c", marker="s", ms=7, lw=1.6)]
    fig.legend(handles, ["one pre-specified triplet, $\\alpha$ = 0.05",
                         "screening of all %d triplets" % tt["triplets"]],
               loc="lower center", bbox_to_anchor=(0.765, 0.0), ncol=1, frameon=False, fontsize=14)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
