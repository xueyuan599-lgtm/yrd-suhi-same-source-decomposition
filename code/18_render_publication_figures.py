"""Render publication figures from the frozen same-source analysis tables.

Run after 16_build_revision_outputs.py. No estimates are refitted here.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "outputs" / "daily"
FIGS = BASE / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

BLUE = "#0072B2"
ORANGE = "#D55E00"
PURPLE = "#9467BD"
GRAY = "#8C959D"
LIGHT = "#D7E5ED"
DARK = "#216B83"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 10,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 8,
    "pdf.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def interval(ax, mean, low, high, y, color, marker="o", size=5):
    if not (np.isfinite(mean) and np.isfinite(low) and np.isfinite(high)):
        return
    ax.errorbar(
        mean, y, xerr=[[mean - low], [high - mean]], fmt=marker,
        color=color, ecolor=color, markersize=size, capsize=2.2,
        linewidth=1.5, markeredgewidth=0, zorder=4,
    )


def style_contrast_axis(ax, title, labels, xlim=(-0.86, 0.80)):
    ax.axvline(0, color="#59636D", linewidth=0.8, linestyle="--", zorder=1)
    ax.set_xlim(*xlim)
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.invert_yaxis()
    ax.set_title(title, loc="left", fontweight="bold", pad=8)
    ax.grid(axis="x", color="#E6EAED", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)


def render_figure2():
    boot = pd.read_csv(OUT / "same_source_cluster_bootstrap.csv")
    sensitivity = pd.read_csv(OUT / "same_source_revision_sensitivity_table.csv")
    expected = {"A_minus_B", "B_minus_C", "C_minus_D"}
    assert set(boot.metric) == expected
    assert set(boot["mode"]) == {"day", "night"}

    fig = plt.figure(figsize=(8.6, 6.7), constrained_layout=False)
    # 左列刻度标签较长（`≥40% mask` 等）。left=.13 / wspace=.27 时，左列标签被画布左缘截断，
    # 右列标签又会压到左图的数据区上（2026-10-04 目视实测）。这里把左边距与列间距留足。
    grid = fig.add_gridspec(2, 2, left=.175, right=.985, top=.875, bottom=.13,
                           wspace=.46, hspace=.40, height_ratios=[.9, 1.65])
    step_names = ["A−B  dates", "B−C  rebuild", "C−D  native"]
    metric_info = [("A_minus_B", BLUE), ("B_minus_C", ORANGE), ("C_minus_D", PURPLE)]
    scenarios = [
        ("20 px primary", "20px primary"),
        ("≥1 pixel", "≥1 pixel"),
        ("≥40% mask", "≥40% mask pixels"),
        ("≥60% mask", "≥60% mask pixels"),
        ("2020 mask", "2020 fixed mask"),
        ("10 paired days", "10 paired days/month"),
    ]

    for col, (mode, title) in enumerate((("day", "Daytime"), ("night", "Nighttime"))):
        top = fig.add_subplot(grid[0, col])
        style_contrast_axis(top, f"{chr(97 + col)}  {title}: primary steps", step_names)
        for j, (metric, color) in enumerate(metric_info):
            row = boot[(boot["mode"] == mode) & (boot.metric == metric)].iloc[0]
            interval(top, row.mean_K, row.cluster_bootstrap_95_lo_K,
                     row.cluster_bootstrap_95_hi_K, j, color, size=5.5)
        top.set_xlabel("Direct difference (K)")
        # 城市数移到面板内右下空白：原位置 (1, 1.08) 会与总标题及分面板标题互相压叠。
        top.text(.985, .06, f"{int(boot[boot['mode'] == mode].cities.iloc[0])} cities",
                 transform=top.transAxes, ha="right", va="bottom",
                 fontsize=7.5, color="#4A5259")

        lower = fig.add_subplot(grid[1, col])
        tick_labels = []
        rows = []
        for label, scenario in scenarios:
            row = sensitivity[(sensitivity["mode"] == mode) &
                              (sensitivity.scenario == scenario)].iloc[0]
            rows.append(row)
            tick_labels.append(label)
        style_contrast_axis(lower, f"{chr(99 + col)}  {title}: observation rules", tick_labels)
        for i, row in enumerate(rows):
            if i in (2, 3):
                lower.axhspan(i - .49, i + .49, color="#F4F6F7", zorder=-2)
            if int(row.city_years) == 0:
                # 无样本行：把说明与样本量合成一个标签居中放置。分开画时 "not estimable"
                # 会与右端的 n=0 叠印（2026-10-04 目视实测），故不再单独画 n=。
                lower.text(.72, i, "not estimable (n = 0)", ha="center", va="center",
                           transform=lower.get_yaxis_transform(),
                           fontsize=7.8, color="#68737D", style="italic")
                continue
            # 支持量写在面板内侧右端；并进刻度标签会把较长标签挤出画布。
            lower.text(.985, i, f"n={int(row.city_years)}",
                       transform=lower.get_yaxis_transform(), ha="right", va="center",
                       fontsize=7, color="#4A5259", zorder=5)
            for prefix, color, marker, offset in (
                ("delta_AB_alt_K", BLUE, "o", -.13),
                ("delta_BC_alt_K", ORANGE, "s", .13),
            ):
                interval(lower, row[prefix + "_mean"], row[prefix + "_lo"],
                         row[prefix + "_hi"], i + offset, color, marker, size=4.5)
        lower.set_xlabel("Matched monthly level difference (K)")

    handles = [Line2D([0], [0], color=BLUE, marker="o", lw=1.4, label="A−B: date sets"),
               Line2D([0], [0], color=ORANGE, marker="s", lw=1.4,
                      label="B−C: 8-day reconstruction"),
               Line2D([0], [0], color=PURPLE, marker="o", lw=1.4,
                      label="C−D: paired product output")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(.58, .025), fontsize=7.6)
    # 旧标题 "Measured changes…" 容易被读成"实测变化"；改为与内容一致的估计量差异。
    fig.text(.175, .99, "Estimator differences and their observation support",
             ha="left", va="top", fontsize=12, weight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"fig2_00_same_source_decomposition.{ext}",
                    dpi=400 if ext == "png" else None)
    plt.close(fig)


def render_figure3():
    annual = pd.read_csv(OUT / "same_source_decomposition_city_year.csv",
                         usecols=["basis", "mode", "year", "city_id", "pixel_rule",
                                  "min_daily_days", "min_8day_periods", "months_available"])
    q = annual[(annual.basis == "dynamic") & (annual.pixel_rule == "20px") &
               (annual.min_daily_days == 6) & (annual.min_8day_periods == 2)].copy()
    assert not q.duplicated(["mode", "year", "city_id"]).any()
    years = np.arange(2001, 2026)
    fig, axes = plt.subplots(2, 1, figsize=(8.1, 4.8), sharex=True,
                             gridspec_kw={"hspace": .20})
    cats = [("0–5 months", lambda x: x <= 5, GRAY),
            ("6–11 months", lambda x: (x >= 6) & (x <= 11), LIGHT),
            ("12 months", lambda x: x == 12, DARK)]
    for ax, mode, title in zip(axes, ("day", "night"), ("Daytime", "Nighttime")):
        g = q[q["mode"] == mode]
        counts = g.groupby("year").size().reindex(years, fill_value=0)
        assert counts.min() > 0
        bottom = np.zeros(len(years))
        for label, test, color in cats:
            k = g.assign(in_cat=test(g.months_available)).groupby("year").in_cat.sum()
            frac = (100 * k.reindex(years, fill_value=0) / counts).to_numpy()
            ax.bar(years, frac, bottom=bottom, width=.80, color=color,
                   edgecolor="white", linewidth=.30, label=label)
            bottom += frac
        assert np.allclose(bottom, 100)
        full = int((g.months_available == 12).sum())
        assert full == (253 if mode == "day" else 103)
        ax.set_ylim(0, 103)
        ax.set_yticks([0, 25, 50, 75, 100], ["0", "25", "50", "75", "100%"])
        ax.set_ylabel(title, fontweight="bold", labelpad=9)
        ax.text(.99, 1.03, f"12 months: {full}/{len(g)} city-years ({full/len(g):.1%})",
                transform=ax.transAxes, ha="right", fontsize=8)
        ax.grid(axis="y", color="#E4E8EB", linewidth=.6)
        ax.set_axisbelow(True)
    axes[-1].set_xticks(years[::2])
    axes[-1].set_xticklabels(years[::2], rotation=0)
    axes[-1].set_xlabel("Year")
    axes[0].legend(handles=[Patch(facecolor=c, edgecolor="none", label=l)
                            for l, _, c in cats], ncol=3, frameon=False,
                   loc="upper left", bbox_to_anchor=(0, 1.25), fontsize=7.7)
    fig.subplots_adjust(left=.095, right=.99, top=.88, bottom=.14)
    fig.text(.095, .985, "Annual observation support is incomplete",
             ha="left", va="top", fontsize=12, weight="bold")
    fig.text(.095, .03, "Each bar is the share of candidate city-years with that many estimable months; yearly denominators vary from 24 to 37.",
             ha="left", fontsize=7.4, color="#4A5259")
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"fig3_00_complete_month_support.{ext}",
                    dpi=400 if ext == "png" else None)
    plt.close(fig)


if __name__ == "__main__":
    render_figure2()
    render_figure3()
    print("Rendered Figure 2 and Figure 3 from frozen analysis CSVs.")
