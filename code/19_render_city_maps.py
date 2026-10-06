"""Render Figure 4 — city-level maps of the same-source monthly level differences.

Run after 16_build_revision_outputs.py. No estimates are refitted here: every plotted
value is the within-city mean of a field already written by 16 (city-years 2001–2025,
primary rule = dynamic basis / 20 px / 6 paired days / 2 periods).

Aggregation matches Table 2: average repeated city-years inside each city first, then
treat cities equally. The panel means printed at the end are asserted against
same_source_cluster_bootstrap.csv, so the map cannot drift away from the main table.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from matplotlib.colors import ListedColormap


BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "outputs" / "daily"
FIGS = BASE / "figures"
FIGS.mkdir(parents=True, exist_ok=True)
CITY_TIF = BASE / "data" / "gee_export" / "uhi_cityid.tif"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.4,
    "ytick.labelsize": 7.4,
    "pdf.fonttype": 42,
})

CMAP = plt.get_cmap("RdBu_r").copy()
CMAP.set_bad(alpha=0.0)
BASE_GRAY = ListedColormap(["#E7EAEC"])


def load_city_means() -> pd.DataFrame:
    cols = ["basis", "year", "mode", "city_id", "pixel_rule", "min_daily_days",
            "min_8day_periods", "delta_A_to_B_common_months_K",
            "delta_B_to_C_common_months_K"]
    df = pd.read_csv(OUT / "same_source_decomposition_city_year.csv", usecols=cols)
    q = df[(df.basis == "dynamic") & (df.pixel_rule == "20px") &
           (df.min_daily_days == 6) & (df.min_8day_periods == 2)].copy()
    # 与 Table 2 相同的支持口径：该表只对"四种估计量共有的月份口径差"取均值，
    # 因此夜间 803 行中只有 688 行有可用的共同月份差值（34 城）；白天 803 行全有（37 城）。
    q = q[q.delta_A_to_B_common_months_K.notna() & q.delta_B_to_C_common_months_K.notna()]
    assert len(q[q["mode"] == "day"]) == 803 and q[q["mode"] == "day"].city_id.nunique() == 37
    assert len(q[q["mode"] == "night"]) == 688 and q[q["mode"] == "night"].city_id.nunique() == 34
    return q.groupby(["mode", "city_id"])[
        ["delta_A_to_B_common_months_K", "delta_B_to_C_common_months_K"]].mean()


def check_against_bootstrap(means: pd.DataFrame) -> None:
    boot = pd.read_csv(OUT / "same_source_cluster_bootstrap.csv")
    pairs = {"A_minus_B": "delta_A_to_B_common_months_K",
             "B_minus_C": "delta_B_to_C_common_months_K"}
    for mode in ("day", "night"):
        for metric, col in pairs.items():
            expect = float(boot[(boot["mode"] == mode) & (boot.metric == metric)].mean_K.iloc[0])
            got = float(means.loc[mode][col].mean())
            assert abs(got - expect) < 1e-3, (mode, metric, got, expect)
            print(f"  {mode:5s} {metric}: city-mean {got:+.3f} K "
                  f"(bootstrap table {expect:+.3f} K)")


def decorate_map(ax, extent) -> None:
    """比例尺 + 指北针 + 经纬刻度——该刊对地图类图的硬要求。"""
    left, right, bottom, top = extent
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)
    ax.set_xticks(np.arange(115, 123, 1))
    ax.set_yticks(np.arange(28, 36, 1))
    ax.set_xticklabels([f"{t}°E" for t in np.arange(115, 123, 1)])
    ax.set_yticklabels([f"{t}°N" for t in np.arange(28, 36, 1)])
    ax.tick_params(length=2, pad=1.5)
    for s in ax.spines.values():
        s.set_linewidth(0.5)
        s.set_color("#7A8186")

    # 比例尺：31°N 处 1° 经度 ≈ 111.32·cos(31°) ≈ 95.4 km，图上标明该近似纬度。
    bar_km, bar_lat = 100.0, 31.0
    bar_deg = bar_km / (111.32 * np.cos(np.deg2rad(bar_lat)))
    x0, y0 = left + 0.30, bottom + 0.45
    ax.plot([x0, x0 + bar_deg], [y0, y0], color="#2B2B2B", lw=2.0,
            solid_capstyle="butt", zorder=6)
    for x in (x0, x0 + bar_deg):
        ax.plot([x, x], [y0, y0 + 0.12], color="#2B2B2B", lw=1.0, zorder=6)
    ax.text(x0 + bar_deg / 2, y0 + 0.16, f"{bar_km:.0f} km", ha="center", va="bottom",
            fontsize=6.2, color="#2B2B2B", zorder=6)

    ax.annotate("N", xy=(right - 0.42, bottom + 0.62), xytext=(right - 0.42, bottom + 0.28),
                ha="center", va="center", fontsize=6.4, color="#2B2B2B", zorder=6,
                arrowprops=dict(arrowstyle="-|>", color="#2B2B2B", lw=1.0,
                                shrinkA=1.5, shrinkB=1.5))
    ax.set_aspect("equal")


def main() -> None:
    means = load_city_means()
    print("对照 same_source_cluster_bootstrap.csv：")
    check_against_bootstrap(means)

    with rasterio.open(CITY_TIF) as ds:
        ids = ds.read(1)
        tr = ds.transform
    assert ids.shape == (885, 855), ids.shape
    ext = (tr.c, tr.c + tr.a * ids.shape[1], tr.f + tr.e * ids.shape[0], tr.f)
    city_base = np.where(ids > 0, 1.0, np.nan)

    panels = [("day", "delta_A_to_B_common_months_K", "(a)  Daytime: A−B"),
              ("night", "delta_A_to_B_common_months_K", "(b)  Nighttime: A−B"),
              ("day", "delta_B_to_C_common_months_K", "(c)  Daytime: B−C"),
              ("night", "delta_B_to_C_common_months_K", "(d)  Nighttime: B−C")]

    maps = {}
    for mode, col, title in panels:
        v = np.full(ids.shape, np.nan)
        for cid, val in means.loc[mode][col].items():
            v[ids == cid] = val
        maps[(mode, col)] = v

    pooled = np.concatenate([m[np.isfinite(m)] for m in maps.values()])
    vmax = np.ceil(max(np.percentile(np.abs(pooled), 98), 0.5) / 0.25) * 0.25
    print(f"共用色标：0 ± {vmax:.2f} K（覆盖 |值| 的 98 分位）")

    fig = plt.figure(figsize=(7.2, 6.8))
    # top=.925 / bottom=.075：面板标题画在轴上方，留出总标题与脚注的行高，
    # 否则总标题会压住 (a)(b) 的标题（2026-10-04 第一次出图即踩到）。
    grid = fig.add_gridspec(2, 2, left=.07, right=.885, top=.925, bottom=.075,
                            wspace=.16, hspace=.22)
    for k, (mode, col, title) in enumerate(panels):
        ax = fig.add_subplot(grid[k // 2, k % 2])
        ax.imshow(city_base, extent=ext, origin="upper", cmap=BASE_GRAY,
                  vmin=0, vmax=1, interpolation="nearest", zorder=1)
        im = ax.imshow(maps[(mode, col)], extent=ext, origin="upper", cmap=CMAP,
                       vmin=-vmax, vmax=vmax, interpolation="nearest", zorder=2)
        decorate_map(ax, ext)
        ax.set_title(title, loc="left", fontweight="bold", pad=4)
        n_cities = int(means.loc[mode][col].notna().sum())
        ax.text(.99, .99, f"{n_cities} cities", transform=ax.transAxes, ha="right",
                va="top", fontsize=6.4, color="#4A5259", zorder=6)

    cax = fig.add_axes([.905, .20, .017, .55])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("City-mean level difference (K)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    cb.outline.set_linewidth(0.5)
    fig.text(.07, .99, "Same-source level differences by city, 2001–2025",
             ha="left", va="top", fontsize=11.5, weight="bold")
    # 脚注分两行：合成一行会超出画布右缘（2026-10-04 实测）。
    fig.text(.07, .030,
             "Blue = negative, red = positive; grey = mask with no estimable city-year; "
             "white = outside the 41 city masks.",
             ha="left", va="bottom", fontsize=6.5, color="#4A5259")
    fig.text(.07, .006,
             "A−B blue: independent-date estimate above same-day pairing. "
             "B−C red: composite below the daily estimate.",
             ha="left", va="bottom", fontsize=6.5, color="#4A5259")

    for e in ("png", "pdf"):
        fig.savefig(FIGS / f"fig4_00_city_level_maps.{e}",
                    dpi=400 if e == "png" else None)
    plt.close(fig)
    print("已写 fig4_00_city_level_maps.png/.pdf")


if __name__ == "__main__":
    main()
