"""Render Figure 5 — one worked example of why A−B is not zero.

Run after 09_build_daily_records.py and 16_build_revision_outputs.py. The city-year is
chosen by a rule written before looking at the figure (see pick_example), and every
number drawn is read from the frozen CSVs — none is typed in by hand.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "outputs" / "daily"
DATA = BASE / "data"
FIGS = BASE / "figures"
FIGS.mkdir(parents=True, exist_ok=True)
MIN_PX = 20

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 9,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "pdf.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

BLUE = "#0072B2"
ORANGE = "#D55E00"
GRAY = "#8C959D"


def pick_example(daily: pd.DataFrame, names: dict) -> tuple:
    """选例规则（先写死再看图）：

    夜间城市–年的 7 月，要求该月 ① 城乡同日有效日 ≥ 6；② 存在仅单侧有效的日期。
    在候选里优先南京、其次合肥；同一优先档内取**同日有效日最多**的年份（并列取较早年份）。
    "取观测最完整的 7 月"是事前可写的规则，不是事后挑一个差值最好看的案例。
    """
    jul = daily[daily.month == 7].copy()
    jul["paired"] = (jul.u_n >= MIN_PX) & (jul.r_n >= MIN_PX)
    jul["single_side"] = (jul.u_n >= MIN_PX) ^ (jul.r_n >= MIN_PX)
    g = jul.groupby(["city_id", "year"]).agg(n_paired=("paired", "sum"),
                                             n_single=("single_side", "sum"))
    ok = g[(g.n_paired >= 6) & (g.n_single >= 1)].reset_index()
    assert len(ok) > 0, "没有城市–年满足夜间 7 月选例规则"
    pref = {"Nanjing": 0, "Hefei": 1}
    ok["priority"] = ok.city_id.map(lambda c: pref.get(names.get(c, "").split("_")[-1], 2))
    ok = ok.sort_values(["priority", "n_paired", "year", "city_id"],
                        ascending=[True, False, True, True])
    row = ok.iloc[0]
    return int(row.city_id), int(row.year), ok


def main() -> None:
    cov = pd.read_csv(DATA / "city_covariates.csv")
    cov = cov[cov.year == 2020][["city_id", "city"]]
    names = dict(zip(cov.city_id, cov.city))

    daily = pd.read_csv(DATA / "daily" / "daily_records_night.csv",
                        usecols=["city_id", "year", "mode", "date", "month",
                                 "u_mean", "u_n", "r_mean", "r_n"])
    assert set(daily["mode"]) == {"night"}
    city_id, year, cand = pick_example(daily, names)
    print(f"选例：{names[city_id]}（city_id={city_id}）{year} 年 7 月")
    print(f"  候选数 {len(cand)}，前 5 名：")
    for _, r in cand.head(5).iterrows():
        print(f"    {names[int(r.city_id)]:22s} {int(r.year)} "
              f"同日有效性={int(r.n_paired)} 单侧日期={int(r.n_single)}")

    sel = daily[(daily.city_id == city_id) & (daily.year == year) &
                (daily.month == 7)].sort_values("date").copy()
    sel["paired"] = (sel.u_n >= MIN_PX) & (sel.r_n >= MIN_PX)
    sel["date"] = pd.to_datetime(sel["date"])
    sel["u_minus_r"] = sel.u_mean - sel.r_mean
    a_month = float(sel.loc[sel.paired, "u_minus_r"].mean())
    print(f"  7 月同日配对均值 A = {a_month:+.4f} K（{int(sel.paired.sum())} 个同日有效日）")

    # 交叉核对 1：该城市–月必须与冻结的月表逐值一致。
    mon = pd.read_csv(OUT / "same_source_decomposition_monthly.csv",
                      usecols=["basis", "year", "city_id", "month", "pixel_rule",
                               "min_daily_days", "min_8day_periods", "mode",
                               "A_paired_daily_K", "A_n_days"])
    mrow = mon[(mon.basis == "dynamic") & (mon.year == year) & (mon.city_id == city_id) &
               # 注意：列名 `mode` 与 DataFrame.mode() 同名，必须用 [] 取值，不能写 mon.mode。
               (mon["month"] == 7) & (mon["mode"] == "night") &
               (mon.pixel_rule == "20px") & (mon.min_daily_days == 6) &
               (mon.min_8day_periods == 2)]
    assert len(mrow) == 1, mrow
    a_frozen = float(mrow.A_paired_daily_K.iloc[0])
    assert abs(a_frozen - a_month) < 1e-4, (a_frozen, a_month)
    print(f"  冻结月表 A = {a_frozen:+.4f} K（{int(mrow.A_n_days.iloc[0])} 天）→ 一致")

    # 交叉核对 2：该城市–年必须存在于冻结的城市–年表。
    cy = pd.read_csv(OUT / "same_source_decomposition_city_year.csv",
                     usecols=["basis", "year", "city_id", "mode", "pixel_rule",
                              "min_daily_days", "min_8day_periods",
                              "A_paired_daily_common_months_mean_K",
                              "B_independent_dates_common_months_mean_K",
                              "delta_A_to_B_common_months_K"])
    cyrow = cy[(cy.basis == "dynamic") & (cy.year == year) & (cy.city_id == city_id) &
               (cy["mode"] == "night") & (cy.pixel_rule == "20px") &
               (cy.min_daily_days == 6) & (cy.min_8day_periods == 2)]
    assert len(cyrow) == 1, cyrow
    print(f"  冻结城市–年 A = {float(cyrow.A_paired_daily_common_months_mean_K.iloc[0]):+.4f} K · "
          f"B = {float(cyrow.B_independent_dates_common_months_mean_K.iloc[0]):+.4f} K · "
          f"A−B = {float(cyrow.delta_A_to_B_common_months_K.iloc[0]):+.4f} K")

    eight = pd.read_csv(OUT / "same_source_decomposition_8day_city_period.csv")
    eight = eight[(eight.basis == "dynamic") & (eight["mode"] == "night") &
                  (eight.year == year) & (eight.city_id == city_id)].copy()
    assert len(eight) > 0, "该城市–年无 8 天窗口记录"
    eight["date"] = pd.to_datetime(eight["date"])
    eight["C"] = eight.c_u_mean - eight.c_r_mean
    eight["D"] = eight.d_u_mean - eight.d_r_mean
    win = eight[(eight.date >= f"{year}-06-01") & (eight.date <= f"{year}-09-01")].sort_values("date")
    assert len(win) >= 3, len(win)

    fig, axes = plt.subplots(3, 1, figsize=(6.0, 7.6),
                             gridspec_kw={"hspace": .42})
    lab = f"{names[city_id].split('_')[-1]} {year}"

    ax = axes[0]
    ax.plot(sel.date, sel.u_n, "-o", color=ORANGE, ms=3.2, lw=1.0, label="Urban pixels")
    ax.plot(sel.date, sel.r_n, "-s", color=BLUE, ms=3.2, lw=1.0, label="Rural pixels")
    ax.axhline(MIN_PX, color="#59636D", ls="--", lw=0.9)
    ax.text(sel.date.iloc[-1], MIN_PX * 1.2, f"{MIN_PX}-pixel screen",
            ha="right", fontsize=7, color="#59636D")
    ax.set_yscale("log")
    ax.set_ylabel("Valid 1-km pixels")
    # 左下有数据点，图例改放右下并加白底，避免压线（首次出图时被压住）。
    ax.legend(frameon=True, facecolor="white", framealpha=.9, edgecolor="none",
              fontsize=7.3, ncol=2, loc="lower right")
    ax.set_title(f"(a)  Valid pixels per day, {lab} July", loc="left", fontweight="bold")
    ax.tick_params(axis="x", labelbottom=False)

    ax = axes[1]
    single = sel[~sel.paired]
    pair = sel[sel.paired]
    ax.plot(pair.date, pair.u_minus_r, "o", color="#216B83", ms=4.2,
            label="Both zones valid (paired day)")
    ax.plot(single.date, single.u_minus_r, "x", color=GRAY, ms=4.2,
            label="Only one zone valid")
    ax.axhline(a_month, color=ORANGE, lw=1.4)
    ax.text(sel.date.iloc[0], a_month, f"  A = {a_month:+.2f} K", va="bottom",
            fontsize=7.4, color=ORANGE)
    ax.axhline(0, color="#59636D", ls=":", lw=0.9)
    ax.set_ylabel("Daily U − R (K)")
    ax.legend(frameon=False, fontsize=7.3, loc="upper left")
    ax.set_title("(b)  Daily urban minus rural contrast", loc="left", fontweight="bold")
    ax.tick_params(axis="x", labelbottom=False)

    ax = axes[2]
    ax.plot(win.date, win.C, "-o", color=ORANGE, ms=3.4, lw=1.1,
            label="C: rebuilt from screened daily pixels")
    ax.plot(win.date, win.D, "-s", color="#6A3D9A", ms=3.4, lw=1.1,
            label="D: native MOD11A2")
    ax.axhline(0, color="#59636D", ls=":", lw=0.9)
    ax.set_ylabel("8-day U − R (K)")
    ax.set_xlabel("Date")
    ax.legend(frameon=False, fontsize=7.3, loc="best")
    ax.set_title("(c)  Rebuilt vs native 8-day windows", loc="left", fontweight="bold")

    for ax in axes:
        ax.grid(axis="y", color="#E6EAED", lw=0.6)
        ax.set_axisbelow(True)
    fig.suptitle("Why paired-date and independent-date contrasts differ",
                 x=.012, y=.995, ha="left", va="top", fontsize=10.5, weight="bold")
    fig.subplots_adjust(left=.135, right=.985, top=.935, bottom=.06)
    for e in ("png", "pdf"):
        fig.savefig(FIGS / f"fig5_00_mechanism_example.{e}",
                    dpi=400 if e == "png" else None)
    plt.close(fig)
    print("已写 fig5_00_mechanism_example.png/.pdf")


if __name__ == "__main__":
    main()
