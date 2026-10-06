"""QC 档位敏感性：`lte(1)`（主口径）vs `eq(0)`（收严）——昼夜分开做（2026-10-02）。

## 为什么单独一个脚本

主口径 `03_build_daily_panel.py` 只读 `D*` 波段（已按 `bitwiseAnd(3).lte(1)` 掩膜）。
收严档要**另起一套有效像元判定**，无法从 `D*` 反推，故 02 导出时额外带了原始
`Q*` 波段（每景波段数 ×2）。本脚本读 `Q*`，在**同一份栅格**上重算，掩膜仍取
`masks_matched.npz`（12 的唯一产物）——不引入第二份城乡判定。

## 它检验什么（用户的原始顾虑）

> 白天 QC 有「优」（位 00）档，夜间没有。**必须昼夜一起换**。

主口径是 `lte(1)`（优+次优）。本脚本给出 `eq(0)`（**只要优**）下的同一组数，
回答两件事：

1. 收严后**还剩多少有效像元/日**——夜间若几乎归零，则「严档同号」这条稳健性
   检验在夜间**结构上不可执行**，正文必须写明，而不是假装做过。
2. 收严后**加权效应的符号与量级**是否维持。

## 与主口径的可比性边界

收严集是主口径集的**子集**（`(Q&3)==0` ⊂ `(Q&3)<=1`），故两者的差只来自掩膜，
不来自任何其他口径差异。但**保留的像元不是随机子集**——它是算法更自信的那部分，
本身可能与城乡温差相关。故本脚本的结论只能读作「换档会不会翻号」，
**不能**读作「主口径偏了多少」。

用法：
    python 07_qc_tier_sensitivity.py --mode night --years 2020 2023
    python 07_qc_tier_sensitivity.py --mode day --years 2020 2023
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")   # 只改 stdout 时 stderr 仍走 GBK，中文报错会乱码

parser = argparse.ArgumentParser(description="QC 档位敏感性（逐日）")
parser.add_argument("--mode", choices=("night", "day"), default="night")
parser.add_argument("--years", type=int, nargs="+", default=[2020, 2023])
parser.add_argument("--glob", default=None,
                    help="显式指定导出文件通配（相对 data/gee_export/），且**全部视为 "
                         "`--years` 里的年份**；默认按 uhi_{mode}daily_{year}_*.tif 逐年找。"
                         "用于在正式批次到盘前拿探针文件（qcprobe_*）先跑通。")
parser.add_argument("--tag", default=None,
                    help="输出文件名后缀。让**子集/探针**运行与默认全量运行并存。"
                         "没有它时两者写同一路径：探针结果（如 2020 夜间那批）会被"
                         "全量运行静默冲掉，正文与状态记录引用的数字就失去出处。"
                         "默认：给了 `--glob` 时为 `_glob`，否则无后缀。")
ARGS = parser.parse_args()
MODE, TAG = ARGS.mode, ARGS.mode
# 输出后缀。**不能复用 `TAG`**——那是模式名，用于拼 `uhi_{TAG}daily_...`。
# 统一补前导下划线，使 `--tag x` 与自动的 `_glob` 命名风格一致。
OSUF = ("_" + ARGS.tag.lstrip("_")) if ARGS.tag else ("_glob" if ARGS.glob else "")
LABEL = "夜间" if MODE == "night" else "白天"

FILL = -9999
SCALE_K = 0.02
MIN_DAYS = 6
# 主口径 = bitwiseAnd(3).lte(1)（优+次优）；收严 = bitwiseAnd(3).eq(0)（只要优）。
# 与 08/16/16b/02 的第一条逐字相同，改动须四处一起改。
TIERS = [("lte1", lambda q: (q & 3) <= 1), ("eq0", lambda q: (q & 3) == 0)]

BASE = Path(__file__).resolve().parent.parent.parent      # analysis/
DATA = BASE / "data"
EXPORT = DATA / "gee_export"
OUTDIR = BASE / "outputs" / "daily"
TABDIR = BASE / "tables"
OUTDIR.mkdir(parents=True, exist_ok=True)
TABDIR.mkdir(parents=True, exist_ok=True)

masks = np.load(DATA / "masks_matched.npz")

print("=" * 78)
print(f"{LABEL} QC 档位敏感性   [mode={MODE} · 年 {ARGS.years}]")
print("=" * 78)
print("  主口径 lte(1)（优+次优） vs 收严 eq(0)（只要优）")
print(f"  掩膜 {DATA.name}/masks_matched.npz（12 的唯一产物）\n")


if ARGS.glob and len(ARGS.years) > 1:
    # `--glob` 命中的文件不带年份，若指定多年会被**逐个年份重复读入**，
    # 把同一份数据当成两个年份的观测——静默得出假结果。故直接拒绝。
    print(f"--glob 只能配单一年份（现给 {ARGS.years}）：所匹配的文件不含年份，"
          f"多给会被重复计成多个年份", file=sys.stderr)
    raise SystemExit(1)


def files_for(year):
    """该年份要读哪些文件。`--glob` 时全部文件都算给该年份（探针文件不含年份）。"""
    if ARGS.glob:
        return sorted(EXPORT.glob(ARGS.glob))
    return sorted(EXPORT.glob(f"uhi_{TAG}daily_{year}_*.tif"))


def load_pairs(year):
    """→ {日期: (D 波段号, Q 波段号)}，按日期序。Q 缺失的日期直接丢弃并计数。"""
    files = files_for(year)
    got = {}
    for f in files:
        with rasterio.open(f) as s:
            desc = list(s.descriptions)
            idx = {d: i + 1 for i, d in enumerate(desc) if d}
            for d in desc:
                if d and d.startswith("D") and f"Q{d[1:]}" in idx:
                    date = d[1:]
                    if date in got:
                        raise RuntimeError(f"{year}: 日期 {date} 在两个文件里重复出现")
                    got[date] = (idx[d], idx[f"Q{d[1:]}"])
    return dict(sorted(got.items()))


rows = []
for year in ARGS.years:
    pairs = load_pairs(year)
    if not pairs:
        print(f"{year}: 缺带 QC 的导出文件，跳过（02 须带 --with-qc 重导）")
        continue
    cities = sorted(int(k.split("|")[1]) for k in masks.files
                    if k.startswith(f"u{year}|") and f"r{year}|{k.split('|')[1]}" in masks.files)
    u_px = {c: len(masks[f"u{year}|{c}"]) for c in cities}
    r_px = {c: len(masks[f"r{year}|{c}"]) for c in cities}

    # 「像元 → 城市序号」查找表（0 = 不属于任何城），用途与 03 同：把逐城
    # fancy-indexing（2 档 × 2 层 × nc 城 = 144 次/波段）换成 `bincount` 一次遍历。
    # 未改前跑全量 25 年是小时级。口径不变：仍是非填充且过该档 QC 的像元均值。
    nc = len(cities)
    ci_u = np.zeros(int(masks["H"]) * int(masks["W"]), dtype=np.int32)
    ci_r = np.zeros(int(masks["H"]) * int(masks["W"]), dtype=np.int32)
    for j, c in enumerate(cities, start=1):
        ci_u[masks[f"u{year}|{c}"]] = j
        ci_r[masks[f"r{year}|{c}"]] = j

    # acc[tier][city] = {"per_day": [(um, un, rm, rn)], "u_sum","u_n","r_sum","r_n"}
    acc = {t: {c: {"per_day": [], "u_sum": 0.0, "u_n": 0, "r_sum": 0.0, "r_n": 0}
               for c in cities} for t, _ in TIERS}

    for f in files_for(year):
        with rasterio.open(f) as src:
            desc = list(src.descriptions)
            idx = {d: i + 1 for i, d in enumerate(desc) if d}
            pairs = [(d, f"Q{d[1:]}") for d in desc
                     if d and d.startswith("D") and f"Q{d[1:]}" in idx]
            if not pairs:
                continue
            # **一次性读全部波段**，不要逐波段 `src.read(band)`——理由与 08 逐字
            # 相同：导出是 `interleave=pixel`，读单个波段要把整块 256×256 瓦片
            # 解压后丢掉 331/332，单波段读的代价≈整幅文件。2026-10-03 实测
            # **108×**（193 s/文件 → 1.8 s/文件）。不修则本脚本全量 23 年是小时级。
            arr = src.read()
            for d, qn in pairs:
                lst = arr[idx[d] - 1].astype(np.float32)
                lst[lst == FILL] = np.nan
                lst *= SCALE_K
                qc = arr[idx[qn] - 1]
                flat_lst, flat_qc = lst.ravel(), qc.ravel()
                # `lst` 已把 FILL 置 nan，故 `isfinite` 即「非填充」。原写法把两步
                # （QC 筛 + 非填充筛）分开做，此处合并——结果集合相同，少两次全幅拷贝。
                lst_ok = np.isfinite(flat_lst)
                for tier, keep in TIERS:
                    ok = lst_ok & keep(flat_qc)
                    fv = flat_lst[ok]
                    iu, ir = ci_u[ok], ci_r[ok]
                    cnt_u = np.bincount(iu, minlength=nc + 1)
                    sum_u = np.bincount(iu, weights=fv, minlength=nc + 1)
                    cnt_r = np.bincount(ir, minlength=nc + 1)
                    sum_r = np.bincount(ir, weights=fv, minlength=nc + 1)
                    for j, c in enumerate(cities, start=1):
                        un, rn = int(cnt_u[j]), int(cnt_r[j])
                        a = acc[tier][c]
                        a["u_sum"] += float(sum_u[j]); a["u_n"] += un
                        a["r_sum"] += float(sum_r[j]); a["r_n"] += rn
                        a["per_day"].append(
                            (float(sum_u[j] / un) if un else np.nan, un,
                             float(sum_r[j] / rn) if rn else np.nan, rn))
            del arr          # 下一个文件前释放 ~0.5 GB，别让峰值叠加

    base_t = TIERS[0][0]
    for tier, _ in TIERS:
        pw_l, eq_l, noc = [], [], 0
        day_ret = []                           # 可算日数中位（受门限约束，属效应侧）
        for c in cities:
            a = acc[tier][c]
            if a["u_n"] == 0 or a["r_n"] == 0:
                noc += 1
                continue
            pw = a["u_sum"] / a["u_n"] - a["r_sum"] / a["r_n"]
            um_all = [um for um, un, _, _ in a["per_day"] if un > 0]
            rm_all = [rm for _, _, rm, rn in a["per_day"] if rn > 0]
            if len(um_all) < MIN_DAYS or len(rm_all) < MIN_DAYS:
                noc += 1
                continue
            eq = float(np.mean(um_all)) - float(np.mean(rm_all))
            day_ret.append(len(um_all))
            pw_l.append(pw); eq_l.append(eq)

        # ---- 数据量（**独立于 MIN_DAYS 门限**）--------------------------------
        # 「这档还剩多少数据」与「数据够不够估一个效应」是两个问题。把数据量挂在
        # 门限之下会出现最坏的组合：收严档越没数据（正是最该看数据量的场合），
        # 数字越报 NaN——2026-10-02 实跑夜间 eq0 就是这个症状（n_cities=0，
        # 保留率全 NaN，等于该档的关键信息一条都没留下）。
        #
        # `u_n` 在 `for f in files`→`for d in desc`（每天）内累加，是**像元·日**总数，
        # 故 u_n/u_px 是「每像元平均有效日数」（量级 10–60），**不是**保留率。
        # 保留率只能在**同一批像元**上比两档的像元·日（相除时 u_px 约掉）。
        u_dpp, r_dpp, u_share, r_share = [], [], [], []
        for c in cities:
            a = acc[tier][c]
            if a["u_n"] > 0:
                u_dpp.append(a["u_n"] / u_px[c])
            if a["r_n"] > 0:
                r_dpp.append(a["r_n"] / r_px[c])
            if tier != base_t:
                b_u, b_r = acc[base_t][c]["u_n"], acc[base_t][c]["r_n"]
                if b_u and b_r:
                    u_share.append(a["u_n"] / b_u)
                    r_share.append(a["r_n"] / b_r)

        # 绝对像元·日。占比（下两列）是**逐城比值的均值**，在收严档会退化成
        # 「极小数 ÷ 大数」——城区 0.0% 与乡村 0.3% 并排时读者容易读成
        # 「优档偏乡村」，但那可能只是 3 : 300 这种量级。故把分子分母都留在表里，
        # 让「差距是真实的还是噪声」可以直接判断，不靠脑补。
        tot_u = sum(acc[tier][c]["u_n"] for c in cities)
        tot_r = sum(acc[tier][c]["r_n"] for c in cities)

        rows.append({
            "mode": MODE, "year": year, "tier": tier, "n_cities": len(pw_l),
            "dropped_cities": noc,
            "urban_pixel_days": tot_u, "rural_pixel_days": tot_r,
            "pixel_weighted_K": round(float(np.mean(pw_l)), 3) if pw_l else np.nan,
            "equal_any_K": round(float(np.mean(eq_l)), 3) if eq_l else np.nan,
            "weighting_effect_K": round(float(np.mean(eq_l) - np.mean(pw_l)), 3) if pw_l else np.nan,
            "urban_days_per_px": round(float(np.mean(u_dpp)), 1) if u_dpp else np.nan,
            "rural_days_per_px": round(float(np.mean(r_dpp)), 1) if r_dpp else np.nan,
            "urban_retention_vs_main": round(float(np.mean(u_share)), 3) if u_share else np.nan,
            "rural_retention_vs_main": round(float(np.mean(r_share)), 3) if r_share else np.nan,
            "median_urban_days": float(np.median(day_ret)) if day_ret else np.nan,
        })

if not rows:
    print("没有可算的档位（导出文件缺 Q 波段？）", file=sys.stderr)
    raise SystemExit(1)

df = pd.DataFrame(rows)
out = TABDIR / f"qc_tier_sensitivity_{MODE}{OSUF}.csv"
df.to_csv(out, index=False, encoding="utf-8-sig")
print(df.to_string(index=False))
print(f"\n已写 {out.relative_to(BASE)}")

# ---------------------------------------------------------------- 判读
lines = [f"# QC 档位敏感性（逐日 {LABEL}）", "",
         "主口径 `bitwiseAnd(3).lte(1)`（优+次优）vs 收严 `bitwiseAnd(3).eq(0)`（只要优）。",
         "掩膜同为 `masks_matched.npz`，两档只差 QC 判定。", "", df.to_markdown(index=False), ""]
for year in ARGS.years:
    sub = df[df.year == year]
    if len(sub) < 2:
        continue
    base = sub[sub.tier == "lte1"].iloc[0]
    strict = sub[sub.tier == "eq0"].iloc[0]
    # 两件事分开报：**数据量**（总可报）与**能否估效应**（可能为 0 城）。
    lines += [f"## {year}", "",
              f"- 主口径 `lte1` 效应：{base.n_cities} 城可估 · 加权 "
              f"{base.pixel_weighted_K:+.3f} · 等权 {base.equal_any_K:+.3f} · "
              f"**效应 {base.weighting_effect_K:+.3f} K**",
              f"- 收严 `eq0` 数据量：像元·日 城区 "
              f"{int(strict.urban_pixel_days):,} · 乡村 {int(strict.rural_pixel_days):,}"
              f"（对主口径之比 城区 {strict.urban_retention_vs_main:.4f} · 乡村 "
              f"{strict.rural_retention_vs_main:.4f}）；每像元有效日数 城区 "
              f"{strict.urban_days_per_px:.3f} · 乡村 {strict.rural_days_per_px:.3f}",
              "  （城乡两列的差距**不要**读作「优档偏乡村」的分布证据：收严档在"
              "夜间的总量本身极小——见同目录 `QC_PREVALENCE_*.md`——逐城比值不稳，"
              "此处只报数，不作分布判断。）"]
    if strict.n_cities == 0:
        lines += [f"- 收严 `eq0` **一个城都估不出效应**（{int(strict.dropped_cities)} 城"
                  "因有效日数不足被剔除）。**这不是脚本缺陷，是数据量本身不够**——"
                  "占比见同目录 `QC_PREVALENCE_*.md`。", ""]
    else:
        lines.append(f"- 收严 `eq0` 效应：{strict.n_cities} 城可估 · 加权 "
                     f"{strict.pixel_weighted_K:+.3f} · 等权 "
                     f"{strict.equal_any_K:+.3f} · "
                     f"**效应 {strict.weighting_effect_K:+.3f} K**")
        same = (base.weighting_effect_K * strict.weighting_effect_K) > 0
        lines += [f"→ {'**同号**' if same else '**翻号**'}；"
                  + ("收严档仍有数据，稳健性检验**可执行**。"
                     if strict.n_cities >= 0.8 * base.n_cities
                     else f"收严后仅剩 {strict.n_cities}/{base.n_cities} 城，"
                          "**该检验在此档位不可执行**，正文须如实写明。"), ""]

rep = OUTDIR / f"QC_TIER_{MODE.upper()}{OSUF}.md"
rep.write_text("\n".join(lines), encoding="utf-8")
print(f"已写 {rep.relative_to(BASE)}")
