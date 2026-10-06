"""QC 档位**流行率**：每档到底占多少像元（2026-10-02）。

## 为什么需要这个脚本

J9-b 定的主口径是 `bitwiseAnd(3).lte(1)`（优+次优），依据是用户提出的顾虑
「白天 QC 有『优』(位 00) 档，**夜间没有**」。这个顾虑当时只是口头判断。

2026-10-02 用 Drive 上遗留的 `qcprobe_night_2020_*.tif` 原始 `Q` 波段实测得：
夜间 `(Q&3)==0` 仅占 **0.235%**、`==1` 占 99.765%。但这属于**一次性内联命令**，
数字没有可复现的出处——违反本项目「数字须回溯到产出那一行」的纪律。
本脚本把那次测量**固化**：任何带 `Q` 波段的导出都能跑，输出落 `tables/`。

同时，一旦 `--with-qc` 的批次到盘，白天侧的同一张表可**同法**产出，
「昼夜 QC 是否对称」就从口头变成两侧都有数的对照。

## 它输出什么

- 每档（`(Q&3)==0` / `==1` / 其它）的像元数与占比（在**有效 LST 像元**内，`D != -9999`）
- 含该档像元的日数、按日的该档占比中位
- 该档像元数超过 `--day-threshold`（默认 1000）的日数

## 它**不**输出什么

不下任何效应或斜率结论。本脚本只回答「这档有多少数据」，不回答「用这档会不会翻号」
（后者是 `07_qc_tier_sensitivity.py`）。两者分开，因为流行率是**产品属性**，
与本文的城市掩膜、年份都无关。

用法：
    python 08_qc_prevalence.py                      # 自动找盘上带 Q 的导出
    python 08_qc_prevalence.py --mode night
    python 08_qc_prevalence.py --mode night --glob "qcprobe_night_2020_*.tif" --tag 2020probe

输出命名（2026-10-03 加）：默认全量运行写 `qc_prevalence_{mode}.csv`；
带 `--glob` 的探针运行自动写 `..._glob.csv`，或用 `--tag` 指定后缀。
**这是为了防止探针结果被全量运行冲掉**——2020 夜间那批（`46,204 / 19,669,507`）
只存在于探针文件里，全量运行不含它，若同名覆盖则该数字失去可复现出处。
"""

import argparse
import collections
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")   # 只改 stdout 时 stderr 仍走 GBK，中文报错会乱码

parser = argparse.ArgumentParser(description="QC 档位流行率")
parser.add_argument("--mode", choices=("night", "day"), default="night")
parser.add_argument("--glob", default=None,
                    help="显式指定导出文件通配（相对 data/gee_export/）；默认自动发现")
parser.add_argument("--day-threshold", type=int, default=1000,
                    help="「该档像元数 > 此值」才算这一天真的有这一档")
parser.add_argument("--tag", default=None,
                    help="输出文件名后缀。让**子集/探针**运行与默认全量运行并存。"
                         "没有它时两者写同一路径：探针结果（如 `qcprobe_night_2020_*` "
                         "那批）会被全量运行静默冲掉，正文与状态记录引用的数字就失去出处。"
                         "默认：给了 `--glob` 时为 `_glob`，否则无后缀。")
ARGS = parser.parse_args()
MODE = ARGS.mode
# 输出后缀。默认全量运行不带后缀；带 `--glob` 的探针运行自动加 `_glob`，
# 显式 `--tag` 优先（如 `--tag 2020probe`）。**统一补前导下划线**——否则
# `--tag 2020probe` 会产出 `..._night2020probe.csv`，与 `_glob` 的风格不一致，
# 也与手工命名（`..._night_2020probe.csv`）撞不上。
OSUF = ("_" + ARGS.tag.lstrip("_")) if ARGS.tag else ("_glob" if ARGS.glob else "")

FILL = -9999
# 与 03/07/08/16/16b 的第一条逐字相同；改动须全部一起改。
# 这里列出全部分档只为统计流行率，不涉及任何掩膜。
TIER_LABEL = {0: "优 (00)", 1: "次优 (01)"}

BASE = Path(__file__).resolve().parent.parent.parent      # analysis/
EXPORT = BASE / "data" / "gee_export"
TABDIR = BASE / "tables"
OUTDIR = BASE / "outputs" / "daily"
TABDIR.mkdir(parents=True, exist_ok=True)
OUTDIR.mkdir(parents=True, exist_ok=True)


def find_files():
    """优先用正式导出（`uhi_{tag}daily_*`），没有再退回探针文件（`qcprobe_*`）。

    正式文件可能**不带 Q 波段**（pilot 那批就不带），故不能只看文件名，
    还要验一下波段名——否则会静默统计出 0 个像元。
    """
    if ARGS.glob:
        return sorted(EXPORT.glob(ARGS.glob)), "显式 --glob"
    pat = f"uhi_{MODE}daily_*.tif"
    cand = sorted(EXPORT.glob(pat))
    if cand:
        return cand, f"正式导出 {pat}"
    probe = sorted(EXPORT.glob(f"qcprobe_{MODE}_*.tif"))
    if probe:
        return probe, f"探针文件 qcprobe_{MODE}_*.tif"
    return [], "（无）"


files, how = find_files()
print("=" * 78)
print(f"QC 档位流行率   [mode={MODE}]")
print("=" * 78)
print(f"  文件来源：{how}   共 {len(files)} 个")
if not files:
    print("data/gee_export/ 下没有可用文件（需带 Q 波段的导出，或先跑 qcprobe）",
          file=sys.stderr)
    raise SystemExit(1)

tot = collections.Counter()          # 档 -> 像元数（有效 LST 内）
n_valid = 0
per_day = []                         # (日期, 有效像元数, {档: 数})
n_bands_seen = 0

for f in files:
    with rasterio.open(f) as s:
        desc = [d for d in s.descriptions if d]
        idx = {d: i + 1 for i, d in enumerate(s.descriptions) if d}
        pairs = [(d, f"Q{d[1:]}") for d in desc
                 if d.startswith("D") and f"Q{d[1:]}" in idx]
        if not pairs:
            continue
        # **一次性读全部波段，不要逐波段 `s.read(band)`。** 这些导出是
        # `interleave=pixel`：一个 256×256 瓦片里交错存着全部 332 个波段的像素，
        # 于是「读第 k 个波段」实际要把整块瓦片解压、再丢掉 331/332——
        # 单波段读的代价约等于整幅文件。
        # 2026-10-03 实测（`qcprobe_night_2020_h1/h2` 两个同构文件，交叉测避免
        # 页缓存污染）：逐波段 193 s/文件，一次性读 1.8 s/文件，**108×**。
        # 这才是 08 全量跑不完的真因——按旧路径昼夜合计约 5.7 h，两次都在
        # 30 min 后台时限被杀。代价是峰值内存 ≈ 波段数 × 756,675 × 2 B（约 0.5 GB）。
        arr = s.read()
        for d, qn in pairs:
            n_bands_seen += 1
            lst = arr[idx[d] - 1]
            ok = lst != FILL                     # 有效 LST 像元 = 非填充
            q = arr[idx[qn] - 1][ok]
            # `np.bincount` 而非 `Counter((q & 3).tolist())`：后者逐像元建 Python
            # int 再哈希，单波段最高 75 万次。同批实测省约 1 s/波段（非主因）。
            cnt = np.bincount((q & 3).astype(np.int8), minlength=4)
            nz = np.nonzero(cnt)[0]
            for k in nz:
                tot[int(k)] += int(cnt[k])
            n_valid += int(ok.sum())
            per_day.append((d[1:], int(ok.sum()),
                            {int(k): int(cnt[k]) for k in nz}))
        del arr          # 下一个文件前先释放 ~0.5 GB，别让峰值叠加

if n_bands_seen == 0:
    print(f"{len(files)} 个文件里没有任何 D/Q 配对波段——"
          "这批导出大概没带 --with-qc", file=sys.stderr)
    raise SystemExit(1)

print(f"  配对波段 {n_bands_seen} 个 · 有效 LST 像元 {n_valid:,}\n")

rows = []
for tier in sorted(tot):
    n = tot[tier]
    shares = [c.get(tier, 0) / v for _, v, c in per_day if v > 0]
    days_any = sum(1 for _, _, c in per_day if c.get(tier, 0) > 0)
    days_big = sum(1 for _, _, c in per_day if c.get(tier, 0) > ARGS.day_threshold)
    rows.append({
        "mode": MODE,
        "tier": tier,
        "label": TIER_LABEL.get(tier, f"其它 ({tier:02b})"),
        "pixels": n,
        "share_of_valid": round(n / n_valid, 6) if n_valid else np.nan,
        "days_with_any": days_any,
        "days_total": len(per_day),
        "days_over_threshold": days_big,
        "median_daily_share": round(float(np.median(shares)), 6) if shares else np.nan,
    })

df = pd.DataFrame(rows)
out = TABDIR / f"qc_prevalence_{MODE}{OSUF}.csv"
df.to_csv(out, index=False, encoding="utf-8-sig")
print(df.to_string(index=False))
print(f"\n已写 {out.relative_to(BASE)}")

# ---------------------------------------------------------------- 判读
lines = [f"# QC 档位流行率（{MODE}）", "",
         f"来源：`data/gee_export/` 的 {len(files)} 个文件（{how}），"
         f"{n_bands_seen} 个 D/Q 配对波段，**有效 LST** 像元 {n_valid:,} 个。", "",
         "「有效 LST」= 非填充（`D != -9999`）。占比的分母是**有效 LST 像元**，",
         "不是全部像元——故本表回答「拿到一个有效观测，它是哪一档」。", "",
         df.to_markdown(index=False), ""]

good = df[df.tier == 0]
other = df[df.tier == 1]
if len(good) and n_valid:
    g = good.iloc[0]
    lines += ["## 读法", "",
              f"- 「优」档 `(Q&3)==0` 占 **{g.share_of_valid:.3%}**"
              f"（{int(g.pixels):,} / {n_valid:,}），按日占比中位 "
              f"**{g.median_daily_share:.3%}**；{int(g.days_total)} 天里有 "
              f"{int(g.days_over_threshold)} 天该档超过 {ARGS.day_threshold} 像元。"]
    if len(other):
        o = other.iloc[0]
        lines.append(f"- 「次优」档 `(Q&3)==1` 占 **{o.share_of_valid:.3%}**"
                     f"（{int(o.pixels):,}）。")
    if g.share_of_valid < 0.01:
        lines += ["",
                  "⇒ **该档在此模式下几乎为空**。若正文声称「只用优档」，"
                  "实际上要么用掉的数据极少、要么是误读了 QC 位。",
                  "  这正是 `07_qc_tier_sensitivity.py` 的 `eq0` 收严档在此模式"
                  "**结构性不可执行**的原因——不是脚本缺陷，是产品属性。"]
    else:
        lines += ["",
                  "⇒ 该档有实质数据，`07_qc_tier_sensitivity.py` 的 `eq0` 收严档"
                  "**可执行**，其结论可以报告。"]
lines += ["", "> 本脚本只报**有多少数据**，不报**效应会不会翻号**——后者见 `07`。",
          "> 流行率是产品属性，与本文的城市掩膜、年份选择无关，故可单独引用。", ""]

rep = OUTDIR / f"QC_PREVALENCE_{MODE.upper()}{OSUF}.md"
rep.write_text("\n".join(lines), encoding="utf-8")
print(f"已写 {rep.relative_to(BASE)}")
