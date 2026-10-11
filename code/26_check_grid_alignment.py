"""Quadrature gate: Aqua rasters must sit on *exactly* the Terra grid (2026-10-11).

Why this exists
---------------
The whole Aqua replication rests on one assumption: `MYD11A1/MYD11A2` were exported
with the same `region` / `scale` / `crs` as the Terra batch, so the **same local
mask** (`masks_matched.npz`) can be applied pixel-for-pixel without re-deriving
anything. If that assumption is off by even one pixel, every zone mean shifts and
**nothing raises an error** — the pipeline happily produces numbers that are wrong.

So this gate compares, per file:

1. `width` / `height` / `transform` / `crs` against the Terra file for the same
   mode+year(+half);
2. the **band descriptions** (`D{YYYYMMDD}`) — reported as *coverage*, not pass/fail:
   the two products genuinely carry different day sets (Terra MOD11A1 is missing
   granules that Aqua has — measured 2026-10-11: 7 days in 2003 H2, 2 in 2008 H2).
   Grids must coincide exactly; band order is irrelevant because the pipeline matches
   on the description. Only a *large* shortfall (see DATE_OVERLAP_MIN) fails.
3. `width`/`height` against the mask's own `H`/`W`, which is the property the local
   masking actually depends on.

Usage
-----
    python 26_check_grid_alignment.py                # aqua vs terra, both modes
    python 26_check_grid_alignment.py --modes day

Exit code 0 = every checked file aligned; 1 = at least one mismatch (listed).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import rasterio

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

BASE = Path(__file__).resolve().parents[2]          # analysis/
EXPORT = BASE / "data" / "gee_export"
DATA = BASE / "data"

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from products import get as get_product  # noqa: E402


def d_bands(path: Path) -> list[str]:
    """`D*` band descriptions, in file order (Q bands are not comparable across products)."""
    with rasterio.open(path) as src:
        return [d for d in src.descriptions if d and d.startswith("D")]


DATE_OVERLAP_MIN = 0.90        # 两平台日期集重叠率低于此才判 FAIL（差异按覆盖报出）
coverage_notes: list[str] = []


def compare_dates(a_name: str, da: list[str], t_name: str, dt: list[str]) -> bool:
    """返回 True=可接受。日期集不同是**预期**（Terra 有缺granule），只在大面积缺失时失败。"""
    sa, st = set(da), set(dt)
    inter, union = len(sa & st), len(sa | st)
    overlap = inter / union if union else 1.0
    if sa != st:
        only_t = sorted(st - sa)
        only_a = sorted(sa - st)
        coverage_notes.append(
            f"{a_name} vs {t_name}: 重叠 {inter}/{union} = {overlap:.1%}"
            f" · 仅 Terra {len(only_t)} 个{(' ' + str(only_t[:3])) if only_t else ''}"
            f" · 仅 Aqua {len(only_a)} 个{(' ' + str(only_a[:3])) if only_a else ''}")
    return overlap >= DATE_OVERLAP_MIN


def grid(path: Path):
    with rasterio.open(path) as src:
        return {"w": src.width, "h": src.height, "t": src.transform, "crs": str(src.crs)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--product", choices=("aqua",), default="aqua",
                    help="只有 aqua 需要本闸门（terra 是基准本身）")
    ap.add_argument("--modes", nargs="+", choices=("day", "night"), default=["day", "night"])
    args = ap.parse_args()

    prod = get_product(args.product)
    terra = get_product("terra")
    masks = np.load(DATA / "masks_matched.npz")
    MH, MW = int(masks["H"]), int(masks["W"])
    print(f"掩膜网格 masks_matched.npz: H={MH} W={MW}")

    problems: list[str] = []
    checked = 0

    for mode in args.modes:
        a_stem = prod["stem_daily"].format(tag=mode)
        t_stem = terra["stem_daily"].format(tag=mode)
        years = sorted({int(p.name.split("_")[2]) for p in EXPORT.glob(f"{a_stem}_*.tif")})
        print(f"\n[{mode}] Aqua 逐日文件年份 {len(years)} 个："
              f"{years[0] if years else '-'}–{years[-1] if years else '-'}")
        for y in years:
            for half in ("h1", "h2"):
                a = EXPORT / f"{a_stem}_{y}_{half}.tif"
                t = EXPORT / f"{t_stem}_{y}_{half}.tif"
                if not a.exists():
                    problems.append(f"缺文件 {a.name}")
                    continue
                ga = grid(a)
                checked += 1
                if (ga["w"], ga["h"]) != (MW, MH):
                    problems.append(f"{a.name}: 网格 {ga['w']}x{ga['h']} != 掩膜 {MW}x{MH}")
                if not t.exists():
                    print(f"  {y} {half}: ⚠ 无 Terra 对照（{t.name} 不在盘上），只做了掩膜尺寸检查")
                    continue
                gt = grid(t)
                if (ga["w"], ga["h"], ga["t"], ga["crs"]) != (gt["w"], gt["h"], gt["t"], gt["crs"]):
                    problems.append(f"{a.name}: 网格/仿射/CRS 与 {t.name} 不等")
                da, dt_ = d_bands(a), d_bands(t)
                if not compare_dates(a.name, da, t.name, dt_):
                    problems.append(f"{a.name}: 与 {t.name} 日期集重叠 <{DATE_OVERLAP_MIN:.0%}"
                                    f"（aqua {len(da)} / terra {len(dt_)}）")
            if (y - years[0]) % 5 == 0:
                print(f"  … 已查至 {y}")

        # 8 天原生序列（估计量 D）
        a_seq = prod["stem_8day"].format(tag=mode)
        t_seq = terra["stem_8day"].format(tag=mode)
        for y in years:
            a = EXPORT / f"{a_seq}_{y}.tif"
            t = EXPORT / f"{t_seq}_{y}.tif"
            if not a.exists():
                problems.append(f"缺文件 {a.name}")
                continue
            ga = grid(a)
            checked += 1
            if (ga["w"], ga["h"]) != (MW, MH):
                problems.append(f"{a.name}: 网格 {ga['w']}x{ga['h']} != 掩膜 {MW}x{MH}")
            if not t.exists():
                continue
            gt = grid(t)
            if (ga["w"], ga["h"], ga["t"], ga["crs"]) != (gt["w"], gt["h"], gt["t"], gt["crs"]):
                problems.append(f"{a.name}: 网格/仿射/CRS 与 {t.name} 不等")
            da, dt_ = d_bands(a), d_bands(t)
            if not compare_dates(a.name, da, t.name, dt_):
                problems.append(f"{a.name}: 与 {t.name} 合成期重叠 <{DATE_OVERLAP_MIN:.0%}"
                                f"（aqua {len(da)} / terra {len(dt_)}）")

    print("\n" + "=" * 78)
    print(f"检查 {checked} 个 Aqua 栅格")
    if coverage_notes:
        print(f"\n覆盖差异（**预期**，非缺陷；Terra MOD11A1 本身有缺 granule 的日期）：")
        for n in coverage_notes:
            print(f"    {n}")
    if problems:
        print(f"\n✗ {len(problems)} 项**网格/结构**不一致 —— **停止**，先查清再进分析：")
        for p in problems[:40]:
            print(f"    {p}")
        return 1
    print("\n✓ 网格/仿射/CRS 与 Terra 逐值相等，且等于掩膜尺寸；日期集差异已如实列出（非缺陷）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
