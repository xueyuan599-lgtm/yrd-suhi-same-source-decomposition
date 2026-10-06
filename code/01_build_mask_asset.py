"""把 `12_build_panel.py` 的城乡掩膜送上 Earth Engine —— **上传段被阻断（2026-10-02 实测）**。

> ⚠ **状态：掩膜 GeoTIFF 的生成与回读校验已通过；`earthengine upload image` 段不可用。**
>
> 实测（本机、本日）：
> 1. 该 CLI 版本**只收 Cloud Storage 源**：传本地路径直接报
>    `Invalid Cloud Storage URL`（`ee/cli/utils.py` 的 `is_gcs_path` 只认 `gs://`）。
> 2. 旧版 CLI 的「本地文件自动暂存到 EE 桶」能力**已移除**；
>    `upload_dir_to_bucket()` 也是 `get_bucket(dest)` —— 目标桶必须已存在。
> 3. 本项目 **没有任何 GCS 桶**（用 EE 持久凭据列桶返回空；凭据本身带
>    `devstorage.full_control` 作用域，是**没有桶**而非没权限）。建桶需计费项目，
>    属对外动作，不在本任务授权内。
>
> 结论：**日下无法创建 asset**。逐日 pilot 因此改走「只对 2020+2023 导逐日栅格、
> 本地用 12 的原掩膜归约」的 B 方案（见 `02_export_daily_seq.py`）。
> 本脚本**保留**：第一段（掩膜 → 与导出网格逐字同源的 GeoTIFF + 回读校验）是对的，
> 日后若有可用桶，补上上传即可复用。
>
> **后续更正（2026-10-02 实测）**：全量栅格体量此前估为 ~28 GB，实测两年 454.5 MB、
> 外推 25 年约 **5.7 GB**，该估算**作废**。故 B 方案本就能铺满 2001–2025，
> asset 从「必须解决」降为**可选优化**（省的是本地磁盘与一遍本地归约，不是可行性）。
> 若日后仍要建 asset：先建桶，上传段可直接复用；或把掩膜在 EE 侧用**已验证与 12
> 逐像元相同**的方式重建（须逐像元比对后方可采信）。

## 为什么必须上传，不能在 EE 里重算

逐日方案走**服务器端 zonal reduce**（不导栅格）。服务器要按城乡分区求均值，
就必须知道**哪些像元属于哪个城市的城区 / 乡村**。

规则（`.claude/rules/…` 与 DESIGN.md §9.1.2）写明：**城乡掩膜的唯一实现是
`12_build_panel.py`**。若在 EE 里照 GHSL/NDVI/SRTM 再算一遍，那就是第二份实现——
两份迟早漂移，而口径不一致是**静默错误**（两边都跑得通，数却对不上）。

故本脚本只做一件事：把 `masks_matched.npz` 里 12 已经算好的下标，
**原样**写成 GeoTIFF 并上传成 asset。掩膜的作者仍然只有 12 一个。

## 网格必须逐字相同（本脚本最容易错的地方）

已实测：导出 tif 的变换是

    0.008983152841195215, 0, 114.840625921839, 0, -0.008983152841195215, 35.106161303

而 `x0/px = 12784.0`、`y0/px = -3908.0` —— **都是整数**。
原点正好落在像元尺寸的整数倍上，说明 EE 吸附的是 **CRS 全局原点栅格**
（倍数格），**不是**本次导出 region 的边界。

> **2026-10-02 更正**：此处原写「`x0/s = 12784.004`、都不是整数 / 按 region 边界吸附」，
> 该数**回溯不到任何产出**，与 `zones_manifest.json` 记录的 `crsTransform` 自算结果
> （恰为整数）矛盾，已按实测改写。工程结论未变，但理由换了。

后果与对策：既然吸附基准是**全局格**，只要 scale 与 CRS 相同，EE 的多次导出
与本脚本用 rasterio 写出的掩膜就落在**同一张格**上。本脚本的职责是把
`--reference` 的 transform **原样**记进 `data/daily/zones_manifest.json`
的 `crsTransform`，下游据此对齐，**不各写一份常数**——手工写常数才是错位的真正来源。

## 编码

每个城市一个波段，取值：1 = 城区 · 2 = 乡村 · 0 = 其他。
波段名写成 `c{city_id}`，但下游**一律按波段号 select**——
EE 对 GeoTIFF 波段描述的保留不保证，按名字选会静默选错。

用法：
    python 01_build_mask_asset.py --dry-run   # 只写 tif + 自检，不上传
    python 01_build_mask_asset.py             # 写 tif 并上传 asset
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import rasterio

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")   # 只改 stdout 时 stderr 仍走 GBK，中文报错会乱码

PROJECT = "evident-ocean-510002-t0"
PILOT_YEARS = (2020, 2023)
ASSET_PREFIX = "yrd_zones"
# 网格基准：与 08/16/16b 导出逐字同源的参考文件（只用来取 transform / crs / 尺寸）
REFERENCE = "uhi_cityid.tif"
EE_CLI = Path(__file__).resolve().parents[2] / ".venv_daily" / "Scripts" / "earthengine.exe"

parser = argparse.ArgumentParser(description="城乡掩膜 → EE asset")
parser.add_argument("--dry-run", action="store_true", help="只写 tif 与自检，不上传")
parser.add_argument("--years", type=int, nargs="+", default=list(PILOT_YEARS))
ARGS = parser.parse_args()

BASE = Path(__file__).resolve().parent.parent.parent      # analysis/
DATA = BASE / "data"
OUTDIR = DATA / "daily"
OUTDIR.mkdir(parents=True, exist_ok=True)

with rasterio.open(DATA / "gee_export" / REFERENCE) as ref:
    REF_TRANSFORM, REF_CRS = ref.transform, ref.crs
    REF_W, REF_H = ref.width, ref.height
masks = np.load(DATA / "masks_matched.npz")
H, W = int(masks["H"]), int(masks["W"])

print("=" * 78)
print("城乡掩膜 → EE asset")
print("=" * 78)
print(f"  参考网格 {REFERENCE}: {REF_W}×{REF_H} · {REF_CRS}")
print(f"  变换 {tuple(round(v, 9) for v in REF_TRANSFORM)[:6]}")
# 断言：参考文件的尺寸必须与掩膜一致，否则平铺下标会错位
assert (REF_W, REF_H) == (W, H), f"参考网格 {REF_W}×{REF_H} ≠ 掩膜 {W}×{H}"
print(f"  掩膜网格 {W}×{H} ✓ 与参考一致")

manifest = {}

for year in ARGS.years:
    cities = sorted(int(k.split("|")[1]) for k in masks.files if k.startswith(f"u{year}|"))
    if not cities:
        print(f"\n{year}: 无掩膜键，跳过")
        continue

    cube = np.zeros((len(cities), H, W), dtype=np.uint8)
    stats, overlap_total = [], 0
    for i, cid in enumerate(cities):
        u = masks[f"u{year}|{cid}"]
        r = masks[f"r{year}|{cid}"]
        # 城乡必须互斥。若不互斥，同一像元会被两种 zone 各算一次，
        # 而 reduce 两侧都把它算进去——SUHII 被稀释却看不出来。
        overlap = np.intersect1d(u, r, assume_unique=False)
        overlap_total += overlap.size
        flat = cube[i].reshape(-1)
        flat[u] = 1
        flat[r] = 2
        stats.append((cid, u.size, r.size))

    if overlap_total:
        raise RuntimeError(f"{year}: 城乡掩膜有 {overlap_total} 个重叠像元，"
                           f"与 12 的构造（rural ⊂ ~urban）矛盾，停止")

    # bbox 取**掩膜自身的**范围，不用 GAUL 边界：掩膜是权威，
    # 用多边形会引入一份与 asset 无关的几何，且 GAUL 版本可变。
    nz = np.flatnonzero(cube.any(axis=0))
    rows, cols = nz // W, nz % W
    x0, y0 = REF_TRANSFORM.c, REF_TRANSFORM.f
    px, py = REF_TRANSFORM.a, -REF_TRANSFORM.e
    bbox = [x0 + cols.min() * px, y0 - (rows.max() + 1) * py,
            x0 + (cols.max() + 1) * px, y0 - rows.min() * py]

    name = f"{ASSET_PREFIX}_{year}"
    tif = OUTDIR / f"{name}.tif"
    with rasterio.open(tif, "w", driver="GTiff", height=H, width=W,
                       count=len(cities), dtype="uint8", crs=REF_CRS,
                       transform=REF_TRANSFORM, compress="DEFLATE",
                       tiled=True, blockxsize=256, blockysize=256) as dst:
        dst.write(cube)
        for i, cid in enumerate(cities):
            dst.set_band_description(i + 1, f"c{cid}")

    u_tot = sum(s[1] for s in stats)
    r_tot = sum(s[2] for s in stats)
    print(f"\n{year}: {len(cities)} 城 · 城区 {u_tot:,} 像元 · 乡村 {r_tot:,} 像元")
    print(f"  乡村/城区 像元比 中位 {np.median([s[2] / s[1] for s in stats]):.2f}")
    print(f"  写出 {tif.relative_to(BASE)}  ({tif.stat().st_size / 1e6:.2f} MB)")
    print(f"  bbox {[round(v, 5) for v in bbox]}")

    # 回读校验：写出的 tif 必须与 npz 逐像元相同。写错方向（行列转置）
    # 不会报错，只会让城乡区域整体错位。
    with rasterio.open(tif) as chk:
        assert chk.transform == REF_TRANSFORM, "写出的 transform 与参考不符"
        back = chk.read()
    for i, cid in enumerate(cities):
        # ⚠ 不能用 np.array_equal(sel, 1)：它要求**形状相同**，标量 1 的形状是 ()，
        # 故对任何非 0 维数组都返回 False（实测 numpy 2.3.5，文档即如此）。
        # 曾因此误判「tif 写错」，实际数据完全正确。一律用显式逐元素比较。
        flat_back = back[i].reshape(-1)
        assert (flat_back[masks[f"u{year}|{cid}"]] == 1).all(), f"{year} city {cid} 城区回读不符"
        assert (flat_back[masks[f"r{year}|{cid}"]] == 2).all(), f"{year} city {cid} 乡村回读不符"
    assert int(back.sum()) == u_tot + 2 * r_tot, "回读校验和与写入不符"
    print(f"  回读校验 ✓ 逐像元与 masks_matched.npz 相同")

    manifest[year] = {"asset": f"projects/{PROJECT}/assets/{name}", "tif": tif.name,
                      "bands": len(cities), "cities": cities, "bbox": bbox,
                      "crsTransform": list(REF_TRANSFORM)[:6],
                      "urban_px": u_tot, "rural_px": r_tot}

(OUTDIR / "zones_manifest.json").write_text(
    json.dumps({"project": PROJECT, "reference": REFERENCE, "years": manifest},
               ensure_ascii=False, indent=2), encoding="utf-8")
print(f"\n清单已写 {(OUTDIR / 'zones_manifest.json').relative_to(BASE)}")

if ARGS.dry_run:
    print("\n[--dry-run] 未上传。")
    raise SystemExit(0)

if not EE_CLI.exists():
    print(f"\n✗ 找不到 earthengine CLI: {EE_CLI}", file=sys.stderr)
    raise SystemExit(2)

print(f"\n{'=' * 78}\n上传\n{'=' * 78}")
for year, info in manifest.items():
    cmd = [str(EE_CLI), "upload", "image",
           "--asset_id", info["asset"],
           "--pyramiding_policy", "mode",
           str(OUTDIR / info["tif"])]
    print(f"  {info['asset']}")
    proc = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    out = (proc.stdout or "").strip()
    err = (proc.stderr or "").strip()
    if proc.returncode != 0:
        # 权限类失败单独标出来：asset 写权限不足与凭据失效的处置不同
        hint = "（疑似无 asset 写权限 → 按约定停止并转 B 方案，不重试）" \
            if any(w in err.lower() for w in ("permission", "forbidden", "not authorized", "403")) else ""
        print(f"    ✗ 上传失败 rc={proc.returncode} {hint}\n    {err[:400]}", file=sys.stderr)
        raise SystemExit(3)
    print(f"    ✓ {out[:200] or '已提交'}")

print(f"\n上传任务已提交。asset 状态查询：")
print(f"  {EE_CLI} task list")
print("asset 就绪后再跑 02_verify_grid_alignment.py —— 网格对齐是硬前置，未过不得导出。")
