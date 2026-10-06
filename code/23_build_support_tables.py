"""Build three manuscript tables that expose the numbers behind existing claims.

Table 5  QA tier counts by year  — the per-year evidence behind the §3.4 asymmetry claim,
         read from the archived `qc_tier_sensitivity_{day,night}.csv` (tier `eq0` rows).
Table 6  City-level steps        — the per-city values behind Figure 4, from
         `same_source_decomposition_city_year.csv` (primary rule, city-equal means).
Table 7  Monthly support structure — when a city-month is estimable at all, from
         `daily_records_{day,night}.csv` and the Section 2.3 estimability rule.

Nothing is re-estimated: every cell is read from a frozen file. The script prints the same
rows again in Chinese for the Chinese manuscript, using a pinyin→中文 city map that is
asserted to cover all 41 units, so the two tables cannot drift apart in shape.
"""

from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parents[2]
TAB = BASE / "tables"
DATA = BASE / "data"
OUT = BASE / "outputs" / "daily"
ARCH = TAB / "_archive_pre_major_revision"
MIN_PX, MIN_DAYS, MIN_PERIODS = 20, 6, 2

ZH_CITY = {
    "Shanghai": "上海", "Nanjing": "南京", "Wuxi": "无锡", "Xuzhou": "徐州", "Changzhou": "常州",
    "Suzhou": "苏州", "Nantong": "南通", "Lianyungang": "连云港", "Huaian": "淮安",
    "Yancheng": "盐城", "Yangzhou": "扬州", "Zhenjiang": "镇江", "Taizhou": "泰州",
    "Suqian": "宿迁", "Hangzhou": "杭州", "Ningbo": "宁波", "Wenzhou": "温州",
    "Taizhou_Zhejiang": "台州",
    "Jiaxing": "嘉兴", "Huzhou": "湖州", "Shaoxing": "绍兴", "Jinhua": "金华",
    "Quzhou": "衢州", "Zhoushan": "舟山", "Lishui": "丽水", "Hefei": "合肥", "Wuhu": "芜湖",
    "Bengbu": "蚌埠", "Huainan": "淮南", "Maanshan": "马鞍山", "Huaibei": "淮北",
    "Tongling": "铜陵", "Anqing": "安庆", "Huangshan": "黄山", "Chuzhou": "滁州",
    "Fuyang": "阜阳", "Suzhou_Anhui": "宿州", "Luan": "六安", "Bozhou": "亳州",
    "Chizhou": "池州", "Xuancheng": "宣城",
}
ZH_PROV = {"Shanghai": "上海", "Jiangsu": "江苏", "Zhejiang": "浙江", "Anhui": "安徽"}
# 英文表统一用标准撇号写法；daily/covariates 里是 Huaian / Luan（无撇号）。
EN_NAME = {"Huaian": "Huai'an", "Luan": "Lu'an"}


def en_name(name: str) -> str:
    return EN_NAME.get(name, name.replace("_", " "))


def fmt(v, nd=2, dash="—"):
    return dash if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:+.{nd}f}"


def pct(v, dash="—"):
    return dash if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{100 * v:.1f}%"


def table5() -> list:
    rows = []
    for mode in ("day", "night"):
        df = pd.read_csv(ARCH / f"qc_tier_sensitivity_{mode}.csv")
        q = df[df.tier == "eq0"].set_index("year")
        for yr in sorted(q.index):
            r = q.loc[yr]
            rows.append((mode, int(yr), int(r.n_cities), r.urban_retention_vs_main,
                         r.weighting_effect_K))
    day = {r[1]: r for r in rows if r[0] == "day"}
    night = {r[1]: r for r in rows if r[0] == "night"}
    body = ["| Year | Day cities | Day urban-pixel retention | Day weighting effect (K) | "
            "Night cities | Night urban-pixel retention | Night weighting effect (K) |",
            "|---|---|---|---|---|---|---|"]
    for yr in sorted(day):
        d, n = day[yr], night[yr]
        body.append(f"| {yr} | {d[2]} | {pct(d[3])} | {fmt(d[4])} | {n[2]} | {pct(n[3])} | {fmt(n[4])} |")
    return body


def city_frame() -> tuple:
    cols = ["basis", "year", "mode", "city_id", "pixel_rule", "min_daily_days",
            "min_8day_periods", "delta_A_to_B_common_months_K", "delta_B_to_C_common_months_K"]
    df = pd.read_csv(OUT / "same_source_decomposition_city_year.csv", usecols=cols)
    q = df[(df.basis == "dynamic") & (df.pixel_rule == "20px") &
           (df.min_daily_days == MIN_DAYS) & (df.min_8day_periods == MIN_PERIODS)]
    q = q[q.delta_A_to_B_common_months_K.notna() & q.delta_B_to_C_common_months_K.notna()]
    m = q.groupby(["mode", "city_id"])[
        ["delta_A_to_B_common_months_K", "delta_B_to_C_common_months_K"]].mean()
    cov = pd.read_csv(DATA / "city_covariates.csv")
    cov = cov[cov.year == 2020][["city_id", "city"]]
    # 值表 m 的列是 MultiIndex（metric, mode）。**元数据单独放一个单层列的 DataFrame**：
    # 往 MultiIndex 列的帧里赋值会生成 ('name','') 这类复合列，随后 r['key'] 取回的是
    # 长度 1 的 Series 而不是标量（2026-10-04 实跑踩到 unhashable type: 'Series'）。
    frame = m.unstack("mode")
    names = cov.set_index("city_id")
    meta = pd.DataFrame({"city": [names.loc[c, "city"] for c in frame.index]}, index=frame.index)
    meta["prov"] = meta.city.str.split("_").str[0]
    meta["name"] = meta.city.str.split("_").str[-1].str.replace("'", "", regex=False)
    # 两处拼音同名必须按省份拆开，否则中文表会出现两行同一个中文名：
    #   江苏苏州 Suzhou / 安徽宿州 Suzhou；江苏泰州 Taizhou / 浙江台州 Taizhou（不同汉字）。
    meta["key"] = meta.name
    meta.loc[(meta.name == "Suzhou") & (meta.prov == "Anhui"), "key"] = "Suzhou_Anhui"
    meta.loc[(meta.name == "Taizhou") & (meta.prov == "Zhejiang"), "key"] = "Taizhou_Zhejiang"
    missing = sorted(set(meta.key) - set(ZH_CITY))
    assert not missing, f"中文城市名缺失：{missing}"
    return frame, meta


def table6(frame: pd.DataFrame, meta: pd.DataFrame) -> tuple:
    body = ["| City | Province | Day A−B (K) | Night A−B (K) | Day B−C (K) | Night B−C (K) |",
            "|---|---|---|---|---|---|"]
    zh = ["| 城市 | 省份 | 白天 A−B（K） | 夜间 A−B（K） | 白天 B−C（K） | 夜间 B−C（K） |",
          "|---|---|---|---|---|---|"]
    for cid, r in meta.sort_values(["prov", "key"]).iterrows():
        def val(metric, mode, f=frame):
            v = f.loc[cid, (metric, mode)]
            return "—" if pd.isna(v) else f"{v:+.2f}"
        body.append(f"| {en_name(r['name'])} | {r['prov']} | "
                    f"{val('delta_A_to_B_common_months_K', 'day')} | "
                    f"{val('delta_A_to_B_common_months_K', 'night')} | "
                    f"{val('delta_B_to_C_common_months_K', 'day')} | "
                    f"{val('delta_B_to_C_common_months_K', 'night')} |")
        zh.append(f"| {ZH_CITY[r['key']]} | {ZH_PROV[r['prov']]} | "
                  f"{val('delta_A_to_B_common_months_K', 'day')} | "
                  f"{val('delta_A_to_B_common_months_K', 'night')} | "
                  f"{val('delta_B_to_C_common_months_K', 'day')} | "
                  f"{val('delta_B_to_C_common_months_K', 'night')} |")
    return body, zh


def fixed_base(basis: str) -> tuple:
    """同一张城市级表，换成另一套掩膜口径——用来检验图 4 的符号格局是否随掩膜时相改变。"""
    cols = ["basis", "year", "mode", "city_id", "pixel_rule", "min_daily_days",
            "min_8day_periods", "delta_A_to_B_common_months_K", "delta_B_to_C_common_months_K"]
    df = pd.read_csv(OUT / "same_source_decomposition_city_year.csv", usecols=cols)
    q = df[(df.basis == basis) & (df.pixel_rule == "20px") &
           (df.min_daily_days == MIN_DAYS) & (df.min_8day_periods == MIN_PERIODS)]
    q = q[q.delta_A_to_B_common_months_K.notna() & q.delta_B_to_C_common_months_K.notna()]
    m = q.groupby(["mode", "city_id"])[
        ["delta_A_to_B_common_months_K", "delta_B_to_C_common_months_K"]].mean().unstack("mode")
    cov = pd.read_csv(DATA / "city_covariates.csv")
    cov = cov[cov.year == 2020][["city_id", "city"]]
    names = cov.set_index("city_id")
    meta = pd.DataFrame({"city": [names.loc[c, "city"] for c in m.index]}, index=m.index)
    meta["prov"] = meta.city.str.split("_").str[0]
    meta["name"] = meta.city.str.split("_").str[-1].str.replace("'", "", regex=False)
    meta["key"] = meta.name
    meta.loc[(meta.name == "Suzhou") & (meta.prov == "Anhui"), "key"] = "Suzhou_Anhui"
    meta.loc[(meta.name == "Taizhou") & (meta.prov == "Zhejiang"), "key"] = "Taizhou_Zhejiang"
    return m, meta


def table8(m: pd.DataFrame, meta: pd.DataFrame) -> tuple:
    body = ["| City | Province | Day A−B (K) | Night A−B (K) | Day B−C (K) | Night B−C (K) |",
            "|---|---|---|---|---|---|"]
    zh = ["| 城市 | 省份 | 白天 A−B（K） | 夜间 A−B（K） | 白天 B−C（K） | 夜间 B−C（K） |",
          "|---|---|---|---|---|---|"]
    for cid, r in meta.sort_values(["prov", "key"]).iterrows():
        def val(metric, mode):
            v = m.loc[cid, (metric, mode)]
            return "—" if pd.isna(v) else f"{v:+.2f}"
        cells = [val("delta_A_to_B_common_months_K", "day"), val("delta_A_to_B_common_months_K", "night"),
                 val("delta_B_to_C_common_months_K", "day"), val("delta_B_to_C_common_months_K", "night")]
        body.append(f"| {en_name(r['name'])} | {r['prov']} | " + " | ".join(cells) + " |")
        zh.append(f"| {ZH_CITY[r['key']]} | {ZH_PROV[r['prov']]} | " + " | ".join(cells) + " |")
    return body, zh


def main() -> None:
    b5 = table5()
    frame, meta = city_frame()
    b6, z6 = table6(frame, meta)
    # 表号与文件名的对应（2026-10-04 重构后）：正文表 6 = 逐年严格质量；表 8 = 城市级步骤；
    # 表 10 = 固定掩膜城市级；月度支持结构移到 24_support_diagnostics.py。
    specs = [("table6_00_qa_tier_by_year.md", "**Table 6.** Strict-quality support by year.",
              b5, None),
             ("table8_00_city_level_steps.md", "**Table 8.** City-level same-source steps.",
              b6, z6),
             ]
    for fname, cap, body, _ in specs:
        note = {
            "table6_00_qa_tier_by_year.md":
                "Note. Rows list only the years with a raw daily quality band. Cities are those "
                "retaining at least six shared months under the strict tier; retention is the "
                "strict-to-main ratio of urban pixel-days; the weighting effect is the strict "
                "minus main difference in the equal-weight level, reported per year.",
            "table8_00_city_level_steps.md":
                "Note. City-equal means over 2001–2025 under the primary rule; a dash marks a "
                "city with no estimable city-year for that overpass. The sign convention follows "
                "the text: day A−B is negative in all 37 cities with a daytime value.",
        }[fname]
        text = "\n".join([cap, ""] + body + ["", note, ""])
        (TAB / fname).write_text(text, encoding="utf-8")
        print(f"已写 tables/{fname} · 行 {len(body) - 2}")
    # 表 6 的单元格是纯数值，中英逐字相同；只有表头需要译，故这里只译表头再复用行。
    z5 = ["| 年份 | 白天城市数 | 白天城区像元保留率 | 白天等权效应（K） | 夜间城市数 | "
          "夜间城区像元保留率 | 夜间等权效应（K） |", "|---|---|---|---|---|---|---|"] + b5[2:]
    print("\n=== 中文表 6（逐年严格质量，正文表 6） ===")
    print("\n".join(z5))
    fm, fmeta = fixed_base("fixed2020")
    b8, z8 = table8(fm, fmeta)
    (TAB / "table10_00_city_level_fixed_mask.md").write_text(
        "\n".join(["**Table 10.** City-level same-source steps under the fixed 2020 mask."] + [""] + b8 +
                  ["", "Note. Same construction as Table 6, but with the urban and rural masks held "
                       "at their 2020 epoch for every year instead of following the five-year GHSL "
                       "epochs.", ""]), encoding="utf-8")
    print(f"\n已写 tables/table10_00_city_level_fixed_mask.md · 行 {len(b8) - 2}")
    print("\n=== 中文表 8（动态掩膜，正文表 8） ===")
    print("\n".join(z6))
    print("\n=== 中文表 10（固定 2020 掩膜，正文表 10） ===")
    print("\n".join(z8))
    agree = []
    for mode in ("day", "night"):
        for metric, sign in (("delta_A_to_B_common_months_K", -1),
                             ("delta_B_to_C_common_months_K", +1)):
            a = frame[(metric, mode)].dropna()
            b = fm[(metric, mode)].reindex(a.index).dropna()
            same = int((np.sign(a.reindex(b.index)) == np.sign(b)).sum())
            agree.append((mode, metric, len(b), same))
            print(f"  口径一致性 {mode:5s} {metric[:10]:10s} 可比城市 {len(b):2d} · 同号 {same:2d}")
    day = frame[("delta_A_to_B_common_months_K", "day")]
    print(f"\n核对：白天有值城市 {day.notna().sum()} · 全负 {bool((day.dropna() < 0).all())} · "
          f"最小 {day.min():.2f}")


if __name__ == "__main__":
    main()
