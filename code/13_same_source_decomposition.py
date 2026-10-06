"""Same-source daily versus 8-day decomposition for the Remote Sensing revision.

Primary rules fixed before reading new outcome tables:
- MOD11A1 QC bits 0-1 <= 1 (the daily exports are already masked this way).
- Main spatial gate: >=20 valid 1-km grid cells on each side on a valid date/period.
- A/B monthly gate: >=6 days (paired for A; each side independently for B).
- C/D monthly gate: >=2 paired 8-day periods; annual estimate requires all 12 months.
- Monthly estimates receive equal weight in the annual city estimate.
The existing fixed 2020 urban/rural masks are a prespecified mask sensitivity.
"""
from __future__ import annotations
import argparse, json, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from scipy.stats import theilslopes

BASE = Path(__file__).resolve().parents[2]  # analysis/
EXPORT = BASE / "data" / "gee_export"
DATA = BASE / "data"
OUT = BASE / "outputs" / "daily"
FILL = -9999
PIXEL_RULES = {"any": (1, 0.0), "20px": (20, 0.0), "40pct": (10, .40), "60pct": (10, .60)}


def get_masks(year, basis, masks):
    mask_year = year if basis == "dynamic" else 2020
    out = {}
    for k in masks.files:
        m = re.fullmatch(r"([ur])(\d{4})\|(\d+)", k)
        if not m or int(m.group(2)) != mask_year:
            continue
        side, _, cid = m.groups()
        out.setdefault(int(cid), {})[side] = masks[k].astype(np.int64, copy=False)
    return {cid: v for cid, v in out.items() if "u" in v and "r" in v}


def open_checked(path, ref):
    ds = rasterio.open(path)
    if (ds.width, ds.height, ds.crs, ds.transform) != (ref.width, ref.height, ref.crs, ref.transform):
        ds.close()
        raise ValueError(f"Grid mismatch: {path.name}")
    return ds


def zone_stats(arr, masks):
    flat = arr.reshape(-1)
    out = {}
    for cid, pair in masks.items():
        vals = {}
        for side in ("u", "r"):
            pix = flat[pair[side]]
            good = pix != FILL
            n = int(good.sum())
            vals[side + "_n"] = n
            vals[side + "_mean"] = float(pix[good].mean() * .02) if n else np.nan
        out[cid] = vals
    return out


def threshold(n, mask_n, rule):
    minimum, fraction = PIXEL_RULES[rule]
    return n >= max(minimum, fraction * mask_n)


def daily_zone_rows(mode, year, mask_sets):
    """Read each source raster once, then aggregate under the requested dynamic/fixed masks."""
    rows, c_sums, c_counts = [], {}, {}
    ref_path = EXPORT / f"uhi_{mode}daily_{year}_h1.tif"
    with rasterio.open(ref_path) as ref:
        comp_path = EXPORT / f"uhi_{mode}seq_{year}.tif"
        with open_checked(comp_path, ref) as comp:
            comp_dates = [pd.to_datetime(d[1:]) for d in comp.descriptions]
            jan1 = pd.Timestamp(year, 1, 1)
            if not comp_dates or comp_dates[0] != jan1 or any((d-jan1).days % 8 for d in comp_dates):
                raise ValueError(f"Unexpected 8-day date labels for {mode} {year}")
            if any((b-a).days != 8 for a,b in zip(comp_dates,comp_dates[1:])):
                print(f"WARNING missing native composite period(s): {mode} {year}",flush=True)
            c_sums = {d: np.zeros((ref.height, ref.width), np.float32) for d in comp_dates}
            c_counts = {d: np.zeros((ref.height, ref.width), np.uint8) for d in comp_dates}
            native = comp.read()  # Native MOD11A2 raw DN; zone_stats applies the 0.02 K scale.
            for half in ("h1", "h2"):
                p = EXPORT / f"uhi_{mode}daily_{year}_{half}.tif"
                with open_checked(p, ref) as daily:
                    data = daily.read()
                    for i, desc in enumerate(daily.descriptions):
                        if not desc or not desc.startswith("D"):
                            continue  # Q* bands are already represented by the server-side D* mask.
                        date = pd.to_datetime(desc[1:])
                        ci = (date-jan1).days // 8
                        start = jan1 + pd.Timedelta(days=8*ci)
                        a = data[i]
                        valid = a != FILL
                        if start in c_sums:  # Keep missing native MOD11A2 windows missing.
                            np.add(c_sums[start], a, out=c_sums[start], where=valid, casting="unsafe")
                            c_counts[start][valid] += 1
                        for basis, masks in mask_sets.items():
                            for cid, v in zone_stats(a, masks).items():
                                rows.append({"basis":basis,"year":year,"mode":mode,"city_id":cid,
                                             "date":date,**v})
            comp_rows = []
            for j, start in enumerate(comp_dates):
                cnt = c_counts[start]
                c_arr = np.full(cnt.shape, FILL, np.int16)
                good = cnt > 0
                c_arr[good] = np.rint(c_sums[start][good] / cnt[good]).astype(np.int16)
                for basis, masks in mask_sets.items():
                    c_stats, d_stats = zone_stats(c_arr, masks), zone_stats(native[j], masks)
                    for cid in set(c_stats) & set(d_stats):
                        comp_rows.append({"basis":basis,"year":year,"mode":mode,"city_id":cid,
                                          "date":start,**{"c_"+k:v for k,v in c_stats[cid].items()},
                                          **{"d_"+k:v for k,v in d_stats[cid].items()}})
    return pd.DataFrame(rows), pd.DataFrame(comp_rows)

def summarize_monthly(daily, comp, mask_sizes, basis, min_days, min_comps):
    rows = []
    for (year, cid, month), day in daily.assign(month=daily.date.dt.month).groupby(["year", "city_id", "month"]):
        sizes = mask_sizes.get((year, cid), (np.nan, np.nan))
        for rule in PIXEL_RULES:
            ud = day.apply(lambda r: threshold(r.u_n, sizes[0], rule), axis=1)
            rd = day.apply(lambda r: threshold(r.r_n, sizes[1], rule), axis=1)
            both = ud & rd
            a_n = int(both.sum()); u_n = int(ud.sum()); r_n = int(rd.sum())
            a = float((day.loc[both, "u_mean"] - day.loc[both, "r_mean"]).mean()) if a_n >= min_days else np.nan
            b = (float(day.loc[ud, "u_mean"].mean() - day.loc[rd, "r_mean"].mean())
                 if u_n >= min_days and r_n >= min_days else np.nan)
            cr = comp[(comp.year == year) & (comp.city_id == cid) & (comp.date.dt.month == month)]
            cu = cr.apply(lambda r: threshold(r.c_u_n, sizes[0], rule), axis=1)
            crr = cr.apply(lambda r: threshold(r.c_r_n, sizes[1], rule), axis=1)
            du = cr.apply(lambda r: threshold(r.d_u_n, sizes[0], rule), axis=1)
            drr = cr.apply(lambda r: threshold(r.d_r_n, sizes[1], rule), axis=1)
            cpair, dpair = cu & crr, du & drr
            c_n, d_n = int(cpair.sum()), int(dpair.sum())
            c = float((cr.loc[cpair, "c_u_mean"] - cr.loc[cpair, "c_r_mean"]).mean()) if c_n >= min_comps else np.nan
            d = float((cr.loc[dpair, "d_u_mean"] - cr.loc[dpair, "d_r_mean"]).mean()) if d_n >= min_comps else np.nan
            cd_pair = cpair & dpair
            cd_n = int(cd_pair.sum())
            cd = (float(((cr.loc[cd_pair, "c_u_mean"]-cr.loc[cd_pair, "c_r_mean"])-
                          (cr.loc[cd_pair, "d_u_mean"]-cr.loc[cd_pair, "d_r_mean"])).mean())
                  if cd_n >= min_comps else np.nan)
            rows.append({"basis": basis, "year": int(year), "city_id": int(cid), "month": int(month),
                         "pixel_rule": rule, "min_daily_days": min_days, "min_8day_periods": min_comps,
                         "A_paired_daily_K": a, "A_n_days": a_n, "B_independent_dates_K": b,
                         "B_n_urban_days": u_n, "B_n_rural_days": r_n, "C_daily_rebuilt_8day_K": c,
                         "C_n_periods": c_n, "D_native_MOD11A2_K": d, "D_n_periods": d_n,
                         "C_minus_D_paired_K": cd, "CD_n_paired_periods": cd_n,
                         "delta_A_to_B_month_K": a-b if np.isfinite(a) and np.isfinite(b) else np.nan,
                         "delta_B_to_C_month_K": b-c if np.isfinite(b) and np.isfinite(c) else np.nan})
    return pd.DataFrame(rows)


def annualize(monthly):
    idc = ["basis", "year", "city_id", "pixel_rule", "min_daily_days", "min_8day_periods"]
    cols = ["A_paired_daily_K", "B_independent_dates_K", "C_daily_rebuilt_8day_K", "D_native_MOD11A2_K", "C_minus_D_paired_K"]
    rows=[]
    for key,g in monthly.groupby(idc,sort=True):
        rec=dict(zip(idc,key))
        rec["months_available"] = int(g[cols[:4]].notna().all(axis=1).sum())
        for col in cols:
            n=int(g[col].notna().sum()); rec["n_months_"+col]=n
            rec[col.replace("_K","_annual_K")] = float(g[col].mean()) if n==12 else np.nan
        rec["delta_A_to_B_K"] = rec.get("A_paired_daily_annual_K",np.nan)-rec.get("B_independent_dates_annual_K",np.nan)
        rec["delta_B_to_C_K"] = rec.get("B_independent_dates_annual_K",np.nan)-rec.get("C_daily_rebuilt_8day_annual_K",np.nan)
        rec["delta_C_to_D_K"] = rec.get("C_minus_D_paired_annual_K",np.nan)
        common=g.dropna(subset=cols[:4])
        rec["common_months_all_methods"] = int(len(common))
        for col in cols[:4]:
            rec[col.replace("_K","_common_months_mean_K")] = float(common[col].mean()) if len(common)>=6 else np.nan
        rec["delta_A_to_B_common_months_K"] = float(common[cols[0]].sub(common[cols[1]]).mean()) if len(common)>=6 else np.nan
        rec["delta_B_to_C_common_months_K"] = float(common[cols[1]].sub(common[cols[2]]).mean()) if len(common)>=6 else np.nan
        rec["delta_C_to_D_common_months_K"] = float(common[cols[4]].mean()) if len(common)>=6 else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)

def trend_rows(annual):
    metrics=["A_paired_daily_annual_K","B_independent_dates_annual_K","C_daily_rebuilt_8day_annual_K",
             "D_native_MOD11A2_annual_K","delta_A_to_B_K","delta_B_to_C_K","delta_C_to_D_K"]
    out=[]
    group_cols=["basis","mode","pixel_rule","min_daily_days","min_8day_periods"]
    for keys,g in annual.groupby(group_cols):
        basis,mode,rule,md,mc=keys
        complete=g[g[metrics].notna().all(axis=1)]
        year_counts=complete.groupby("year").city_id.nunique()
        years=sorted(year_counts[year_counts>=5].index.tolist())
        cities=sorted(c for c in complete.city_id.unique()
                      if set(years).issubset(set(complete.loc[complete.city_id==c,"year"]))) if years else []
        if len(cities)<5 or len(years)<8:
            out.append({"basis":basis,"mode":mode,"pixel_rule":rule,"min_daily_days":md,
                        "min_8day_periods":mc,"metric":"all","balanced_cities":len(cities),
                        "balanced_years":len(years),"Sen_K_per_decade":np.nan,"CI95_lo":np.nan,
                        "CI95_hi":np.nan,"HR_MK_p":np.nan})
            continue
        for metric in metrics:
            x=(complete[complete.city_id.isin(cities)&complete.year.isin(years)]
               .groupby("year")[metric].mean().sort_index())
            if len(x)<8: continue
            s=theilslopes(x.values,x.index.values,alpha=.95)
            pval=np.nan
            try:
                sys.path.insert(0,str(BASE/".deps")); import pymannkendall as mk
                pval=float(mk.hamed_rao_modification_test(x.values).p)
            except Exception: pass
            out.append({"basis":basis,"mode":mode,"pixel_rule":rule,"min_daily_days":md,
                        "min_8day_periods":mc,"metric":metric,"balanced_cities":len(cities),
                        "balanced_years":len(x),"Sen_K_per_decade":float(s.slope*10),
                        "CI95_lo":float(s.low_slope*10),"CI95_hi":float(s.high_slope*10),"HR_MK_p":pval})
    return pd.DataFrame(out)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--years",nargs="+",type=int,default=[2020,2023])
    ap.add_argument("--modes",nargs="+",choices=["night","day"],default=["night","day"])
    ap.add_argument("--bases",nargs="+",choices=["dynamic","fixed2020"],default=["dynamic","fixed2020"])
    ap.add_argument("--min-days",nargs="+",type=int,default=[3,6,10])
    ap.add_argument("--min-comps",nargs="+",type=int,default=[1,2,3])
    args=ap.parse_args(); OUT.mkdir(parents=True,exist_ok=True)
    masks=np.load(DATA/"masks_matched.npz")
    p_month=OUT/"same_source_decomposition_monthly.csv"
    p_annual=OUT/"same_source_decomposition_city_year.csv"
    p_comp=OUT/"same_source_decomposition_8day_city_period.csv"
    first_month=first_annual=first_comp=True
    for mode in args.modes:
        for year in args.years:
            mask_sets={b:get_masks(year,b,masks) for b in args.bases}
            mask_sets={b:v for b,v in mask_sets.items() if v}
            if not mask_sets: continue
            d,c=daily_zone_rows(mode,year,mask_sets)
            c.to_csv(p_comp,index=False,mode="w" if first_comp else "a",header=first_comp)
            first_comp=False
            month_chunks=[]; annual_chunks=[]
            for basis,m in mask_sets.items():
                db=d[d.basis.eq(basis)]; cb=c[c.basis.eq(basis)]
                sizes={(year,cid):(len(v["u"]),len(v["r"])) for cid,v in m.items()}
                for md in args.min_days:
                    for mc in args.min_comps:
                        mon=summarize_monthly(db,cb,sizes,basis,md,mc); mon["mode"]=mode
                        ann=annualize(mon); ann["mode"]=mode
                        month_chunks.append(mon); annual_chunks.append(ann)
            if month_chunks:
                pd.concat(month_chunks,ignore_index=True).to_csv(p_month,index=False,
                    mode="w" if first_month else "a",header=first_month); first_month=False
            if annual_chunks:
                pd.concat(annual_chunks,ignore_index=True).to_csv(p_annual,index=False,
                    mode="w" if first_annual else "a",header=first_annual); first_annual=False
            print(f"completed {mode} {year}",flush=True)
    if not first_annual:
        annual=pd.read_csv(p_annual)
        primary=annual[annual.basis.eq("dynamic")&annual.pixel_rule.eq("20px")&
                       annual.min_daily_days.eq(6)&annual.min_8day_periods.eq(2)]
        trend_rows(primary).to_csv(OUT/"same_source_decomposition_trends.csv",index=False)
    print(f"outputs: {p_month}, {p_annual}, {p_comp}",flush=True)
if __name__=="__main__": main()









