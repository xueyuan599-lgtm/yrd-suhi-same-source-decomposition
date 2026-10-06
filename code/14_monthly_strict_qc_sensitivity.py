"""Monthly equal-weighted paired MOD11A1 sensitivity to the strict QC tier.
Uses saved daily city summaries; years without raw Q bands (2020, 2023) are excluded.
"""
import argparse, re
from pathlib import Path
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parents[2]
DATA=BASE/"data"
OUT=BASE/"outputs"/"daily"
RULES={"any":(1,0.),"20px":(20,0.),"40pct":(10,.4),"60pct":(10,.6)}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--min-days",nargs="+",type=int,default=[3,6,10])
    args=ap.parse_args(); masks=np.load(DATA/"masks_matched.npz")
    months=[]; cityyears=[]
    for mode in ("night","day"):
        d=pd.read_csv(DATA/"daily"/f"daily_records_{mode}.csv",dtype={"date":str})
        d["date"]=pd.to_datetime(d.date,format="%Y%m%d",errors="coerce")
        d=d[d.date.notna()].copy(); d["month"]=d.date.dt.month
        d["qc_available"]=d.u_mean_strict.notna() & d.r_mean_strict.notna()
        d=d[d.qc_available]
        for (year,cid,month),g in d.groupby(["year","city_id","month"],sort=True):
            ku=f"u{int(year)}|{int(cid)}"; kr=f"r{int(year)}|{int(cid)}"
            if ku not in masks or kr not in masks: continue
            nu,nr=len(masks[ku]),len(masks[kr])
            for rule,(minpx,frac) in RULES.items():
                tu=max(minpx,frac*nu); tr=max(minpx,frac*nr)
                for min_days in args.min_days:
                    main_u=g.u_n.ge(tu); main_r=g.r_n.ge(tr); main_pair=main_u & main_r
                    strict_u=g.u_n_strict.ge(tu); strict_r=g.r_n_strict.ge(tr); strict_pair=strict_u & strict_r
                    nm=int(main_pair.sum()); ns=int(strict_pair.sum())
                    main=(float((g.loc[main_pair,"u_mean"]-g.loc[main_pair,"r_mean"]).mean()) if nm>=min_days else np.nan)
                    strict=(float((g.loc[strict_pair,"u_mean_strict"]-g.loc[strict_pair,"r_mean_strict"]).mean()) if ns>=min_days else np.nan)
                    months.append({"mode":mode,"year":int(year),"city_id":int(cid),"month":int(month),
                                   "pixel_rule":rule,"min_days":min_days,"main_paired_days":nm,
                                   "strict_paired_days":ns,"main_A_K":main,"strict_A_K":strict,
                                   "strict_minus_main_K":strict-main if np.isfinite(strict) and np.isfinite(main) else np.nan})
    m=pd.DataFrame(months)
    annual=[]
    for key,g in m.groupby(["mode","year","city_id","pixel_rule","min_days"],sort=True):
        mode,year,cid,rule,md=key
        common=g.dropna(subset=["main_A_K","strict_A_K"])
        n=int(len(common)); nmain=int(g.main_A_K.notna().sum()); nstrict=int(g.strict_A_K.notna().sum())
        annual.append({"mode":mode,"year":year,"city_id":cid,"pixel_rule":rule,"min_days":md,
                       "main_months":nmain,"strict_months":nstrict,"common_months":n,
                       "main_12month_K":float(g.main_A_K.mean()) if nmain==12 else np.nan,
                       "strict_12month_K":float(g.strict_A_K.mean()) if nstrict==12 else np.nan,
                       "strict_minus_main_12month_K":float(g.strict_A_K.mean()-g.main_A_K.mean()) if nmain==12 and nstrict==12 else np.nan,
                       "main_common_months6_K":float(common.main_A_K.mean()) if n>=6 else np.nan,
                       "strict_common_months6_K":float(common.strict_A_K.mean()) if n>=6 else np.nan,
                       "strict_minus_main_common_months6_K":float(common.strict_minus_main_K.mean()) if n>=6 else np.nan})
    a=pd.DataFrame(annual)
    OUT.mkdir(parents=True,exist_ok=True)
    m.to_csv(OUT/"strict_qc_monthly_paired.csv",index=False)
    a.to_csv(OUT/"strict_qc_city_year.csv",index=False)
    summ=(a[(a.pixel_rule=="20px")&(a.min_days==6)]
          .groupby(["mode","year"],as_index=False)
          .agg(city_years=("city_id","nunique"),complete12=("strict_12month_K",lambda x:int(x.notna().sum())),
               partial6=("strict_minus_main_common_months6_K",lambda x:int(x.notna().sum())),
               strict_mean_12=("strict_12month_K","mean"),main_mean_12=("main_12month_K","mean"),
               strict_minus_main_12=("strict_minus_main_12month_K","mean"),
               strict_minus_main_partial6=("strict_minus_main_common_months6_K","mean")))
    summ.to_csv(OUT/"strict_qc_year_summary.csv",index=False)
    print("monthly",len(m),"city-years",len(a),"QC-complete years",sorted(a.year.unique()),"outputs",OUT)
if __name__=="__main__": main()
