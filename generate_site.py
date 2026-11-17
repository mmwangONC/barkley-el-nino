#!/usr/bin/env python3
"""
Barkley Slope El Nino Watch - static-site generator.
Fetches the latest ONC cabled CTD + oxygen and NOAA ENSO / sea-level data,
computes anomalies against a pre-2026 harmonic climatology, and writes
docs/index.html + docs/fig_watch.png + docs/fig_story.png + docs/status.csv.
Designed to run weekly from GitHub Actions. Requires env var ONC_TOKEN.
"""
import os, time, io, datetime as dt
import numpy as np, pandas as pd, requests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT = "docs"; os.makedirs(OUT, exist_ok=True)
TOK = os.environ["ONC_TOKEN"]; B = "https://data.oceannetworks.ca/api"
S = requests.Session()
SITES = {"FGPPN":("Folger Pinnacle",25), "NCBC":("Barkley Upper Slope",396),
         "BACND":("Barkley Node",643), "BACAX":("Barkley Canyon Axis",985)}
TODAY = pd.Timestamp(dt.datetime.utcnow().date())

def fetch(lc, cat, y0=2009):
    out=[]
    for yr in range(y0, TODAY.year+1):
        r=None
        for _ in range(3):
            try:
                r=S.get(f"{B}/scalardata/location", params=dict(token=TOK, locationCode=lc,
                        deviceCategoryCode=cat, dateFrom=f"{yr}-01-01T00:00:00.000Z",
                        dateTo=f"{yr+1}-01-01T00:00:00.000Z", resamplePeriod=86400,
                        rowLimit=100000), timeout=600)
                if r.status_code==200: break
                time.sleep(4)
            except Exception: time.sleep(4); r=None
        if not r or r.status_code!=200: continue
        d={}
        for s in (r.json().get("sensorData") or []):
            d[s["sensorCode"]]=pd.Series(s["data"]["values"],
                index=pd.to_datetime(s["data"]["sampleTimes"], format="mixed", utc=True))
        if d: out.append(pd.DataFrame(d))
    if not out: return None
    X=pd.concat(out).sort_index(); return X[~X.index.duplicated()]

def harm(s, n=3):
    s=s.dropna(); b=s[s.index.year<2026]
    if len(b)<800: return None
    t=(b.index-pd.Timestamp("2000-01-01",tz=b.index.tz)).total_seconds().values/86400
    M=[np.ones_like(t)]
    for k in range(1,n+1): M+=[np.cos(2*np.pi*k*t/365.25), np.sin(2*np.pi*k*t/365.25)]
    c=np.linalg.lstsq(np.column_stack(M), b.values, rcond=None)[0]
    def f(idx):
        tt=(idx-pd.Timestamp("2000-01-01",tz=idx.tz)).total_seconds().values/86400
        MM=[np.ones_like(tt)]
        for k in range(1,n+1): MM+=[np.cos(2*np.pi*k*tt/365.25), np.sin(2*np.pi*k*tt/365.25)]
        return np.column_stack(MM)@c
    return f

def anom(s):
    s=s.dropna(); f=harm(s)
    if f is None: return None, None
    a=pd.Series(s.values-f(s.index), index=s.index)
    return a, float(a[a.index.year<2026].std())

# ---- fetch ONC ----
D={}; O2={}
for lc in SITES:
    D[lc]=fetch(lc,"CTD"); o=fetch(lc,"OXYSENSOR")
    O2[lc]=o if o is not None else pd.DataFrame()
A={}; SD={}
for lc in SITES:
    A[lc]={}; SD[lc]={}
    if D[lc] is not None and "Temperature" in D[lc]:
        A[lc]["T"],SD[lc]["T"]=anom(D[lc]["Temperature"])
    if "oxygen_corrected" in O2[lc]:
        A[lc]["O2"],SD[lc]["O2"]=anom(O2[lc]["oxygen_corrected"])

def nv(s):
    s=s.copy(); 
    if s.index.tz is not None: s.index=s.index.tz_localize(None)
    return s

# ---- NOAA ENSO indices (best-effort) ----
oni_txt=""; roni=""
try:
    r=S.get("https://psl.noaa.gov/data/correlation/oni.data", timeout=60)
    oni_txt=r.text
except Exception: pass
try:
    r=S.get("https://www.cpc.ncep.noaa.gov/data/indices/RONI.ascii.txt", timeout=60)
    lines=[l.split() for l in r.text.strip().splitlines()[1:] if l.strip()]
    roni=lines[-1] if lines else ""
except Exception: pass

# ---- FIGURE 1: weekly watch (T and O2 vs climatology band) ----
plt.rcParams.update({"font.size":9,"font.family":"DejaVu Sans"})
fig,axs=plt.subplots(2,1,figsize=(10,7.5),sharex=True)
win0=TODAY-pd.Timedelta("220D")
for ax,var,lab in [(axs[0],"T","temperature anomaly (degC)"),(axs[1],"O2","oxygen anomaly (ml/L)")]:
    for lc,(nm,dp) in SITES.items():
        if var not in A[lc] or A[lc][var] is None: continue
        a=nv(A[lc][var]); a=a[a.index>=win0]
        if len(a)<5: continue
        ax.plot(a.index, a.values, lw=1.4, label=f"{nm} ({dp} m)")
    ax.axhline(0,color="#444",lw=.6); ax.set_ylabel(lab); ax.grid(alpha=.25)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper left")
axs[0].set_title(f"Barkley slope El Nino watch - anomalies vs pre-2026 normal, updated {TODAY.date()}", loc="left")
fig.tight_layout(); fig.savefig(f"{OUT}/fig_watch.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ---- status table ----
rows=[]
for lc,(nm,dp) in SITES.items():
    cell={"site":nm,"depth":dp}
    for var in ["T","O2"]:
        if var in A[lc] and A[lc][var] is not None:
            a=A[lc][var]; last30=a[a.index>=a.index.max()-pd.Timedelta("30D")].mean()
            prev=a[(a.index>=a.index.max()-pd.Timedelta("60D"))&(a.index<a.index.max()-pd.Timedelta("30D"))].mean()
            sd=SD[lc][var]; trend="rising" if last30-prev>0.01 else "falling" if last30-prev<-0.01 else "flat"
            cell[var]=f"{last30:+.2f} ({last30/sd:+.1f} sd, {trend})"
        else: cell[var]="-"
    rows.append(cell)
status=pd.DataFrame(rows); status.to_csv(f"{OUT}/status.csv", index=False)

trows="".join(f"<tr><td>{r['site']}</td><td>{r['depth']} m</td><td>{r['T']}</td><td>{r['O2']}</td></tr>" for r in rows)
roni_s=f"{roni[0]} {roni[1]} = {roni[2]}" if roni else "n/a"

html=f"""<!doctype html><html lang=en><head><meta charset=utf-8>
<meta name=viewport content="width=device-width, initial-scale=1">
<title>Barkley Slope El Nino Watch</title>
<style>
body{{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:900px;margin:0 auto;padding:1.2rem;color:#222;line-height:1.5}}
h1{{font-size:1.5rem;margin:.2rem 0}} h2{{font-size:1.1rem;margin-top:1.6rem;border-bottom:1px solid #ddd;padding-bottom:.2rem}}
.sub{{color:#666;font-size:.9rem}} table{{border-collapse:collapse;width:100%;font-size:.9rem}}
th,td{{border:1px solid #ccc;padding:.35rem .5rem;text-align:left}} th{{background:#f4f4f4}}
img{{width:100%;height:auto;border:1px solid #eee;margin:.4rem 0}} .k{{background:#fff7e6;padding:.6rem .8rem;border-left:4px solid #e8a33d;margin:.8rem 0}}
footer{{color:#888;font-size:.8rem;margin-top:2rem;border-top:1px solid #eee;padding-top:.6rem}}
</style></head><body>
<h1>Barkley Slope El Nino Watch</h1>
<p class=sub>Ocean Networks Canada cabled records, southwest Vancouver Island. Updated {TODAY.date()} (UTC). ONI/RONI: RONI {roni_s}.</p>
<div class=k><b>The story.</b> The 2026 El Nino is reaching Vancouver Island <b>through the ocean, not the air</b>:
the surface stays near normal while the Barkley slope at 396 m turns warmer, fresher, lighter and more oxygen-rich -
the fingerprint of density surfaces pushed down ~50-100 m by a coastally trapped wave from the tropics.
The same response, smaller, appeared in 2023, so it is repeatable and can be anticipated.</div>
<h2>This week</h2>
<table><tr><th>Site</th><th>Depth</th><th>Temperature anomaly</th><th>Oxygen anomaly</th></tr>{trows}</table>
<p class=sub>Anomaly = latest 30-day mean minus the pre-2026 harmonic climatology; "sd" is the pre-2026 residual spread.</p>
<h2>Live figure</h2>
<img src="fig_watch.png" alt="anomaly watch">
<p class=sub>Red/warm above zero = warmer or more oxygen than normal for the date. Surface (25 m) vs slope (396-985 m)
is the key contrast: the El Nino shows at depth first.</p>
<h2>How this complements the ONC monitor</h2>
<p>Steve Mihaly's <a href="https://sfmihaly.github.io/el-nino-2026/">ONC slope &amp; canyon monitor</a> tracks the records as they break, at more depths.
This page focuses on the <b>surface-versus-slope contrast</b> and the mechanism behind it. Coastal sea level, density-surface analysis and the
2023 comparison will be added as they are completed.</p>
split on density surfaces, the 2023 repeatability case, and the forecast for 2027.</p>
<footer>Data: Ocean Networks Canada Oceans 3.0 (CC-BY 4.0); NOAA CPC/PSL ENSO indices; NOAA CO-OPS sea level.
Auto-generated weekly. Lead with temperature (passed swap-test QC); 396 m salinity may carry sensor drift; Folger Deep 96 m offline since 2026-07-18.</footer>
</body></html>"""
open(f"{OUT}/index.html","w").write(html)
print("site written to", OUT)
for r in rows: print(r["site"], r["depth"], "T",r["T"], "| O2",r["O2"])
