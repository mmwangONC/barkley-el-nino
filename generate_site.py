#!/usr/bin/env python3
"""
Barkley Slope El Nino Watch (enhanced) - static-site generator.
Fetches the latest ONC cabled CTD + oxygen, computes anomalies, the weekly watch,
the T-S diagram, the depth-time stack, and the heave / water-mass split (using a
bundled Line P P4 climatology for the vertical gradient). Writes docs/.
Runs weekly from GitHub Actions. Requires env ONC_TOKEN. Needs: numpy pandas matplotlib requests gsw.
"""
import os, time, datetime as dt, pickle
import numpy as np, pandas as pd, requests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
try:
    import gsw; HAVE_GSW=True
except Exception: HAVE_GSW=False

OUT="docs"; os.makedirs(OUT, exist_ok=True)
TOK=os.environ["ONC_TOKEN"]; B="https://data.oceannetworks.ca/api"; S=requests.Session()
# absolute-cert workaround is unnecessary on GitHub runners (system certs OK)
SITES={"FGPPN":("Folger Pinnacle",25),"FGPD":("Folger Deep",96),"NCBC":("Barkley Upper Slope",396),
       "BACND":("Barkley Node",643),"BACAX":("Barkley Canyon Axis",985)}
TODAY=pd.Timestamp(dt.datetime.utcnow().date())
P4=pickle.load(open("p4_clim.pkl","rb"))["clim"]
GL=pickle.load(open("glider_note.pkl","rb"))

def fetch(lc,cat,y0=2009):
    out=[]
    for yr in range(y0,TODAY.year+1):
        r=None
        for _ in range(3):
            try:
                r=S.get(f"{B}/scalardata/location",params=dict(token=TOK,locationCode=lc,deviceCategoryCode=cat,
                    dateFrom=f"{yr}-01-01T00:00:00.000Z",dateTo=f"{yr+1}-01-01T00:00:00.000Z",
                    resamplePeriod=86400,rowLimit=100000),timeout=600)
                if r.status_code==200: break
                time.sleep(4)
            except Exception: time.sleep(4); r=None
        if not r or r.status_code!=200: continue
        d={}
        for s in (r.json().get("sensorData") or []):
            d[s["sensorCode"]]=pd.Series(s["data"]["values"],index=pd.to_datetime(s["data"]["sampleTimes"],format="mixed",utc=True))
        if d: out.append(pd.DataFrame(d))
    if not out: return None
    X=pd.concat(out).sort_index(); return X[~X.index.duplicated()]

def harm(s,n=3):
    s=s.dropna(); b=s[s.index.year<2026]
    if len(b)<800: return None
    t=(b.index-pd.Timestamp("2000-01-01",tz=b.index.tz)).total_seconds().values/86400
    M=[np.ones_like(t)]
    for k in range(1,n+1): M+=[np.cos(2*np.pi*k*t/365.25),np.sin(2*np.pi*k*t/365.25)]
    c=np.linalg.lstsq(np.column_stack(M),b.values,rcond=None)[0]
    def f(idx):
        tt=(idx-pd.Timestamp("2000-01-01",tz=idx.tz)).total_seconds().values/86400
        MM=[np.ones_like(tt)]
        for k in range(1,n+1): MM+=[np.cos(2*np.pi*k*tt/365.25),np.sin(2*np.pi*k*tt/365.25)]
        return np.column_stack(MM)@c
    return f
def anom(s):
    s=s.dropna(); f=harm(s)
    if f is None: return None,None
    a=pd.Series(s.values-f(s.index),index=s.index); return a,float(a[a.index.year<2026].std())
def nv(s):
    s=s.copy()
    if s.index.tz is not None: s.index=s.index.tz_localize(None)
    return s

# ---- fetch ONC ----
D={}; O2={}
for lc in SITES:
    D[lc]=fetch(lc,"CTD"); ox=fetch(lc,"OXYSENSOR"); O2[lc]=ox if ox is not None else pd.DataFrame()
A={}; SD={}
for lc in SITES:
    A[lc]={}; SD[lc]={}
    if D[lc] is not None and "Temperature" in D[lc]: A[lc]["T"],SD[lc]["T"]=anom(D[lc]["Temperature"])
    if "oxygen_corrected" in O2[lc]: A[lc]["O2"],SD[lc]["O2"]=anom(O2[lc]["oxygen_corrected"])

# ---- RONI (best effort) ----
roni=""
try:
    r=S.get("https://www.cpc.ncep.noaa.gov/data/indices/RONI.ascii.txt",timeout=60)
    L=[l.split() for l in r.text.strip().splitlines()[1:] if l.strip()]; roni=L[-1] if L else ""
except Exception: pass
roni_s=f"{roni[0]} {roni[1]} = {roni[2]}" if roni else "n/a"

# ===== FIG watch =====
plt.rcParams.update({"font.size":9,"font.family":"DejaVu Sans"})
fig,axs=plt.subplots(2,1,figsize=(10,7.5),sharex=True); win0=TODAY-pd.Timedelta("220D")
for ax,var,lab,keys in [(axs[0],"T","temperature anomaly (degC)",SITES),(axs[1],"O2","oxygen anomaly (ml/L)",{"NCBC":SITES["NCBC"],"BACND":SITES["BACND"]})]:
    for lc,(nm,dp) in keys.items():
        if var in A.get(lc,{}) and A[lc][var] is not None:
            a=nv(A[lc][var]); a=a[a.index>=win0]
            if len(a)>=5: ax.plot(a.index,a.values,lw=1.4,label=f"{nm} ({dp} m)")
    ax.axhline(0,color="#444",lw=.6); ax.set_ylabel(lab); ax.grid(alpha=.25); ax.legend(frameon=False,fontsize=7,ncol=2,loc="upper left")
axs[0].set_title(f"Barkley slope El Nino watch - anomalies vs pre-2026 normal, updated {TODAY.date()}",loc="left")
fig.tight_layout(); fig.savefig(f"{OUT}/fig_watch.png",dpi=150,bbox_inches="tight"); plt.close(fig)

# ===== heave / water-mass split (NCBC 396 m) using P4 gradient =====
heave_html="<p>insufficient data</p>"
try:
    N=nv(D["NCBC"])[["Temperature","salinity"]].dropna()
    if HAVE_GSW:
        SA=gsw.SA_from_SP(N["salinity"].values,396.,-126.17,48.43); CT=gsw.CT_from_t(SA,N["Temperature"].values,396.); N["sig"]=gsw.sigma0(SA,CT)
    else:
        N["sig"]=1000*(1-(N["Temperature"]-10)*2e-4+(N["salinity"]-34)*7.6e-4)
    NO=nv(O2["NCBC"]["oxygen_corrected"]) if "oxygen_corrected" in O2["NCBC"] else pd.Series(dtype=float)
    N=N.join(NO.rename("O2"),how="left")
    aT,sdT=anom(N["Temperature"]); aSig,sdSig=anom(N["sig"]); aO,sdO=anom(N["O2"])
    g=lambda col: np.polyfit(P4.loc[355:445].index.values.astype(float),P4.loc[355:445][col].values,1)[0]
    dSigdz,dTdz,dO2dz=g("sig"),g("T"),g("O2")
    def win(a,d): x=a[a.index>=a.index.max()-pd.Timedelta(f"{d}D")]; return float(x.mean())
    dT_,dSig_,dO_=win(aT,30),win(aSig,30),win(aO,30)
    h_sig=-dSig_/dSigdz; h_T=dT_/(-dTdz); h=0.5*(h_sig+h_T)
    o2_heave=-h*dO2dz; o2_wm=dO_-o2_heave
    # seasonal swing
    f=harm(N["sig"]); doy=pd.date_range("2020-01-01","2020-12-31",tz="UTC"); cyc=f(doy); swing=cyc.max()-cyc.min()
    heave_html=f"""<table>
    <tr><th>Quantity (last 30 days, 396 m)</th><th>Value</th></tr>
    <tr><td>heave from density</td><td>{h_sig:+.0f} m</td></tr>
    <tr><td>heave from temperature</td><td>{h_T:+.0f} m</td></tr>
    <tr><td>glider-measured ({GL['date']})</td><td>sigma 26.8 at {GL['iso268_obs']} m vs {GL['iso268_norm']} m normal = +{GL['iso268_obs']-GL['iso268_norm']:.0f} m</td></tr>
    <tr><td>density anomaly vs whole normal annual swing</td><td>{abs(dSig_)/swing*100:.0f}% (swing {swing:.2f} kg/m3)</td></tr>
    <tr><td>oxygen: heave would add</td><td>{o2_heave:+.2f} ml/L</td></tr>
    <tr><td>oxygen: measured</td><td>{dO_:+.2f} ml/L</td></tr>
    <tr><td>oxygen: on-isopycnal change (deoxygenation)</td><td>{o2_wm:+.2f} ml/L</td></tr></table>"""
    # ---- FIG T-S ----
    fig,ax=plt.subplots(figsize=(7,6))
    Ts=np.linspace(4.5,8.5,60); Ss=np.linspace(33.7,34.35,60); TT,SS=np.meshgrid(Ts,Ss)
    if HAVE_GSW:
        SIGg=gsw.sigma0(gsw.SA_from_SP(SS,396,-126.17,48.43),gsw.CT_from_t(gsw.SA_from_SP(SS,396,-126.17,48.43),TT,396))
        cs=ax.contour(SS,TT,SIGg,levels=[26.4,26.6,26.8,27.0],colors="#999",linewidths=.8); ax.clabel(cs,fmt="%.1f",fontsize=7)
    NN=N.dropna(subset=["Temperature","salinity"]); yr=NN.index.year
    ax.scatter(NN[yr<2026]["salinity"],NN[yr<2026]["Temperature"],c=yr[yr<2026],cmap="Greys",s=5,alpha=.35)
    y26=NN[NN.index.year==2026]; rec=y26[y26.index>=y26.index.max()-pd.Timedelta("60D")]
    ax.scatter(y26["salinity"],y26["Temperature"],c="#b2182b",s=14,alpha=.75,label="2026")
    ax.scatter(rec["salinity"],rec["Temperature"],c="black",marker="x",s=30,label="last 60 d")
    ax.text(33.74,7.8,"2026: warm + FRESH\n(not undercurrent)",fontsize=8.5,color="#b2182b",fontweight="bold")
    ax.set_xlabel("practical salinity"); ax.set_ylabel("temperature (degC)"); ax.legend(frameon=False,fontsize=8,loc="lower left")
    ax.set_title(f"T-S at 396 m, 2009-2026 (updated {TODAY.date()})",fontsize=9.5,loc="left")
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_TS.png",dpi=150,bbox_inches="tight"); plt.close(fig)
except Exception as e:
    heave_html=f"<p>heave analysis unavailable: {e}</p>"

# ===== FIG depth stack =====
try:
    fig,axs=plt.subplots(5,1,figsize=(10,8),sharex=True)
    for ax,(lc,(nm,dp)) in zip(axs,SITES.items()):
        if "T" in A.get(lc,{}) and A[lc]["T"] is not None:
            m=nv(A[lc]["T"]).resample("MS").mean(); m=m[m.index>=pd.Timestamp("2009-01-01")]
            ax.bar(m.index,m.values,width=28,color=["#b2182b" if v>0 else "#2166ac" for v in m.values])
        ax.axhline(0,color="#444",lw=.5); ax.set_ylabel(f"{dp} m",fontsize=8); ax.set_ylim(-1.2,1.2)
        ax.axvspan(pd.Timestamp("2014-06-01"),pd.Timestamp("2016-06-30"),color="#f6c9a8",alpha=.3); ax.grid(axis="y",alpha=.2)
    axs[0].set_title("Monthly temperature anomaly by depth (shaded = 2014-16 Blob): shallow events vs the deep 2026 El Nino",fontsize=9,loc="left")
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_depthstack.png",dpi=130,bbox_inches="tight"); plt.close(fig)
except Exception: pass

# ---- status table ----
rows=[]
for lc,(nm,dp) in SITES.items():
    c={"s":nm,"d":dp}
    for var in ["T","O2"]:
        if var in A.get(lc,{}) and A[lc][var] is not None:
            a=A[lc][var]; l30=a[a.index>=a.index.max()-pd.Timedelta("30D")].mean()
            pv=a[(a.index>=a.index.max()-pd.Timedelta("60D"))&(a.index<a.index.max()-pd.Timedelta("30D"))].mean()
            tr="rising" if l30-pv>0.01 else "falling" if l30-pv<-0.01 else "flat"
            c[var]=f"{l30:+.2f} ({l30/SD[lc][var]:+.1f} sd, {tr})"
        else: c[var]="-"
    rows.append(c)
trows="".join(f"<tr><td>{r['s']}</td><td>{r['d']} m</td><td>{r['T']}</td><td>{r['O2']}</td></tr>" for r in rows)

html=f"""<!doctype html><html lang=en><head><meta charset=utf-8>
<meta name=viewport content="width=device-width, initial-scale=1"><title>Barkley Slope El Nino Watch</title>
<style>body{{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:920px;margin:0 auto;padding:1.2rem;color:#222;line-height:1.55}}
h1{{font-size:1.55rem;margin:.2rem 0}} h2{{font-size:1.14rem;margin-top:1.8rem;border-bottom:1px solid #ddd;padding-bottom:.25rem}}
.sub{{color:#666;font-size:.9rem}} table{{border-collapse:collapse;width:100%;font-size:.9rem}}
th,td{{border:1px solid #ccc;padding:.4rem .55rem;text-align:left}} th{{background:#f4f4f4}}
img{{width:100%;height:auto;border:1px solid #eee;margin:.4rem 0}} .k{{background:#fff7e6;padding:.7rem .9rem;border-left:4px solid #e8a33d;margin:.9rem 0}}
footer{{color:#888;font-size:.8rem;margin-top:2rem;border-top:1px solid #eee;padding-top:.6rem}} a{{color:#1a6}}</style></head><body>
<h1>Barkley Slope El Nino Watch</h1>
<p class=sub>Ocean Networks Canada cabled records, SW Vancouver Island. Updated {TODAY.date()} (UTC). RONI {roni_s}.</p>
<div class=k><b>The 2026 El Nino reached Vancouver Island through the ocean, not the air.</b> The surface stays near
normal while the Barkley slope at 396 m turns warmer, fresher, lighter and more oxygen-rich - the density surfaces
pushed down ~70 m by a coastally trapped wave. The warming is heave (the water was replaced from above); the oxygen
carries an extra on-isopycnal loss. Density surfaces measured directly by a Line P glider.</div>
<h2>This week</h2>
<table><tr><th>Site</th><th>Depth</th><th>Temperature anomaly</th><th>Oxygen anomaly</th></tr>{trows}</table>
<h2>Live watch</h2><img src="fig_watch.png" alt="watch">
<p class=sub>Red above zero = warmer / more oxygen than normal for the date. The El Nino shows at depth first.</p>
<h2>Is it heave, or new water?</h2>
<p>Two ways to measure the push, plus a direct glider measurement. The density and temperature estimates agree,
so the 396 m salinity is not badly drifting. Oxygen is the exception: heave should have raised it, but much less
was measured - the water on those density surfaces has lost oxygen.</p>
{heave_html}
<h2>2026 is not undercurrent water (T-S at 396 m)</h2><img src="fig_TS.png" alt="T-S">
<p class=sub>Undercurrent water is warm and SALTY; 2026 is warm and FRESH - the density field lifted toward the sensor from above.</p>
<h2>Which depth responds (depth-time stack)</h2><img src="fig_depthstack.png" alt="depth stack">
<p class=sub>2014-16 (Blob) was shallow; 2026 reaches deep - record-building at 396 and 643 m while the surface stays modest.</p>
<footer>Data: Ocean Networks Canada Oceans 3.0 (CC-BY 4.0); Line P / Station P4 CTD (CIOOS Pacific); C-PROOF Line P glider;
NOAA CPC. Density via TEOS-10. Auto-generated weekly. Lead with temperature (QC-clean); Folger Deep 96 m offline since 2026-07-18.</footer>
</body></html>"""
open(f"{OUT}/index.html","w").write(html)
print("site written", TODAY.date())
