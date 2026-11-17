# Barkley Slope El Nino Watch

An auto-updating page tracking how the 2026-27 El Nino develops in Ocean Networks Canada
cabled records on the southwest Vancouver Island slope. Complements Steve Mihaly's ONC
slope/canyon monitor by adding the mechanism, sea-level evidence, the 2023 analog and a forecast.

## Deploy (one time, ~5 minutes)

1. Create a new GitHub repository (public), e.g. `barkley-el-nino`.
2. Upload the contents of this folder (keep the structure: `generate_site.py`,
   `docs/`, `.github/workflows/update.yml`, this README).
3. In the repo: **Settings -> Secrets and variables -> Actions -> New repository secret**.
   Name it `ONC_TOKEN`, paste your Ocean Networks Canada API token.
4. In the repo: **Settings -> Pages**. Source: "Deploy from a branch", branch `main`, folder `/docs`.
   Your page appears at `https://<your-username>.github.io/barkley-el-nino/`.
5. The page rebuilds every Monday 16:00 UTC automatically. To rebuild now:
   **Actions -> weekly-update -> Run workflow**.

## What it does

`generate_site.py` fetches the latest ONC CTD + oxygen (Folger 25 m, Upper Slope 396 m,
Node 643 m, Canyon Axis 985 m) and NOAA ENSO / sea-level data, computes anomalies against a
pre-2026 harmonic climatology, and writes `docs/index.html`, `docs/fig_watch.png`,
`docs/fig_story.png` and `docs/status.csv`.

The figures in `docs/` shipped with this bundle are the 7 November 2026 snapshot; the first
Action run replaces them with live data.

## Notes
- The ONC fetch of ~18 years for four sites takes several minutes; that is normal in the Action.
- Lead with temperature (it passes instrument-swap QC); 396 m salinity may carry sensor drift.
- Data: Ocean Networks Canada Oceans 3.0 (CC-BY 4.0); NOAA CPC/PSL; NOAA CO-OPS.
