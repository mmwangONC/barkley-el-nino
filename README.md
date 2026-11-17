# Barkley El Niño Tracker

A weekly-updated page tracking how the 2026–27 El Niño shows up on the continental
slope off southwest Vancouver Island, using Ocean Networks Canada (ONC) cabled
observatory data.

**Live page:** https://mmwangonc.github.io/barkley-el-nino/

## What this shows

The 2026–27 El Niño is one of the strongest on record in the tropical Pacific. Off
Vancouver Island, its signal has so far appeared at depth, not at the surface:

- On the Barkley Canyon upper slope (396 m), temperature, salinity, density and
  dissolved oxygen anomalies point to density surfaces being pushed down.
- Near-surface temperature (ONC Folger Pinnacle, 25 m; NOAA Neah Bay buoy) stayed
  near normal or cold through autumn 2026.

The page plots daily and monthly anomalies at several depths against each site's
own climatology, alongside the NOAA Oceanic Niño Index.

## Sites

| Site | Code | Depth |
|---|---|---|
| Folger Pinnacle | FGPPN | 25 m |
| Folger Deep | FGPD | 96 m |
| Barkley Upper Slope | NCBC | 396 m |
| Barkley Node | BACND | 643 m |
| Barkley Canyon Axis | BACAX | 985 m |

## How it updates

A GitHub Actions workflow (`.github/workflows/weekly-update.yml`) runs once a week.
It pulls the latest data from the ONC Oceans 3.0 API, recomputes anomalies, and
rebuilds the page in `/docs`. To force an update: **Actions → weekly-update → Run
workflow**.

## Setup (for a copy of this repo)

1. Add your ONC API token as a repository secret named `ONC_TOKEN`
   (Settings → Secrets and variables → Actions).
2. Allow the workflow to write: Settings → Actions → General → Workflow
   permissions → *Read and write*.
3. Enable Pages: Settings → Pages → Deploy from branch `main`, folder `/docs`.

Get an ONC token from your Oceans 3.0 profile (Web Services API tab).

## Caveats

- **Preliminary.** Recent data have passed ONC's automated QC but not full manual
  review. The latest few weeks can change.
- Anomalies are relative to each site's pre-2026 record (8–17 years depending on
  the site), so record lengths differ between depths.
- Partial months are provisional.
- Salinity and oxygen sensor offsets between deployments are still being checked.
- Folger Deep (96 m) has had no data since July 2026.

## Data sources

- Ocean Networks Canada, Oceans 3.0: https://data.oceannetworks.ca
- NOAA Climate Prediction Center, Oceanic Niño Index
- NOAA National Data Buoy Center, station 46087 (Neah Bay)

Please cite Ocean Networks Canada when using these data, following the ONC data
policy.

## Contact

Manman, Ocean Networks Canada


