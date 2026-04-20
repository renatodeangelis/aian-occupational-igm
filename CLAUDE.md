# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based empirical research project analyzing intergenerational occupational mobility (IGM) of Native American (AIAN) men using linked U.S. Census data from 1900–1940. It produces a paper estimating weighted transition matrices and mobility statistics.

## Running the Code

Scripts must be run sequentially — each produces output consumed by the next:

```r
# Step 1: Clean and merge raw census data
source("code/cleaning-script.R")   # → data/aian_merged.csv

# Step 2: Estimate propensity score weights
source("code/weighting.R")          # → data/aian_weighted.csv

# Step 3: Compute transition matrices and generate figures
source("code/transition_matrices_weighted.R")
```

No build automation exists. Run scripts from RStudio or via `Rscript` from the repo root.

## Architecture

**Sequential ETL pipeline:**

1. **`cleaning-script.R`** — Loads raw linked father-son census pairs from Dropbox, deduplicate father records, applies age filters, classifies occupations via `classify_meso()` and `classify_macro()`, outputs `aian_merged.csv`.

2. **`weighting.R`** — Loads `aian_merged.csv` plus a full AIAN census extract, fits a logistic PS model (`linked ~ birthyr_son + region + education`), computes ATC weights, outputs `aian_weighted.csv`.

3. **`transition_matrices_weighted.R`** — The main analysis script (1100+ lines). Computes weighted 4×4 transition matrices (meso → meso), mobility measures (OM, EM, SM, upward/downward ratios, convergence to stationary distribution), performs bootstrap inference, and generates all publication figures.

**Occupation classification hierarchy:**
- Raw `occ1950` codes → **meso** categories (in `cleaning-script.R` / `weighting.R`): `farmer`, `farmworker`, `prof`, `clerical`, `crafts`, `unskilled`, `unemp`
- `transition_matrices_weighted.R` collapses `prof` + `clerical` → `white_collar` at load time, yielding 6 meso categories: `farmer`, `farmworker`, `white_collar`, `crafts`, `unskilled`, `unemp`
- Meso → **macro** categories: `farming`, `manual`, `nonmanual`, `unemp`

**Regional stratification:** 12 U.S. regions (basin, cali, lakes, midwest, nc, ne, nw, ok, plains, prairie, south, sw).

## General Conduct

- Do not offer historical or factual claims as supporting context unless they can be verified from the code, data, or user-provided materials — if uncertain, say so rather than confabulating

## Known Methodological Issues

`code/methods-notes.md` documents 20+ outstanding issues. The most consequential:

- **2.5 (critical)**: Bootstrap reuses weights estimated on the full sample — SEs are understated. Fix requires re-estimating PS model on each bootstrap resample.
- **3.1**: Within-region weights are not renormalized before regional statistics are computed.

Before modifying occupation classification logic, check both `cleaning-script.R` (lines 119–141) and `weighting.R` (lines 29–51) — they must stay in sync until issue 1.8 is resolved.

## Data

- `data/aian_merged.csv` — Cleaned father-son linked pairs (output of step 1)
- `data/aian_weighted.csv` — Weighted analysis dataset (output of step 2)
- `data/res_counties.csv` — Reservation county reference for regional classification
- Raw data is loaded from a Dropbox URL in `cleaning-script.R` and is not versioned here
