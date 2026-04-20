# Methods Notes: Outstanding Issues and Solutions

This document records methodological concerns identified during code review, along with
proposed fixes where they exist. Organized by script/stage of the pipeline.

---

## 1. Weighting (`weighting.R`)

### 1.1 Bootstrap does not account for weight estimation uncertainty

**This is the most consequential methodological issue in the paper.**

All bootstrap inference in `transition_matrices_weighted.R` resamples observations
but reuses the weights estimated on the full sample. This treats the propensity scores
as known population quantities rather than estimates, understating all standard errors
and confidence intervals. The bias is worst for regional estimates (small n, high
weight variance) and for any statistic sensitive to extreme weights.

**Partial progress**: `compute_weights()` is now in `utils.R` and takes `df_linked`
and `df_full` as arguments, making it callable on bootstrap resamples. But the actual
bootstrap functions (`boot_measures_by_t`, `boot_pmatrix_ci`, `boot_im_by_t`,
`mobility_curve_with_boot`) still resample the pre-weighted `data` object and do not
call `compute_weights()` on each draw. The infrastructure is ready; the wiring is not.

**Fix**:

Pass `aian_merged` (the unweighted linked sample) and `aian_full` through to the
bootstrap loops and call `compute_weights(df_linked[idx,], df_full)` on each resample:

```r
boot_with_reweighting = function(df_linked, df_full, stat_fn, R = 500, .seed = NULL) {
  if (!is.null(.seed)) set.seed(.seed)
  N = nrow(df_linked)

  est = stat_fn(compute_weights(df_linked, df_full))

  boots = replicate(R, {
    idx  = sample.int(N, N, replace = TRUE)
    d_b  = compute_weights(df_linked[idx, ], df_full)
    stat_fn(d_b)
  })

  list(estimate = est, se = sd(boots, na.rm = TRUE), draws = boots)
}
```

**Practical note**: R=200–500 draws is usually sufficient for SE estimation.
The full weight bootstrap should be run for final reported estimates; the fast
fixed-weight bootstrap is acceptable for exploratory analysis.

---

### 1.2 Covariate balance diagnostics are not reported in the paper

`weighting.R` now saves `output/balance_table.csv` (cohort, region, education, urban)
and `output/balance_table_state.csv` after weighting. Remaining decision: whether to
include a balance table in the appendix or cite the max SMD inline as a footnote.

---

## 2. Regional Analysis (`transition_matrices_weighted.R`)

### 2.1 Regional maps should flag California's smaller sample

The analysis uses 7 regions. Six regions have n > 1,200, which is sufficient to
estimate a 4×4 macro transition matrix with reasonable precision. California has
n = 688 and is retained for substantive reasons (distinct labor market and policy
context for AIAN men in this period). All regional figures should visually flag
California — either with a lighter fill, a border marker, or an asterisk — so readers
can discount its estimates appropriately without suppressing them.

ESS should still be computed per region (as noted in the previous version of this
section) and reported in a supplement or footnote. The `results_region` data frame
already carries `n`; ESS can be joined via:

```r
ess_by_region = data |>
  group_by(region) |>
  summarise(
    n_raw = n(),
    ess   = sum(w_atc_norm)^2 / sum(w_atc_norm^2),
    .groups = "drop"
  )
```

---

### 2.2 Regional summary statistics: replace upward/downward ratio

The current `p_upward / p_downward` ratio is not interpretable: numerator and
denominator are computed over different base populations, and the two measures mix
meso and macro categories asymmetrically. The ratio should be dropped entirely.

**Recommended replacement: five statistics at the macro level.**

Use macro categories (farming, manual, nonmanual, nonemp) for all regional map
statistics. Reserve meso categories for the global transition matrix heatmaps and
IM curves, where the richer classification adds value without the cross-region
comparability problem.

The five regional statistics, and their logic:

| Statistic | Interpretation |
|---|---|
| **SM** (structural mobility) | Share of observed movement forced by shifts in the marginal occupational distribution. High SM means sons couldn't stay in their fathers' occupations because those occupations were disappearing — the direct correlate of land dispossession and reservation confinement. |
| **EM** (exchange mobility) | Share of observed movement representing genuine rank-switching, net of structural change. EM and SM together decompose OM: if EM is near zero and SM is high, the paper's argument is that observed mobility is forced rather than achieved. |
| **P(son = manual \| father = farming)** | The proletarianization rate: farming families whose sons ended up in manual wage labor. The headline cell for the paper's structural narrative. Well-estimated with n > 1,200. |
| **P(son = farming \| father = farming)** | Farming persistence. Regional variation in this cell maps onto variation in land loss intensity under allotment. Paired with the preceding statistic, it shows what happened to farming families: some reproduced, most proletarianized. |
| **P(son = nonemp \| father = farming)** | Labor force exit from the farming sector. Provides the third branch of the farming-family story: exit into nonemployment rather than wage labor. |

Statistics 3–5 together give a complete decomposition of farming-family outcomes
across regions. Regional variation in their relative magnitudes is the core empirical
contribution of the regional analysis.

**Honorable mention**: P(son ≠ nonemp | father = nonemp) — exit from nonemployment — is
well-defined and has adequate n, but covaries strongly with the nonemp exit rate above.
Include in a regional table; probably not on the map.

**Excluded statistics and reasons**:

- *OM alone*: conflates direction and magnitude; sensitive to marginal distributions
  in ways that make raw cross-region comparison misleading. Report only as the sum
  EM + SM.
- *d'(1) on regional maps*: see issue 2.3 below.
- *Upward/downward rates (meso)*: asymmetric base populations; ordinality between
  farmer and farmworker is contested; equivalent story told more cleanly by the macro
  statistics above.
- *IM(t, i) curves by region*: still feasible with n > 1,200 for major origin
  categories; consider including in an appendix for completeness.


## 3. Employment Classification (`cleaning-script.R`)

### 3.1 Variable overview

Four variables jointly bear on son's employment classification. Coding for 1940:

| Variable | Coding |
|---|---|
| `labforce_1940` | 1 = not in labor force; 2 = in labor force |
| `wkswork1_1940` | weeks worked in **prior year** (1939), 0–52 |
| `empstatd_1940` | 10 = at work; 11 = has job, not at work; 12 = has job, not at work (emergency); 21 = seeking work (unemployed) |
| `occ_son` (occ1950_1940) | occupation at enumeration (April 1940); > 970 = nonemp/unemp |

Note: `wkswork1` measures 1939 activity; `occ` and `empstatd` measure April 1940 status.
These are different points in time and can legitimately conflict.

Note: `labforce` and `empstatd` are consistent for most observations. `labforce` is
a cruder summary; `empstatd` is more informative where they differ.

---

### 3.2 Current classification approach

`classify_meso(occ_son)` in `utils.R` is purely occ-first:
- `occ <= 970` → occupation category (farmer, farmworker, nonmanual, crafts, unskilled)
- `occ > 970` → nonemp

`labforce_1940`, `wkswork1_1940`, and `empstatd_1940` are carried through to
`aian_merged.csv` but play no role in classification. `nonemp` and `unemp` (occ 980–998)
are collapsed into a single category — intentional, as non-employment is kept as a
substantively meaningful state for this population.

**This is consistent with the occ-first conclusion below.**

---

### 3.3 Gray-area cases (aian_merged, n = 11,890)

Analysis run on the cleaned linked sample. Gray areas are cases where the three signals
conflict.

**Resolved cases — no classification change needed:**

| Pattern | n | Classification |
|---|---|---|
| in LF + valid occ + >0 wks | ~8,778 | Employ by occ |
| not in LF + occ > 970 + 0 wks | ~1,045 | nonemp |
| in LF + occ 980–998 + 0 wks | 226 | nonemp (collapsed) |
| all occ > 970, any LF/wks | ~1,583 | nonemp |

**Gray areas and recommendations:**

| Pattern | n | Recommendation | Certainty |
|---|---|---|---|
| in LF + valid occ + 0 wks, empstatd 10/11/12 | 867 | Classify by occ — currently employed, 0 wks reflects recent entry | High |
| in LF + valid occ + 0 wks, empstatd 21 (seeking work) | 346 | Classify by occ (usual/last occ) — consistent with occ-first rule | Medium — see note |
| not in LF + valid occ + >0 wks | 212 | Classify by occ | Medium |
| not in LF + valid occ + 0 wks | 104 | Classify by occ | Low — most uncertain case |
| in LF + occ=999 + 0 wks | 86 | nonemp | High |
| in LF + occ=999 + >0 wks | 50 | nonemp (can't assign category) | High |
| not in LF + occ=999 + >0 wks | 109 | nonemp (current status) | Medium |

**Note on empstatd=21 (346 seeking work, valid occ, 0 weeks):** These are currently
unemployed men whose occupation code reflects usual/last occupation. Classifying them
by occ treats occupational identity as persistent through unemployment spells, which
is the standard assumption in occupational mobility research. The alternative — sending
them to nonemp — would mean unemployment and true non-employment are treated identically,
which is arguably a greater distortion given that `nonemp` is intended as a meaningful
category here. However, this group is a candidate for the `wkswork_flag` robustness
check (see 3.5).

**Note on occ distribution of gray-area cases:** ~94% of gray-area cases fall in the
three bottom meso categories (farmworker, unskilled, farmer), which is consistent
with intermittent employment in casualized occupations.

---

### 3.4 Pending robustness flag

A `wkswork_flag` analogous to `spread_flag` has not yet been implemented. Candidate
definition: `empstatd_1940 == 21 & wkswork1_1940 == 0` (346 seeking-work cases with
0 prior-year weeks). This would allow re-running the main estimates excluding this
group without changing the primary classification. The 104 `not in LF + valid occ +
0 wks` cases could be included in the flag as a stricter variant.

---

### 3.5 Two-matrix proposal (occupational identity vs. labor market attachment)

**Concept:** Build two macro-level transition matrices using the same categories
(farming / manual / nonmanual / nonemp):

- **Identity matrix** (current): occ-first, as implemented.
- **Attachment matrix**: men with a valid occ code but no labor market activity
  (best candidate rule: `empstatd == 21 & wkswork1 == 0`) reclassified to nonemp.
  Applied symmetrically to fathers and sons.

Framing: not a robustness check, but a secondary analysis examining whether occupational
persistence holds under a stricter definition of labor market participation. Motivated
by the argument that for AIAN men at this period, formal labor market attachment was
inconsistent and occupational identity may not have reflected active employment.

**Prerequisites before pursuing:**

1. **Verify wkswork1_pop_* columns exist in the raw extract.** The father-side
   classification requires wkswork1 data for fathers. The original exclusion
   (`starts_with("wkswork1")`) was removed, so these columns are now in `aian_clean`
   if the extract includes them. Check before investing further.

2. **Father-side rule:** use `wkswork1` from the same census year as the picked
   occupation (`picked_year`, currently dropped from `modal_occ_pick` — would need
   to be retained). Apply: valid occ + wkswork1 = 0 in that year → nonemp in attachment
   matrix. `empstatd` equivalent may not be consistently available for fathers across
   census years; the father rule is therefore blunter than the son rule.

3. **Gate on empirical result:** implement both matrices and compare before committing
   to framing. The contested cases are ~10% of the sample and concentrated in the
   bottom categories. If the matrices are nearly identical, drop the comparison.
   If they differ in the nonemp row/column in ways that speak to economic marginalization,
   develop the framing.

**Key risks:**

- Class-correlated reclassification: the attachment matrix drains farmworker/unskilled
  disproportionately (partly the point, but needs to be narrated carefully).
- Father-son asymmetry in classification rules (empstatd available for sons, not fathers)
  should be disclosed.
- Threshold choice (0 weeks + seeking work) should be fixed a priori and not varied
  to optimize results.

---

## 4. Paper Framing and Scope

### 4.1 Core framing decision

This is a descriptive baseline paper. Its contribution is to produce the first
systematic quantitative evidence on the structure and magnitude of intergenerational
occupational mobility among AIAN men in the late 19th–early 20th century. The correct
framing is: "We document the occupational structure and its intergenerational
transmission for a population whose economic history has been studied qualitatively
but not quantified at this scale."

The paper should not present itself as identifying causes of mobility or making
causal comparisons across groups. The Markov chain machinery is a descriptive tool,
not a structural model.

---

### 4.2 Paper structure

**Main text:**

1. Weighted macro transition matrices (P, pi_0, pi*) — the primary descriptive object
2. EM/SM decomposition — the main analytical contribution; the ratio SM/(SM+EM) is
   the paper's headline finding
3. Regional analysis: five statistics (SM, EM, P(manual|farming), P(farming|farming),
   P(nonemp|farming)) as maps and a regional table
4. Cohort stationarity results — brief table, not a supplement

**Appendix:**

- Meso-level transition matrices (support the macro results)
- d(t), d'(t), IM curves (characterize the Markov chain mathematically; informative
  but not the story)
- Log-linear UNIDIFF analysis by region (tests whether regional variation is a matter
  of degree vs. kind; use macro 4×4 to avoid sparse cells)
- Robustness checks: spread_flag exclusions, wkswork_flag (once implemented), 25–44
  age restriction
- Balance table from `output/balance_table.csv`

**Cut or defer:**

- Multi-generation projections as a narrative device — one sentence on the implied
  steady state is enough
- Two-matrix proposal (section 3.5) — defer to future work; do not add to this paper

---

### 4.3 Abstract and introduction framing

The abstract should open with the gap and the data, not the methods:

> "No systematic quantitative evidence exists on intergenerational occupational
> mobility for Native American men in the early 20th century. We use linked 1900–1940
> U.S. Census records to document the structure of father-son occupational mobility
> for this population. We find that observed mobility was overwhelmingly structural —
> driven by the collapse of the farming sector — rather than achieved through
> individual rank-switching, consistent with the historical record of allotment-era
> dispossession."

The introduction should spend more words on historical context (allotment, land
dispossession, reservation confinement) and less on the Markov chain framework.
Readers need to understand why this period and this population matter before
encountering the first equation.

Three limitations should be named clearly in the introduction, once, and not repeated:

1. **Linkage selection**: estimates represent men who can be linked; if geographic
   mobility and occupational mobility are positively correlated, these are lower
   bounds on true mobility.
2. **No within-sample comparison group**: magnitudes are contextualised via parallel
   work on other populations.
3. **Occupation measurement**: occ1950 codes applied to pre-1940 AIAN labor markets
   introduce measurement noise, particularly for casual and subsistence work.

---

### 4.4 Methods section framing

Order: data → weighting → mobility measures. Do not lead with the Markov chain
framing; lead with the transition matrix as a descriptive object.

Present PS weighting as a correction for known selection, not as a technical
contribution. One paragraph on the selection problem, one on ATC weighting as the
solution for observable characteristics, one sentence acknowledging unobservable
selection remains and noting the likely direction of bias.

Retire "Markov chain" from the front half of the paper. Use "weighted transition
matrix" as the primary term; the Markov interpretation can appear when discussing
the convergence measures in the appendix.

---

### 4.5 Results section framing

Lead with the transition matrix figure — give the reader the full picture first.
Then decompose into EM and SM. Narrative structure:

1. Here is the overall pattern (macro P matrix + pi_0 + pi*)
2. Here is how much of it is structural displacement vs. rank switching (EM/SM)
3. Here is how it varies regionally (maps + table)

Do not present the regional matrices in the main text (7 × 4×4 = 112 cells; no
reader will absorb these inline). The regional section is maps and a summary table.

---

### 4.6 Code changes required for reframing

**`compute_mobility_stats()` in `transition_matrices_weighted.R`** currently computes
`p_upward`, `p_downward`, `d_prime_1`, `om_weighted`, and `om_unweighted`. Replace
with SM, EM, and the three conditional probabilities:

```r
compute_mobility_stats = function(df) {
  df = df |> mutate(w_atc_norm = w_atc_norm / sum(w_atc_norm) * n())

  P      = p_matrix(df, macro_pop, macro_son, weighted = TRUE)
  pi0    = pi_0(df, macro_pop)
  pistar = pi_star(P)

  om_val = om(P, pi0, t = 0)
  sm_val = sum(abs(as.numeric(pi0) - as.numeric(pistar))) / 2
  em_val = om_val - sm_val

  # Farming-family decomposition
  farming_rows = df |> filter(macro_pop == "farming")
  p_manual_given_farming = weighted.mean(
    farming_rows$macro_son == "manual", farming_rows$w_atc_norm, na.rm = TRUE)
  p_farming_given_farming = weighted.mean(
    farming_rows$macro_son == "farming", farming_rows$w_atc_norm, na.rm = TRUE)
  p_nonemp_given_farming = weighted.mean(
    farming_rows$macro_son == "nonemp", farming_rows$w_atc_norm, na.rm = TRUE)

  tibble(
    sm                   = round(sm_val, 3),
    em                   = round(em_val, 3),
    om                   = round(om_val, 3),
    p_manual_fm_farming  = round(p_manual_given_farming, 3),
    p_farming_fm_farming = round(p_farming_given_farming, 3),
    p_nonemp_fm_farming    = round(p_nonemp_given_farming, 3)
  )
}
```

Note: `p_manual_fm_farming + p_farming_fm_farming + p_nonemp_fm_farming + p_nonmanual_fm_farming = 1`
by construction; the three statistics above plus the (small) nonmanual rate fully
decompose farming-father outcomes.

**Regional maps**: replace the three current map calls (`om_1_plot`, `d_1_prime_plot`,
`upward_downward_plot`) with maps of SM, EM, and the three conditional probabilities.

**Cohort stationarity (section 7 of the script)**: move results to be printed as a
formatted table for inclusion in the main text, not just `cat()` output.

---

### 4.7 Birth cohort treatment

The 20–44 age range for sons is retained as the main analysis. Three reasons:

1. The dominant occupation categories (farmer, farmworker, unskilled, nonemp) are
   not subject to the early-career instability that motivates age restrictions in
   white-collar samples. A 22-year-old AIAN farmworker in 1940 is in a stable
   occupational category, not a career placeholder.
2. Restricting to 25–44 loses ~44% of the sample (12,000 → 6,700) and would require
   pooling regions further, reducing the regional analysis.
3. The cohort stationarity test is the appropriate tool for addressing the reviewer
   concern about life cycle and cohort effects — it tests directly whether pooling
   across birth cohorts is warranted.

**Robustness check**: run the main transition matrices on the 25–44 subsample and
report the comparison in a footnote or appendix table. If results are substantively
similar, the age restriction does not change the findings and the full sample is
preferable on power grounds. Report the direction of any differences.

The cohort stationarity section (section 7 of `transition_matrices_weighted.R`)
should be elevated to the main text as a brief table showing OM(1) and d'(1) by
birth cohort, with a one-sentence interpretation of whether the pooled P matrix
is warranted.

---

## 5. Work Plan

### 5.1 Prerequisites

Confirm `data/aian_merged.csv` is stale (cleaning changes not yet rerun) and that
R is pointed at the project root. Scripts must run in pipeline order:
`cleaning-script.R` → `weighting.R` → `transition_matrices_weighted.R`.

---

### 5.2 Block 1 — Run cleaning and read diagnostics (45 min)

**Goal:** Verify new cleaning changes (MAD spread, son-age tiebreaker) and understand
where observations are lost.

1. `source("code/cleaning-script.R")` — read console output carefully
2. Record the six counts: raw → age/sex → dedup → occ_pop recovery → spread filter
   → school filter → final
3. Check: is MAD ≤ 4 dropping more or fewer than the old range ≤ 10? If
   substantially more, adjust threshold
4. Spot check: `table(aian_merged$meso_pop, aian_merged$meso_son)` for plausibility
5. Bonus if time: `table(aian_merged$occ_son[aian_merged$meso_son == "nonemp"])` to
   assess occ == 999 contamination in the nonemp row (see issue 3.2)

If cleaning looks clean, run `source("code/weighting.R")`.

---

### 5.3 Block 2 — Fix `compute_mobility_stats()` (1.5 hours)

**Goal:** Replace broken upward/downward logic with SM, EM, and three conditional
probabilities (see section 4.6 for full replacement code).

After implementing, verify:
- `sm + em ≈ om` for each region
- Three conditional probabilities sum to < 1 (nonmanual is the fourth branch)
- `results_region` carries: `sm`, `em`, `om`, `p_manual_fm_farming`,
  `p_farming_fm_farming`, `p_nonemp_fm_farming`

---

### 5.4 Block 3 — Update regional maps (45 min)

**Goal:** Replace the three old map calls with five new statistics; flag California.

1. Replace `om_1_plot`, `d_1_prime_plot`, `upward_downward_plot` with five new map
   calls using updated `results_region` columns
2. Flag California: add `low_n = (n < 800)` to `centroids` and append an asterisk
   to the label via `geom_sf_text`
3. Visually inspect: SM should be high in Oklahoma/Plains where allotment was most
   aggressive; farming persistence should be higher in regions with stronger land
   retention

---

### 5.5 Block 4 — Wire bootstrap reweighting (2 hours)

**Goal:** Fix issue 1.1. Infrastructure exists in `utils.R`; wiring is missing.

Scope: wire `boot_with_reweighting()` (template in section 1.1) into
`boot_pmatrix_ci` first as proof of concept, then `boot_measures_by_t`. The other
two bootstrap functions (`boot_im_by_t`, `mobility_curve_with_boot`) can follow the
same pattern in a second session if time runs short.

Requires `aian_merged` (unweighted) and `aian_full` (full AIAN extract from
`weighting.R`) to be available in the environment. Verify `aian_full` persists after
sourcing `weighting.R`, or save it to a file.

**Highest risk of overrunning.** If blocked past one hour, note the issue and move on.

---

### 5.6 Block 5 — Cohort stationarity table (30 min)

**Goal:** Elevate section 7 of `transition_matrices_weighted.R` from `cat()` output
to a formatted table for the main text.

Replace the loop with `purrr::map_dfr` returning a tidy data frame. Table should
show: cohort, n, ESS, OM(1), d'(1) — four columns, three rows.

---

### 5.7 Block 6 — 2×2 LF participation matrix (45 min)

**Goal:** Supplementary analysis of intergenerational LF participation persistence
(motivated by nonemp→nonemp measurement issues; see section 3.2).

1. Confirm `labforce_pop_*` columns are in `aian_merged`
2. Define father's LF status from the labforce value in the census year matching
   `occ_pop` (same year as the picked occupation)
3. Build weighted 2×2: father LF (1/2) × son LF (1/2)
4. Output: four cells with bootstrap CIs — a table, not a heatmap

---

### 5.8 End of day — Full pipeline run (30 min)

Source all three scripts in order. Verify outputs exist and all figures render.
Note errors or unexpected results for the next session.

---

### 5.9 Not covered in this session

- Log-linear UNIDIFF (appendix — second session)
- `wkswork_flag` implementation (section 3.4 — straightforward, not urgent)
- 25–44 age robustness check (one-liner once main analysis is stable; see section 4.7)
- Full bootstrap reweighting for `boot_im_by_t` and `mobility_curve_with_boot`
  (continue from block 4)

---

## 6. Remaining Action Items (from April 2026 session)

This section records items identified during the attachment-variable construction session
that were discussed but not yet implemented.

---

### 6.2 `weighting.R`: Pass alt columns through and save `aian_full`

Two tasks:

1. **Verify alt columns pass through.** Confirm that `occ_pop_alt`, `meso_pop_alt`,
   `macro_pop_alt`, `occ_son_alt`, `meso_son_alt`, `macro_son_alt`, `attachment_level_son`,
   `attachment_level_pop`, `picked_year_alt` survive into `aian_weighted.csv`. These
   columns were added in `cleaning-script.R` but `weighting.R` may have a `select()`
   that drops them. Check and fix any such drops.

2. **Save `aian_full` as an RDS file.** Add at the end of `weighting.R`:

```r
saveRDS(aian_full, "data/aian_full.rds")
```

This is required by the bootstrap reweighting fix (section 1.1 / item 6.4 below):
`transition_matrices_weighted.R` needs access to the full unlinked AIAN extract
to re-estimate PS weights on each bootstrap resample. Currently `aian_full` is
created inside `weighting.R` but not persisted.

---

### 6.3 `transition_matrices_weighted.R`: Alt matrix analysis

After alt columns pass through weighting (item 6.2), implement a side-by-side
comparison of the main matrix and the attachment matrix.

**Steps:**

1. After loading `aian_weighted.csv`, re-apply factor levels to alt columns:
```r
data = data |>
  mutate(
    macro_pop_alt = factor(macro_pop_alt, levels = macro_order, ordered = TRUE),
    macro_son_alt = factor(macro_son_alt, levels = macro_order, ordered = TRUE),
    meso_pop_alt  = factor(meso_pop_alt,  levels = meso_order,  ordered = TRUE),
    meso_son_alt  = factor(meso_son_alt,  levels = meso_order,  ordered = TRUE)
  )
```

2. Compute the alt macro transition matrix:
```r
P_alt = p_matrix(data, macro_pop_alt, macro_son_alt, weighted = TRUE)
```

3. Run `compute_mobility_stats` on both the main and alt datasets (using `occ_pop_alt`
   / `occ_son_alt` columns to construct the alt version), and print the side-by-side
   comparison.

4. Run `boot_pmatrix_ci` on `macro_pop_alt` / `macro_son_alt` once bootstrap
   reweighting (item 6.4) is wired.

**Decision gate:** if the main and alt matrices are nearly identical, drop the alt
analysis and note in the paper. If they differ substantially in the nonemp
row/column, develop the framing as a secondary analysis (not a robustness check).

---

### 6.4 Bootstrap reweighting: full wiring (issue 1.1)

Infrastructure is in place (`compute_weights()` in `utils.R`). Three changes needed:

**Step A — `utils.R`: add `boot_with_reweighting()`**

The template is in section 1.1. Add it to `utils.R` after `compute_weights()`:

```r
boot_with_reweighting = function(df_linked, df_full, stat_fn, R = 500, .seed = NULL) {
  if (!is.null(.seed)) set.seed(.seed)
  N   = nrow(df_linked)
  est = stat_fn(compute_weights(df_linked, df_full))
  boots = replicate(R, {
    idx = sample.int(N, N, replace = TRUE)
    d_b = compute_weights(df_linked[idx, ], df_full)
    stat_fn(d_b)
  })
  list(estimate = est, se = sd(boots, na.rm = TRUE), draws = boots)
}
```

**Step B — `utils.R`: update `boot_pmatrix_ci()`**

Add optional `df_linked` / `df_full` arguments. When provided, call
`compute_weights(df_linked[idx,], df_full)` on each draw rather than
resampling the pre-weighted `data` object directly.

**Step C — `transition_matrices_weighted.R`: wire inputs and calls**

At the top of the script, after loading `data`:
```r
aian_full        = readRDS("data/aian_full.rds")
aian_merged_raw  = read_csv("data/aian_merged.csv")  # unweighted linked sample
```

Update the two `boot_pmatrix_ci` calls (lines ~223–224) to pass
`df_linked = aian_merged_raw, df_full = aian_full`.

Add a `boot_em_sm()` call using `boot_with_reweighting` with
`compute_mobility_stats` as the stat function, for global EM/SM CIs.

**Scope boundary:** do NOT update `boot_measures_by_t`, `boot_im_by_t`, or
`mobility_curve_with_boot` in the same pass — these are memory-measure functions
that will require additional refactoring. Handle in a separate session.
