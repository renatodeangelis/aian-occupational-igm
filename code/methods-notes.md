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

### 1.3 ~~Common support check refits a duplicate PS model~~ **FIXED**

`compute_weights()` now returns a named list (`$data`, `$p_hat_full`). `weighting.R`
unpacks this and uses the precomputed PS values for both the linked and unlinked samples
in the common support block. The duplicate `glm()` call has been removed.

---

## 2. Regional Analysis (`transition_matrices_weighted.R`)

### 2.1 Regional maps show no sample size or uncertainty

`om_1_plot`, `d_1_prime_plot`, and `upward_downward_plot` present point estimates for
all 7 regions with equal visual weight. Some regions may have thin linked AIAN samples,
making the 4×4 transition matrix unreliable.

**Fix**: Compute ESS per region (now available via the `n` column added by the
`group_modify` refactor) and either grey out regions below a threshold or overlay raw n.
The `results_region` data frame already carries `n`; join it into `centroids`:

```r
ess_by_region = data |>
  group_by(region) |>
  summarise(
    n_raw = n(),
    ess   = sum(w_atc_norm)^2 / sum(w_atc_norm^2),
    .groups = "drop"
  )
```

Consider suppressing estimates for regions with ESS < 30.

---

### 2.2 The upward/downward mobility ratio is not well defined

`p_upward` and `p_downward` are computed over different sub-populations of fathers,
so their denominators are not comparable. The ratio `p_upward / p_downward` cannot be
interpreted as "X upward movers per downward mover."

Additionally, the downward definition uses meso categories (`meso_pop == "farmer"`,
`meso_pop == "crafts"`) while the upward definition uses macro categories — the two
measures are not constructed symmetrically.

**Fix**: Either (a) define both measures over the same base population, (b) report the
two rates separately rather than as a ratio, or (c) replace with absolute upward
mobility: `P(macro_son == "nonmanual")` unconditional on father's class, which has a
clean interpretation and is directly comparable across regions.

---

## 4. Employment Classification (`cleaning-script.R`)

### 4.1 Variable overview

Four variables jointly bear on son's employment classification. Coding for 1940:

| Variable | Coding |
|---|---|
| `labforce_1940` | 1 = not in labor force; 2 = in labor force |
| `wkswork1_1940` | weeks worked in **prior year** (1939), 0–52 |
| `empstatd_1940` | 10 = at work; 11 = has job, not at work; 12 = has job, not at work (emergency); 21 = seeking work (unemployed) |
| `occ_son` (occ1950_1940) | occupation at enumeration (April 1940); > 970 = nilf/unemp |

Note: `wkswork1` measures 1939 activity; `occ` and `empstatd` measure April 1940 status.
These are different points in time and can legitimately conflict.

Note: `labforce` and `empstatd` are consistent for most observations. `labforce` is
a cruder summary; `empstatd` is more informative where they differ.

---

### 4.2 Current classification approach

`classify_meso(occ_son)` in `utils.R` is purely occ-first:
- `occ <= 970` → occupation category (farmer, farmworker, nonmanual, crafts, unskilled)
- `occ > 970` → nilf

`labforce_1940`, `wkswork1_1940`, and `empstatd_1940` are carried through to
`aian_merged.csv` but play no role in classification. `nilf` and `unemp` (occ 980–998)
are collapsed into a single category — intentional, as non-employment is kept as a
substantively meaningful state for this population.

**This is consistent with the occ-first conclusion below.**

---

### 4.3 Gray-area cases (aian_merged, n = 11,890)

Analysis run on the cleaned linked sample. Gray areas are cases where the three signals
conflict.

**Resolved cases — no classification change needed:**

| Pattern | n | Classification |
|---|---|---|
| in LF + valid occ + >0 wks | ~8,778 | Employ by occ |
| not in LF + occ > 970 + 0 wks | ~1,045 | nilf |
| in LF + occ 980–998 + 0 wks | 226 | nilf (collapsed) |
| all occ > 970, any LF/wks | ~1,583 | nilf |

**Gray areas and recommendations:**

| Pattern | n | Recommendation | Certainty |
|---|---|---|---|
| in LF + valid occ + 0 wks, empstatd 10/11/12 | 867 | Classify by occ — currently employed, 0 wks reflects recent entry | High |
| in LF + valid occ + 0 wks, empstatd 21 (seeking work) | 346 | Classify by occ (usual/last occ) — consistent with occ-first rule | Medium — see note |
| not in LF + valid occ + >0 wks | 212 | Classify by occ | Medium |
| not in LF + valid occ + 0 wks | 104 | Classify by occ | Low — most uncertain case |
| in LF + occ=999 + 0 wks | 86 | nilf | High |
| in LF + occ=999 + >0 wks | 50 | nilf (can't assign category) | High |
| not in LF + occ=999 + >0 wks | 109 | nilf (current status) | Medium |

**Note on empstatd=21 (346 seeking work, valid occ, 0 weeks):** These are currently
unemployed men whose occupation code reflects usual/last occupation. Classifying them
by occ treats occupational identity as persistent through unemployment spells, which
is the standard assumption in occupational mobility research. The alternative — sending
them to nilf — would mean unemployment and true non-employment are treated identically,
which is arguably a greater distortion given that `nilf` is intended as a meaningful
category here. However, this group is a candidate for the `wkswork_flag` robustness
check (see 4.5).

**Note on occ distribution of gray-area cases:** ~94% of gray-area cases fall in the
three bottom meso categories (farmworker, unskilled, farmer), which is consistent
with intermittent employment in casualized occupations.

---

### 4.4 Pending robustness flag

A `wkswork_flag` analogous to `spread_flag` has not yet been implemented. Candidate
definition: `empstatd_1940 == 21 & wkswork1_1940 == 0` (346 seeking-work cases with
0 prior-year weeks). This would allow re-running the main estimates excluding this
group without changing the primary classification. The 104 `not in LF + valid occ +
0 wks` cases could be included in the flag as a stricter variant.

---

### 4.5 Two-matrix proposal (occupational identity vs. labor market attachment)

**Concept:** Build two macro-level transition matrices using the same categories
(farming / manual / nonmanual / nilf):

- **Identity matrix** (current): occ-first, as implemented.
- **Attachment matrix**: men with a valid occ code but no labor market activity
  (best candidate rule: `empstatd == 21 & wkswork1 == 0`) reclassified to nilf.
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
   to be retained). Apply: valid occ + wkswork1 = 0 in that year → nilf in attachment
   matrix. `empstatd` equivalent may not be consistently available for fathers across
   census years; the father rule is therefore blunter than the son rule.

3. **Gate on empirical result:** implement both matrices and compare before committing
   to framing. The contested cases are ~10% of the sample and concentrated in the
   bottom categories. If the matrices are nearly identical, drop the comparison.
   If they differ in the nilf row/column in ways that speak to economic marginalization,
   develop the framing.

**Key risks:**

- Class-correlated reclassification: the attachment matrix drains farmworker/unskilled
  disproportionately (partly the point, but needs to be narrated carefully).
- Father-son asymmetry in classification rules (empstatd available for sons, not fathers)
  should be disclosed.
- Threshold choice (0 weeks + seeking work) should be fixed a priori and not varied
  to optimize results.

---
