# Methods Notes: Outstanding Issues and Solutions

This document records methodological concerns identified during code review, along with
proposed fixes where they exist. Organized by script/stage of the pipeline.

---

## 2. Weighting (`weighting.R`)

### 2.3 Birth year specification is overly flexible

Birth year is entered as a factor with up to 25 levels (born 1896–1920). This eats
degrees of freedom and doesn't exploit the fact that adjacent birth years likely have
similar linkage propensities. The factor specification is unusual in the census-linking
literature and probably overfitting.

| Specification | Pros | Cons |
|---|---|---|
| 5-year cohort bins | Smooth, parsimonious, interpretable | Arbitrary bin edges |
| Quadratic polynomial | 2 df, captures nonlinear age effects | Misses sharp changes |
| Natural spline (`ns(birthyr, df=4)`) | Flexible, smooth, few df | Harder to interpret |
| Factor (current) | Maximally flexible | 24 df, overfits in small regions |

**Fix**: Use 5-year cohort bins (standard in Abramitzky, Boustan, & Eriksson; Bailey
et al.) or a natural spline with 3–4 df:

```r
aian_comb = aian_comb |>
  mutate(cohort = cut(as.numeric(as.character(birthyr_son)),
                      breaks = c(1895, 1900, 1905, 1910, 1915, 1921),
                      labels = c("1896-1900","1901-1905","1906-1910",
                                 "1911-1915","1916-1920")))
```

---

### 2.4 Missing covariates that predict linkage

Several variables available in the 1940 IPUMS census are absent from the PS model but
plausibly predict linkage success, especially for AIAN populations:

| Variable | IPUMS name | Why it predicts linkage |
|---|---|---|
| Urbanization | `urban` / `metro` | Urban AIAN men more mobile, harder to link |
| Literacy | `lit` | Literate men have more consistent records |
| Household size | `famsize` / `nchild` | Larger households harder to disambiguate |
| Marital status | `marst` | Married men more residentially stable |
| Employment status | `empstat` | Employed men have recorded occupations |
| Birthplace | `bpl` | Reservation-born may link differently |

Literacy and urbanization are highest-value: they directly relate to AIAN-specific
mechanisms of linkage failure (reservation vs. urban enumeration, record quality).
`lit_son` is already in `aian_merged` but not in the PS model.

**Fix**: Add urbanization and literacy at minimum. Any variable added to the PS model
must also be available in the comparison sample (`aian_full` from `usa_00021.csv`).
If the extract doesn't include these fields, a new IPUMS extract is needed.

---

### 2.5 No weight trimming or diagnostics (lines 101–105)

Observations with very low `p_hat` receive very large ATC weights. No trimming,
stabilization, or diagnostic output is produced. A small number of extreme weights
can dominate the weighted transition matrix.

**Fix**: Add diagnostics and trim:

```r
aian_ps = aian_ps |>
  mutate(
    w_trim    = pmin(w_atc, quantile(w_atc, 0.99)),
    w_trim_norm = w_trim * n() / sum(w_trim)
  )

# Diagnostics
summary(aian_ps$w_atc_norm)
ess_global = sum(aian_ps$w_atc_norm)^2 / sum(aian_ps$w_atc_norm^2)
```

Report both trimmed and untrimmed estimates and note where they diverge.

---

### 2.6 No covariate balance diagnostics

Weights are estimated but never checked for whether they achieve balance. Standard
practice is to report standardized mean differences (SMDs) for each covariate, weighted
and unweighted. If any SMD exceeds 0.1 after weighting, the PS model is inadequate.

**Fix**:

```r
library(cobalt)
bal.tab(linked ~ cohort * region + education * region + statefip_1940,
        data = aian_comb, weights = "w_atc_norm",
        method = "weighting", estimand = "ATC")
```

**Reporting state balance**: A 48-row state balance table is unwieldy for the paper.
Report balance on the substantive covariates (cohort, education, region, literacy,
urbanization) in the main balance table. State is a nuisance variable — it improves
prediction but the paper makes no state-level claims. For state, report a one-line
summary in the text or appendix:

```r
state_smds = bt$Balance |>
  filter(grepl("statefip", rownames(bt$Balance))) |>
  pull(Diff.Adj)

cat("State balance — max |SMD|:", max(abs(state_smds)),
    " mean |SMD|:", mean(abs(state_smds)), "\n")
```

e.g., "After weighting, the maximum SMD across state indicators is 0.XX." If any
individual state exceeds 0.1–0.2, investigate (likely a state with very few linked
observations receiving extreme weights).

This is non-negotiable for publication — reviewers cannot assess whether the weights
are doing their job without it.

---

### 2.7 No common support check

Observations with very high or very low propensity scores may fall outside the region
of common support between the linked and unlinked samples. These observations receive
extreme weights and contribute disproportionately to estimates.

**Fix**: Plot PS distributions and trim:

```r
ggplot(aian_comb, aes(x = p_hat, fill = factor(linked))) +
  geom_density(alpha = 0.5) +
  labs(x = "Propensity score", fill = "Linked")
```

Trim observations outside the overlap region, or at minimum report the proportion
of observations that fall outside common support.

---

### 2.8 Bootstrap does not account for weight estimation uncertainty

**This is the most consequential methodological issue in the paper.**

All bootstrap inference in `transition_matrices_weighted.R` resamples observations
but reuses the weights estimated on the full sample. This treats the propensity scores
as known population quantities rather than estimates, understating all standard errors
and confidence intervals. The bias is worst for regional estimates (small n, high
weight variance) and for any statistic sensitive to extreme weights.

**Fix — full weight bootstrap**:

Move the weight estimation into a function that takes a dataset and returns a weighted
analysis dataset, then call it on each bootstrap resample. In `weighting.R`:

```r
compute_weights = function(df_linked, df_full) {
  df_comb = bind_rows(
    df_linked |> mutate(linked = 1),
    df_full   |> mutate(linked = 0)
  ) |>
    mutate(
      region     = as.factor(region),
      education  = as.factor(education),
      birthyr_son = as.factor(birthyr_son)
    )

  model = glm(
    linked ~ birthyr_son * region + education * region,
    data   = df_comb,
    family = binomial
  )

  df_linked |>
    mutate(
      p_hat     = predict(model, newdata = df_linked, type = "response"),
      w_atc     = (1 - p_hat) / p_hat,
      w_atc_norm = w_atc * n() / sum(w_atc)
    )
}
```

Then in the bootstrap loops, resample `aian_merged` (the linked sample before
weighting), recompute weights, and pass the reweighted resample to the statistic:

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

**Practical note**: Re-estimating the PS model on each of R=1000 bootstrap draws is
expensive — each draw refits a logit on ~10k+ observations. R=200–500 draws is
usually sufficient for SE estimation. The full weight bootstrap should be run for
final reported estimates; the fast fixed-weight bootstrap is acceptable for
exploratory analysis.

---

## 3. Regional Analysis (`transition_matrices_weighted.R`)

### 3.1 Within-region weights are not renormalized

The global ATC weights sum to N (total linked sample size). When subsetting to a
region, `sum(w_atc_norm)` within that region is some arbitrary fraction of N rather
than n_region. Proportions and expectations computed using these weights are
technically still consistent (because region is in the PS model), but the effective
sample size is miscounted and the weighted proportions do not integrate to 1 within
the region.

**Fix**: Renormalize at the top of `compute_mobility_stats`:

```r
compute_mobility_stats = function(df) {
  df = df |> mutate(w_atc_norm = w_atc_norm / sum(w_atc_norm) * n())
  # ... rest unchanged
}
```

And simplify the regional loop using `group_modify`, which passes the subsetted
dataframe directly:

```r
results_region = data |>
  group_by(region) |>
  group_modify(~ compute_mobility_stats(.x)) |>
  ungroup()
```

---

### 3.2 Regional OM map uses unweighted estimates (line 1087)

`om_1_plot` maps `om_1_unweighted` rather than the weighted `om_1`. All other national
estimates in the paper use ATC-weighted statistics. If this is intentional (e.g., because
regional weight renormalization is unresolved per 3.1), it should be noted in the paper.
If not, it is a bug — change `om_1_unweighted` to `om_1` on lines 1087 and 1089.

---

### 3.3 Regional maps show no sample size or uncertainty

`om_1_plot`, `d_1_prime_plot`, and `upward_downward_plot` present point estimates for
all 12 regions with equal visual weight. Several regions (ne, lakes, midwest) likely
have very few linked AIAN pairs — perhaps fewer than 50 — making the 4x4 transition
matrix unreliable.

**Fix**: Compute effective sample size (ESS) per region and either grey out regions
below a threshold or overlay the raw n:

```r
ess_by_region = data |>
  group_by(region) |>
  summarise(
    n_raw = n(),
    ess   = sum(w_atc_norm)^2 / sum(w_atc_norm^2),
    .groups = "drop"
  )
```

Join `ess_by_region` into `centroids` and overlay on the maps. Consider suppressing
estimates for regions with ESS < 30 (or whatever threshold is defensible).

---

### 3.4 The upward/downward mobility ratio is not well defined

`p_upward` and `p_downward` are computed over different sub-populations of fathers,
so their denominators are not comparable. The ratio `p_upward / p_downward` cannot be
interpreted as "X upward movers per downward mover" because the populations
generating the numerator and denominator are different slices of the data.

Additionally, the downward definition uses meso categories (`meso_pop == "farmer"`,
`meso_pop == "crafts"`) while the upward definition uses macro categories — the two
measures are not constructed symmetrically.

**Fix**: Either (a) define both measures over the same base population (all fathers
except the top and bottom class, respectively), or (b) report the two rates separately
rather than as a ratio, or (c) replace with the standard "absolute upward mobility"
measure: `P(macro_son == "nonmanual")` unconditional on father's class, which has a
clean interpretation and is directly comparable across regions.

---

## 4. Modeling Assumptions

### 4.1 Stationarity across cohorts is assumed but not tested

Sons aged 20–44 in 1940 were born 1896–1920. WWI, the 1920s labour market, and the
early Depression all fall within this window and plausibly affected AIAN occupational
mobility differently across cohorts. The current P matrix pools all cohorts.

**Fix**: Split by birth cohort (e.g., <=1905 / 1906–1915 / >=1916) and compare
transition matrices. Report whether convergence measures differ materially.

---

### 4.2 Ergodicity is assumed but not verified

`pi_star` computes the stationary distribution via eigendecomposition and takes the
eigenvector associated with the eigenvalue closest to 1. If the chain is not ergodic
(e.g., a zero row in P for a sparse regional matrix), there may be multiple unit
eigenvalues and the returned pi* is arbitrary.

**Fix**: After computing P in any context where pi* is used, verify there is exactly
one unit eigenvalue:

```r
verify_ergodic = function(P) {
  eigs = abs(Re(eigen(t(P))$values))
  n_unit = sum(abs(eigs - 1) < 1e-8)
  if (n_unit != 1)
    warning(sprintf("P has %d unit eigenvalues; stationary distribution is not unique.", n_unit))
}
```

---

### 4.3 t=0 on the overall mobility plot is ambiguous

`om(P, pi0, t=0)` = 1 - sum_i pi0_i P_ii is a function of both the diagonal of P and
the initial distribution. It is not a "baseline mobility" in any natural sense — it just
says how many men in the observed initial distribution are in classes with low
persistence. Including it on the OM/EM curve implies a pre-generational baseline that
the data do not support: only t=1 is directly observed; t>=2 are extrapolations under
stationarity.

**Fix**: Either drop t=0 from the plot, or mark t=1 as "observed" and t>=2 as
"projected under stationarity."

---

## 5. Employment Classification (`cleaning-script.R`)

### 5.1 Variable overview

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

### 5.2 Current classification approach

`classify_meso(occ_son)` in `utils.R` is purely occ-first:
- `occ <= 970` → occupation category (farmer, farmworker, nonmanual, crafts, unskilled)
- `occ > 970` → nilf

`labforce_1940`, `wkswork1_1940`, and `empstatd_1940` are carried through to
`aian_merged.csv` but play no role in classification. `nilf` and `unemp` (occ 980–998)
are collapsed into a single category — intentional, as non-employment is kept as a
substantively meaningful state for this population.

**This is consistent with the occ-first conclusion below.**

---

### 5.3 Gray-area cases (aian_merged, n = 11,890)

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
check (see 5.5).

**Note on occ distribution of gray-area cases:** ~94% of gray-area cases fall in the
three bottom meso categories (farmworker, unskilled, farmer), which is consistent
with intermittent employment in casualized occupations.

---

### 5.4 Pending robustness flag

A `wkswork_flag` analogous to `spread_flag` has not yet been implemented. Candidate
definition: `empstatd_1940 == 21 & wkswork1_1940 == 0` (346 seeking-work cases with
0 prior-year weeks). This would allow re-running the main estimates excluding this
group without changing the primary classification. The 104 `not in LF + valid occ +
0 wks` cases could be included in the flag as a stricter variant.

---

### 5.5 Two-matrix proposal (occupational identity vs. labor market attachment)

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
