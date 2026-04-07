# Methods Notes: Outstanding Issues and Solutions

This document records methodological concerns identified during code review, along with
proposed fixes where they exist. Organized by script/stage of the pipeline.

---

## 1. Cleaning (`cleaning-script.R`)

### 1.1 Father birth year is averaged, but the spread is computed and never used (lines 68–73)

```r
summarise(across(starts_with("birthyr"), ~ mean(.x, na.rm = TRUE)))
...
spread = diff(range(...))
```

`birth_median` is taken as the median of birth year estimates across census appearances.
`spread` (the range of those estimates) is computed but discarded. Fathers with large
`spread` (e.g., > 5 years) have inconsistent age reports across censuses and are
plausible false matches. These cases should be flagged or excluded.

**Fix**: Join `spread` back into `aian_merged` and add a robustness check excluding
pairs with `spread > 5`.

---

### 1.2 Multi-valued numeric fields are silently coerced to NA (lines 50–54)

```r
across(24:55, ~ paste(unique(na.omit(.x)), collapse = "; "))
...
across(25:56, as.integer)
```

When a son is linked to two records in the same census wave, numeric fields (including
`occ1950_pop_*`) are collapsed into semicolon-separated strings. `as.integer("100; 200")`
returns `NA` silently. Fathers with multiple occupation records in the same year are
dropped entirely rather than handled explicitly.

**Fix**: Before `as.integer`, take the first value of multi-valued numeric fields:

```r
across(25:56, ~ as.integer(sub(";.*", "", .x)))
```
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
