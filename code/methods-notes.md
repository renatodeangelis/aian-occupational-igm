# Methods Notes: Outstanding Issues and Solutions

This document records methodological concerns identified during code review, along with
proposed fixes where they exist. Organized by script/stage of the pipeline.

---

## 1. Cleaning (`cleaning-script.R`)

### 1.1 Minimum father–son age gap is too loose (line 81)

```r
filter(birthyr_son > birthyr_pop + 15)
```

A 15-year minimum allows fathers who were 16 at the birth of their son. The IGM
literature standard is 25–30 years. At short gaps, the "father" may be an older
brother or other male relative misidentified in the linking — a particular concern for
AIAN households, which were often multi-generational and enumerated inconsistently.

**Fix**: Test sensitivity at 20- and 25-year minimums and report whether estimates change.

---

### 1.2 Father birth year is averaged, but the spread is computed and never used (lines 68–73)

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

### 1.3 Multi-valued numeric fields are silently coerced to NA (lines 50–54)

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

### 1.4 Occupation code 0 (NIU) is classified as "prof" (line 121)

```r
prof_codes = c(0:99, 200:290)
```

In the occ1950 scheme, code 0 is "not in universe" — assigned to people outside the
labour force. These are misclassified as professional workers. This affects sons whose
1940 occupation was not recorded and any fathers with `occ = 0`.

**Fix**:

```r
prof_codes = c(1:99, 200:290)
```

---

### 1.5 "unemp" conflates unemployment, non-employment, and missing (line 131)

```r
occ > 970 ~ "unemp"
```

occ1950 codes 971–999 include: housewife, student, retired/unable to work, inmate,
and 999 (occupation not reported). These are categorically different from unemployed
and seeking work. Treating them as equivalent inflates the "unemp" category and
distorts the transition matrix for men who are simply not in the labour force.

**Fix**: Rename the category to `"nilf"` (not in labour force) and describe its
composition in the data section. Consider splitting 999 (missing) from 971–998
(explicit non-participation) and reporting how many observations fall in each.

---

### 1.6 "prof" collapses professionals with managers and proprietors (line 121)

`prof_codes = c(1:99, 200:290)` (after fix 1.4) groups occ1950's
"professional/technical" (physicians, lawyers, engineers: 1–99) with "managers,
officials, and proprietors" (200–290). These groups differ substantially in education
requirements, intergenerational transmission, and relevance to AIAN economic
history. Small traders and local proprietors (200–290) are not comparable to
licensed professionals.

Now partially mitigated by the `white_collar` collapse at the meso level, but the
macro-level "nonmanual" category still combines them.

**Fix**: Justify in the paper, or split into "prof" (1–99) and "manager" (200–290)
and test whether their transition profiles differ.

---

### 1.7 Crafts codes extracted from 595–970 are undocumented (line 122)

```r
crafts_codes = c(762, 773, 781, 782)
```

Four specific construction trades (painters, plasterers, plumbers, structural metal
workers) are rescued from the "unskilled" range into "crafts." The reason only these
four codes need this treatment — and not other potentially skilled trades in 595–970 —
is unexplained. The selection is arbitrary-seeming without documentation.

**Fix**: Document the source of these codes and check whether other skilled trades in
595–970 (e.g., bricklayers, electricians if miscoded) should also be included.

---

### 1.8 `classify_meso` and `classify_macro` are duplicated (lines 119–141 and `weighting.R:29–51`)

Identical functions appear in both `cleaning-script.R` and `weighting.R`. These have
already silently diverged: `classify_macro` in `weighting.R` still maps
`c("prof", "clerical") ~ "nonmanual"` while `transition_matrices_weighted.R` has
collapsed these into `white_collar` as a post-load mutation. Any future change to the
classification must be made in two places, and the two scripts can silently produce
inconsistent results.

**Fix**: Move both functions to a shared `utils.R` and `source()` it from both scripts.

---

### 1.9 Father age filter at 65 excludes fathers with no estimable birth year (line 110)

```r
filter(dplyr::coalesce(implied_age, Inf) <= 65)
```

`coalesce(NA, Inf)` sets missing implied age to Inf, so fathers whose birth year
cannot be estimated are excluded from the modal occupation selection entirely. This is
an implicit restriction that is not acknowledged. These fathers may be systematically
different (e.g., linked through a single census appearance with no age recorded).

**Fix**: Document explicitly. Consider whether fathers with no estimable birth year
should be retained with occupation drawn from whatever census years are available,
rather than dropped.

---

## 2. Weighting (`weighting.R`)

### 2.1 Missing lower age bound in the full AIAN comparison sample (line 56)

```r
filter(age < 45, school == 1)
```

The unlinked AIAN census extract has no lower age bound. If `school == 1` means
"not attending school" (standard IPUMS 1940 coding), children aged 4+ who are not
enrolled satisfy this filter. These would be classified as "prof" (via the occ=0 bug
above) and included in the PS estimation as controls, severely distorting the model.

**Fix**:

```r
filter(between(age, 20, 44), school == 1)
```

---

### 2.2 PS model does not include region × covariate interactions (line 95)

```r
ps_model = glm(linked ~ birthyr_son + meso_son + region + education,
               data = aian_comb, family = binomial)
```

The model includes region as a main effect, assuming the relationship between linkage
and (birth year, education, occupation) is the same in all regions. This is unlikely —
BIA record-keeping practices, reservation settlement patterns, and name stability
varied substantially by region, meaning education and birth year predict linkage
differently in the South vs. the Basin states.

**Fix**: Add region interactions:

```r
ps_model = glm(
  linked ~ birthyr_son * region + meso_son * region + education * region,
  data   = aian_comb,
  family = binomial
)
```

This is equivalent to region-stratified models but shares data across regions for
stability, avoiding separation problems in small regions. Everything downstream
(`p_hat`, ATC weights, normalization) is unchanged.

---

### 2.3 `meso_son` (an outcome) is a PS covariate (line 95)

The son's occupation is both a covariate in the PS model and the outcome of the
mobility process being studied. This creates circularity: the weights depend on where
the son ended up, which means the weighted transition probabilities are not purely
driven by father-side selection into the linked sample.

**Fix**: Re-estimate without `meso_son` and test whether national and regional
estimates change materially. If they do not, the simpler model is preferable. If they
do, discuss the tradeoff explicitly in the paper.

---

### 2.4 No weight trimming or diagnostics (lines 101–105)

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

### 2.5 Bootstrap does not account for weight estimation uncertainty

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
      meso_son   = as.factor(meso_son),
      region     = as.factor(region),
      education  = as.factor(education),
      birthyr_son = as.factor(birthyr_son)
    )

  model = glm(
    linked ~ birthyr_son * region + meso_son * region + education * region,
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

### 3.2 Regional maps show no sample size or uncertainty

`om_1_plot`, `d_1_prime_plot`, and `upward_downward_plot` present point estimates for
all 13 regions with equal visual weight. Several regions (Northeast, Great Lakes,
Midwest) likely have very few linked AIAN pairs — perhaps fewer than 50 — making the
4×4 transition matrix unreliable.

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

### 3.3 The upward/downward mobility ratio is not well defined

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

### 4.1 First-order Markov assumption is not tested

The entire analysis rests on P(son's class | father's class) being independent of
grandfather's class. For AIAN families with strong dynastic persistence (or strong
dynastic disruption from allotment/assimilation policy), this may not hold. The
linked data may contain grandfather occupations for a subset of pairs (via earlier
census linkages), which would allow a direct test.

**Suggestion**: For the subsample with grandfather information, test whether adding
grandfather's class to the transition model materially changes coefficient estimates.

---

### 4.2 Stationarity across cohorts is assumed but not tested

Sons aged 20–44 in 1940 were born 1896–1920. WWI, the 1920s labour market, and the
early Depression all fall within this window and plausibly affected AIAN occupational
mobility differently across cohorts. The current P matrix pools all cohorts.

**Fix**: Split by birth cohort (e.g., ≤1905 / 1906–1915 / ≥1916) and compare
transition matrices. Report whether convergence measures differ materially.

---

### 4.3 Ergodicity is assumed but not verified

`pi_star` computes the stationary distribution via eigendecomposition and takes the
eigenvector associated with the eigenvalue closest to 1. If the chain is not ergodic
(e.g., a zero row in P for a sparse regional matrix), there may be multiple unit
eigenvalues and the returned π* is arbitrary.

**Fix**: After computing P in any context where π* is used, verify there is exactly
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

### 4.4 t=0 on the overall mobility plot is ambiguous

`om(P, pi0, t=0)` = 1 − Σᵢ π₀ᵢ Pᵢᵢ is a function of both the diagonal of P and the
initial distribution. It is not a "baseline mobility" in any natural sense — it just
says how many men in the observed initial distribution are in classes with low
persistence. Including it on the OM/EM curve implies a pre-generational baseline that
the data do not support: only t=1 is directly observed; t≥2 are extrapolations under
stationarity.

**Fix**: Either drop t=0 from the plot, or mark t=1 as "observed" and t≥2 as
"projected under stationarity."

---

## 5. Gender and Scope

The analysis is restricted to men throughout. This is standard in the IGM literature
(due to low and selective female labour force participation in the early 20th century),
but for AIAN communities the gender dynamics of occupational inheritance may differ
substantially from white or Black patterns. The paper's claims about AIAN
intergenerational mobility are claims about male mobility specifically, and this
boundary condition should be stated explicitly in the introduction.

---

## 6. What Cannot Be Fixed Without More Data

The following limitations are inherent to the data-generating process and cannot be
addressed through reanalysis:

- **Non-random linkage on unobservables.** The PS correction handles observable
  selection. Men who are harder to link due to migration, name changes, or reservation
  residence patterns — and whose mobility trajectories may differ systematically — are
  irretrievably excluded. This is a known limitation of all CLP-based studies.

- **Single observed transition.** P is estimated from one father→son generation. All
  multi-generational projections (t≥2) are extrapolations under the stationarity
  assumption, not observations. A longer panel (grandfather→father→son) would allow
  both testing the Markov assumption and estimating multi-generational persistence
  directly.

- **Female mobility.** Occupational coding for women in early 20th-century censuses is
  not comparably rich or reliable. Female AIAN mobility is simply not estimable from
  this source.

- **Geographic resolution.** Sample sizes preclude county- or reservation-level
  analysis. The 13-region breakdown already pushes against the limits of reliable
  estimation in several regions.
