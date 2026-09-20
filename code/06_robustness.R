################################################################################
# 06_robustness.R
# Robustness checks against the main macro transition matrix estimates.
#
# §1  Cohort stationarity (birth cohorts 1896-1905, 1906-1915, 1916-1920)
# §2  Age window: 25-44 vs main 20-44
# §3  Employed-only: restrict to occ_son <= 970 & occ_pop <= 970
# §4  Unweighted vs weighted transition matrix comparison
#
# Reads:  data/aian_weighted.rds  (via load_global)
#         data/aian_merged.rds    (for occ_pop in §3)
# All output printed to console.
################################################################################

library(dplyr)
library(purrr)
library(expm)

source("code/00_utils.R")

data = load_global()   # sets macro_levels / meso_levels

################################################################################
# §1 COHORT STATIONARITY
# Split into three birth cohorts; compare macro P and scalar mobility.
# Tests whether the pooled P is stable across the observation window.
################################################################################

cat("================================================================================\n")
cat("§1 COHORT STATIONARITY\n")
cat("================================================================================\n")

cohort_labels = c("1891-1895", "1896-1905", "1906-1915", "1916-1920")

data_cohorts = data |>
  mutate(cohort_group = cut(birthyr_son,
                            breaks = c(1890, 1895, 1905, 1915, 1921),
                            labels = cohort_labels))
stopifnot(!any(is.na(data_cohorts$cohort_group)))

cohort_ns = data_cohorts |>
  count(cohort_group) |>
  mutate(ess = map_dbl(cohort_group, function(cg) {
    d = filter(data_cohorts, cohort_group == cg)
    sum(d$w_atc_norm)^2 / sum(d$w_atc_norm^2)
  }))

cat("\n--- Cohort sample sizes ---\n")
print(cohort_ns)

cohort_results = map(cohort_labels, function(cg) {
  d   = filter(data_cohorts, cohort_group == cg) |> renorm()
  P   = p_matrix(d, macro_pop, macro_son, matrix = TRUE)
  pi0 = pi_0(d, macro_pop)
  list(
    cohort = cg,
    P      = P,
    om0    = round(om(P, pi0, t = 0), 3),
    sm0    = round(sm(P, pi0, t = 0), 3),
    d1     = round(dobrushin(P)$d1,   3)
  )
})

cohort_tbl = map_dfr(cohort_results, function(res) {
  tibble(cohort = res$cohort, om_0 = res$om0, sm_0 = res$sm0, log_d1 = res$d1)
})

cat("\n--- Cohort stationarity: OM(0), SM(0), log delta(P) by birth cohort ---\n")
print(cohort_tbl)

cat("\n--- Macro transition matrices by cohort ---\n")
for (res in cohort_results) {
  cat(sprintf("\nCohort: %s\n", res$cohort))
  print(round(res$P, 3))
}

################################################################################
# §2 AGE WINDOW: 25-44 VS MAIN 20-44
################################################################################

cat("\n================================================================================\n")
cat("§2 AGE WINDOW (25-44 vs 20-49)\n")
cat("================================================================================\n")

data_2544 = data |>
  filter((1940 - birthyr_son) >= 25, (1940 - birthyr_son) <= 44) |>
  renorm()

cat(sprintf("\nMain (20-44): n = %d   Restricted (25-44): n = %d   Loss = %.1f%%\n",
            nrow(data), nrow(data_2544),
            (1 - nrow(data_2544) / nrow(data)) * 100))

P_main   = p_matrix(data,     macro_pop, macro_son, matrix = TRUE)
P_2544   = p_matrix(data_2544, macro_pop, macro_son, matrix = TRUE)
pi0_main = pi_0(data,      macro_pop)
pi0_2544 = pi_0(data_2544, macro_pop)

robustness_tbl = tibble(
  sample     = c("main (20-49)", "restricted (25-44)"),
  n          = c(nrow(data), nrow(data_2544)),
  om_0       = round(c(om(P_main, pi0_main, t=0), om(P_2544, pi0_2544, t=0)), 3),
  sm_0       = round(c(sm(P_main, pi0_main, t=0), sm(P_2544, pi0_2544, t=0)), 3),
  log_d1     = round(c(dobrushin(P_main)$d1,       dobrushin(P_2544)$d1),       3),
  farm_farm  = round(c(P_main["farming","farming"], P_2544["farming","farming"]), 3)
)

cat("\n--- Age restriction robustness ---\n")
print(robustness_tbl)

cat("\nRestricted (25-44) macro P:\n"); print(round(P_2544, 3))
cat("\nMain (20-44) macro P:\n");        print(round(P_main,  3))

################################################################################
# §3 EMPLOYED-ONLY
# Restrict to observations where both father and son are employed (OCC <= 970).
# Uses occ_pop from aian_merged.rds (passed through 01_cleaning.R → 02_weighting.R).
################################################################################

cat("\n================================================================================\n")
cat("§3 EMPLOYED-ONLY (occ_son <= 970 & occ_pop <= 970)\n")
cat("================================================================================\n")

if (!"occ_pop" %in% names(data)) {
  cat("NOTE: occ_pop not present in aian_weighted.rds — re-run 01_cleaning.R and 02_weighting.R.\n")
  cat("Skipping §3.\n")
} else {
  data_emp = data |>
    filter(occ_son <= 970, occ_pop <= 970) |>
    renorm()

  cat(sprintf("\nFull: n = %d   Employed-only: n = %d   Loss = %.1f%%\n",
              nrow(data), nrow(data_emp),
              (1 - nrow(data_emp) / nrow(data)) * 100))

  P_emp   = p_matrix(data_emp, macro_pop, macro_son, matrix = TRUE)
  pi0_emp = pi_0(data_emp, macro_pop)

  emp_tbl = tibble(
    sample    = c("main (all)", "employed-only"),
    n         = c(nrow(data), nrow(data_emp)),
    om_0      = round(c(om(P_main, pi0_main, t=0), om(P_emp, pi0_emp, t=0)), 3),
    sm_0      = round(c(sm(P_main, pi0_main, t=0), sm(P_emp, pi0_emp, t=0)), 3),
    log_d1    = round(c(dobrushin(P_main)$d1,       dobrushin(P_emp)$d1),     3),
    farm_farm = round(c(P_main["farming","farming"], P_emp["farming","farming"]), 3)
  )

  cat("\n--- Employed-only robustness ---\n")
  print(emp_tbl)

  cat("\nEmployed-only macro P:\n"); print(round(P_emp, 3))
}

################################################################################
# §4 UNWEIGHTED VS WEIGHTED
################################################################################

cat("\n================================================================================\n")
cat("§4 UNWEIGHTED VS WEIGHTED\n")
cat("================================================================================\n")

P_unw   = p_matrix_unweighted(data, macro_pop, macro_son, matrix = TRUE)
pi0_unw = pi_0_unweighted(data, macro_pop)

cat("\nUnweighted macro P:\n"); print(round(P_unw, 3))
cat("Weighted macro P:\n");    print(round(P_main, 3))

cat("\nMax absolute cell difference:", round(max(abs(P_unw - P_main[rownames(P_unw), colnames(P_unw)])), 3), "\n")

unw_tbl = tibble(
  sample    = c("weighted", "unweighted"),
  om_0      = round(c(om(P_main, pi0_main, t=0), om(P_unw, pi0_unw, t=0)), 3),
  sm_0      = round(c(sm(P_main, pi0_main, t=0), sm(P_unw, pi0_unw, t=0)), 3),
  log_d1    = round(c(dobrushin(P_main)$d1,       dobrushin(P_unw)$d1),     3),
  farm_farm = round(c(P_main["farming","farming"], P_unw["farming","farming"]), 3)
)

cat("\n--- Weighted vs unweighted mobility scalars ---\n")
print(unw_tbl)
