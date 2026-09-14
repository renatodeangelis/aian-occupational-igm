################################################################################
# 03_estimate.R
# Estimate transition matrices, distributions, and mobility scalars.
# ONLY script that calls p_matrix(), pi_0(), pi_star(), boot_pmatrix_ci().
#
# Reads:  data/aian_weighted.rds          (via load_global())
#         data/aian_regional_weighted.rds  (via load_regional())
#         data/aian_full.rds
# Writes: output/estimates.rds
################################################################################

library(dplyr)
library(expm)

source("code/00_utils.R")
source("code/expected_values.R")

data      = load_global()   # sets macro_levels and meso_levels in this frame
aian_full = readRDS("data/aian_full.rds")

################################################################################
# GLOBAL BOOTSTRAP TRANSITION MATRICES
################################################################################

p_mat_macro = boot_pmatrix_ci(data, macro_pop, macro_son,
                               df_linked = data, df_full = aian_full,
                               R = 500, .seed = 123)

p_mat_meso  = boot_pmatrix_ci(data, meso_pop, meso_son,
                               df_linked = data, df_full = aian_full,
                               R = 500, .seed = 123)

################################################################################
# GLOBAL POINT ESTIMATES
################################################################################

P_macro = p_matrix(data, macro_pop, macro_son, matrix = TRUE)
P_meso  = p_matrix(data, meso_pop,  meso_son,  matrix = TRUE)

verify_ergodic(P_macro, "macro global")
verify_ergodic(P_meso,  "meso global")

pi0_macro    = pi_0(data, macro_pop)
pi0_meso     = pi_0(data, meso_pop)
steady_macro = pi_star(P_macro)
steady_meso  = pi_star(P_meso)

cat("\nGlobal macro pi_0:\n");   print(round(pi0_macro,    3))
cat("Global macro pi*:\n");     print(round(steady_macro, 3))
cat("Global macro P:\n");       print(round(P_macro,      3))
cat("\nGlobal meso pi_0:\n");   print(round(pi0_meso,     3))
cat("Global meso pi*:\n");      print(round(steady_meso,  3))

################################################################################
# GLOBAL MOBILITY SCALARS
################################################################################

om_macro = om(P_macro, pi0_macro, t = 0)
sm_macro = sm(P_macro, pi0_macro, t = 0)
em_macro = om_macro - sm_macro

om_meso  = om(P_meso, pi0_meso, t = 0)
sm_meso  = sm(P_meso, pi0_meso, t = 0)
em_meso  = om_meso - sm_meso

cat(sprintf(
  "\nGlobal macro  OM=%.3f  SM=%.3f  EM=%.3f\n",
  om_macro, sm_macro, em_macro))
cat(sprintf(
  "Global meso   OM=%.3f  SM=%.3f  EM=%.3f\n",
  om_meso, sm_meso, em_meso))

################################################################################
# DOBRUSHIN CONTRACTION
################################################################################

dob_mac = dobrushin(P_macro)
dob_mes = dobrushin(P_meso)

cat(sprintf(
  "\nDobrushin d1:  macro=%.4f [%s vs %s]   meso=%.4f [%s vs %s]\n",
  dob_mac$d1, dob_mac$row1, dob_mac$row2,
  dob_mes$d1, dob_mes$row1, dob_mes$row2))

################################################################################
# REGIONAL POINT ESTIMATES
# Uses region-specific PS model data (load_regional()); w_atc_norm already
# trimmed. macro_levels / meso_levels set by load_global() above satisfy TRAP-1.
################################################################################

regions_list = c("sw", "south", "cali", "ok", "plains", "nw", "north")

region_display = c(
  sw     = "Southwest",
  south  = "South",
  cali   = "California",
  ok     = "Oklahoma",
  plains = "Plains",
  nw     = "Northwest",
  north  = "North"
)

compute_regional = function(df_reg) {
  df_reg = renorm(df_reg)
  P      = p_matrix(df_reg, macro_pop, macro_son, matrix = TRUE)
  pis    = pi_star(P)
  pi0    = pi_0(df_reg, macro_pop)

  eig_mods = sort(Mod(eigen(P, only.values = TRUE)$values), decreasing = TRUE)
  lambda2  = eig_mods[2]

  relief = weighted.mean(df_reg$empstatd_1940 == 11, df_reg$w_atc_norm,
                         na.rm = TRUE)

  cell   = function(r, c) if (r %in% rownames(P) && c %in% colnames(P)) P[r, c] else NA_real_
  pisel  = function(k)    if (k %in% names(pis)) pis[k] else NA_real_

  list(
    P          = P,
    pi0        = pi0,
    pistar     = pis,
    n          = nrow(df_reg),
    lambda2    = lambda2,
    om1        = om(P, pi0, t = 0),
    sm1        = sm(P, pi0, t = 0),
    farm_ret   = cell("farming", "farming"),
    pi_farming = pisel("farming"),
    pi_manual  = pisel("manual"),
    pi_nonman  = pisel("nonmanual"),
    pi_nonemp  = pisel("nonemp"),
    relief     = relief
  )
}

regional_data = load_regional()

regional_results = setNames(
  lapply(regions_list, function(r) compute_regional(regional_data[[r]])),
  regions_list
)

for (r in regions_list) {
  res = regional_results[[r]]
  cat(sprintf("  %s  n=%d  farm_ret=%.3f  pi*_farm=%.3f  lambda2=%.3f  relief=%.3f\n",
              region_display[r], res$n, res$farm_ret,
              res$pi_farming, res$lambda2, res$relief))
}

################################################################################
# ASSERTIONS
################################################################################

cat("\n--- Checking global pi_0 ---\n")
pi0_check = pi0_macro[macro_compute_order]
exp_pi0   = EXPECTED$pi0_global_macro
if (!all(abs(pi0_check - exp_pi0) < TOL)) {
  stop(sprintf(
    "Global pi_0 mismatch.\n  got:      %s\n  expected: %s",
    paste(round(pi0_check, 3), collapse = " / "),
    paste(exp_pi0,             collapse = " / ")))
}
cat("Global pi_0 check passed:", paste(round(pi0_check, 3), collapse = " / "), "\n")

cat("\n--- Checking regional estimates ---\n")
num_fields = c("farm_ret", "pi_farming", "pi_manual", "pi_nonman", "lambda2", "relief")

for (r in regions_list) {
  exp = EXPECTED$regional[[r]]
  got = regional_results[[r]]
  if (is.null(exp)) next

  for (fld in num_fields) {
    if (is.na(exp[[fld]])) next
    delta = abs(got[[fld]] - exp[[fld]])
    if (delta > TOL)
      stop(sprintf("[%s] %s: got %.4f, expected %.4f (|delta|=%.4f > TOL=%.3f)",
                   r, fld, got[[fld]], exp[[fld]], delta, TOL))
  }
  if (!is.na(exp$n) && got$n != exp$n)
    stop(sprintf("[%s] n: got %d, expected %d", r, got$n, exp$n))
}
cat("All regional assertions passed.\n")

################################################################################
# SAVE
################################################################################

estimates = list(
  # Global bootstrap tibbles
  p_mat_macro   = p_mat_macro,
  p_mat_meso    = p_mat_meso,
  # Global point-estimate matrices
  P_macro       = P_macro,
  P_meso        = P_meso,
  # Distributions
  pi0_macro     = pi0_macro,
  pi0_meso      = pi0_meso,
  steady_macro  = steady_macro,
  steady_meso   = steady_meso,
  # Dobrushin
  dob_mac       = dob_mac,
  dob_mes       = dob_mes,
  # Scalar mobility (global, t=0)
  om_macro      = om_macro,
  sm_macro      = sm_macro,
  em_macro      = em_macro,
  om_meso       = om_meso,
  sm_meso       = sm_meso,
  em_meso       = em_meso,
  # Regional (named list)
  regional      = regional_results,
  regions_list  = regions_list,
  region_display = region_display
)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
saveRDS(estimates, "output/estimates.rds")
cat("\nWrote output/estimates.rds\n")
