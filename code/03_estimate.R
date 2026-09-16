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

data      = load_global()   # sets macro_levels and meso_levels in this frame
aian_full = readRDS("data/aian_full.rds")

################################################################################
# GLOBAL BOOTSTRAP TRANSITION MATRICES
################################################################################

p_mat_macro = boot_pmatrix_ci(data, macro_pop, macro_son,
                               df_linked = data, df_full = aian_full,
                               R = 2000, .seed = 123)

p_mat_meso  = boot_pmatrix_ci(data, meso_pop, meso_son,
                               df_linked = data, df_full = aian_full,
                               R = 2000, .seed = 123)

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

# Pairwise TV distances — full matrix at macro; sanity-check max == exp(dobrushin$d1)
rd_macro = row_dists(P_macro)
cat("\nGlobal macro pairwise TV distances:\n")
print(round(rd_macro, 3))
stopifnot(abs(max(rd_macro) - exp(dobrushin(P_macro)$d1)) < 1e-10)

################################################################################
# GLOBAL MOBILITY SCALARS
################################################################################

om_macro = om(P_macro, pi0_macro, t = 0)
sm_macro = sm(P_macro, pi0_macro, t = 0)

om_meso  = om(P_meso, pi0_meso, t = 0)
sm_meso  = sm(P_meso, pi0_meso, t = 0)

cat(sprintf("\nGlobal macro  OM=%.3f  SM=%.3f\n", om_macro, sm_macro))
cat(sprintf("Global meso   OM=%.3f  SM=%.3f\n",  om_meso,  sm_meso))

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

region_display = c(
  sw      = "Southwest",
  nplains = "Northern Plains",
  glakes  = "Great Lakes",
  nw      = "Northwest",
  ok      = "Oklahoma",
  cali    = "California",
  nc      = "North Carolina",
  basin   = "Basin and Mountain"
)

compute_regional = function(df_reg) {
  df_reg = renorm(df_reg)
  P      = p_matrix(df_reg, macro_pop, macro_son, matrix = TRUE)
  pis    = pi_star(P)
  pi0    = pi_0(df_reg, macro_pop)
  dob    = dobrushin(P)

  eig_mods = sort(Mod(eigen(P, only.values = TRUE)$values), decreasing = TRUE)
  lambda2  = eig_mods[2]

  relief_lower = weighted.mean(df_reg$empstatd_1940 == 11,  df_reg$w_atc_norm, na.rm = TRUE)
  relief_upper = weighted.mean(df_reg$classwkrd_1940 == 24, df_reg$w_atc_norm, na.rm = TRUE)

  # Farmworker share: proportion of farming category that is farmworker, by generation
  farm_sons = dplyr::filter(df_reg, macro_son == "farming")
  farmwkr_share_son = if (nrow(farm_sons) > 0)
    weighted.mean(farm_sons$meso_son == "farmworker", farm_sons$w_atc_norm, na.rm = TRUE)
  else NA_real_

  farm_dads = dplyr::filter(df_reg, macro_pop == "farming")
  farmwkr_share_pop = if (nrow(farm_dads) > 0)
    weighted.mean(farm_dads$meso_pop == "farmworker", farm_dads$w_atc_norm, na.rm = TRUE)
  else NA_real_

  cell  = function(r, c) if (r %in% rownames(P) && c %in% colnames(P)) P[r, c] else NA_real_
  pisel = function(k)    if (k %in% names(pis)) pis[k] else NA_real_

  list(
    P                 = P,
    pi0               = pi0,
    pistar            = pis,
    n                 = nrow(df_reg),
    lambda2           = lambda2,
    d1                = dob$d1,
    d1_row1           = dob$row1,
    d1_row2           = dob$row2,
    om1               = om(P, pi0, t = 0),
    sm1               = sm(P, pi0, t = 0),
    farm_ret          = cell("farming", "farming"),
    pi0_farming       = if ("farming" %in% names(pi0)) pi0["farming"] else NA_real_,
    farmwkr_share_son = farmwkr_share_son,
    farmwkr_share_pop = farmwkr_share_pop,
    relief_lower      = relief_lower,
    relief_upper      = relief_upper,
    pi_farming        = pisel("farming"),
    pi_manual         = pisel("manual"),
    pi_nonman         = pisel("nonmanual"),
    pi_nonemp         = pisel("nonemp")
  )
}

regional_data = load_regional()

regional_results = setNames(
  lapply(compare_regions, function(r) compute_regional(regional_data[[r]])),
  compare_regions
)

for (r in compare_regions) {
  res = regional_results[[r]]
  cat(sprintf(
    "  %s  n=%d  farm_ret=%.3f  pi0_farm=%.3f  d1=%.3f [%s vs %s]  relief_lo=%.3f  relief_hi=%.3f\n",
    region_display[r], res$n, res$farm_ret, res$pi0_farming,
    res$d1, res$d1_row1, res$d1_row2, res$relief_lower, res$relief_upper))
}

################################################################################
# REGIONAL COUNTERFACTUALS — leave-one-out benchmark
# Benchmark for region k = pool of all other compare_regions.
# pop_n scales each region's weights to its target-population count (option a).
################################################################################

pop_n = table(aian_full$region)

reg_cf = regional_counterfactuals(regional_data, macro_pop, macro_son,
                                   compare = compare_regions, pop_n = pop_n)

cat("\n--- Regional counterfactuals (LOO benchmark) ---\n")
print(reg_cf)

################################################################################
# REGIONAL BOOTSTRAP — benchmark rebuilt in every draw
################################################################################

reg_cf_boot = boot_regional_cf(regional_data, macro_pop, macro_son,
                                compare = compare_regions, pop_n = pop_n,
                                R = 2000, .seed = 456)

cat("\n--- Regional counterfactual bootstrap intervals ---\n")
print(reg_cf_boot)

################################################################################
# SAVE
################################################################################

estimates = list(
  # Global bootstrap (list: P, pi_s, d1)
  p_mat_macro   = p_mat_macro,
  p_mat_meso    = p_mat_meso,
  # Global point-estimate matrices
  P_macro       = P_macro,
  P_meso        = P_meso,
  # Pairwise TV distances
  rd_macro      = rd_macro,
  # Distributions
  pi0_macro     = pi0_macro,
  pi0_meso      = pi0_meso,
  steady_macro  = steady_macro,
  steady_meso   = steady_meso,
  # Dobrushin
  dob_mac       = dob_mac,
  dob_mes       = dob_mes,
  # Scalar mobility (global, t=0: observed father→son generation)
  om_macro      = om_macro,
  sm_macro      = sm_macro,
  om_meso       = om_meso,
  sm_meso       = sm_meso,
  # Regional (named list) and LOO counterfactuals
  regional        = regional_results,
  reg_cf          = reg_cf,
  reg_cf_boot     = reg_cf_boot,
  compare_regions = compare_regions,
  region_display  = region_display
)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
saveRDS(estimates, "output/estimates.rds")
cat("\nWrote output/estimates.rds\n")
