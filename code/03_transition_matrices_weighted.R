################################################################################
############################# 1. SETUP #########################################
################################################################################

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(expm)
library(purrr)
library(maps)
library(sf)
library(weights)
library(jtools)

source("code/utils.R")

data = readRDS("data/aian_weighted.rds") |>
  mutate(w_atc_norm = w_trim_norm)

aian_full = readRDS("data/aian_full.rds")

macro_levels = unique(data$macro_pop)
meso_levels = unique(data$meso_pop)

############################# 2. COMPUTATION ###################################
################################################################################

# Bootstrap cache — set force_rerun = TRUE to discard cache and re-estimate.
# Delete individual files from cache/ to selectively re-run one bootstrap.
# Note: cache is not auto-invalidated if aian_weighted.rds changes; delete
# cache/ manually after re-running weighting.R.
force_rerun = FALSE
cache_dir   = "cache"
dir.create(cache_dir, showWarnings = FALSE)

# Transition matrices with bootstrap CIs
p_mat_macro = cache_load("p_mat_macro", quote(
  boot_pmatrix_ci(data, macro_pop, macro_son,
                  df_linked = data, df_full = aian_full,
                  R = 500, .seed = 123)
), force = force_rerun)
p_mat_meso  = cache_load("p_mat_meso", quote(
  boot_pmatrix_ci(data, meso_pop, meso_son,
                  df_linked = data, df_full = aian_full,
                  R = 500, .seed = 123)
), force = force_rerun)

# Initial and stationary distributions
pi0_vec_macro = pi_0(data, macro_pop)
pi0_vec_meso  = pi_0(data, meso_pop)
P_macro_global = p_matrix(data, macro_pop, macro_son, TRUE)
P_meso_global  = p_matrix(data, meso_pop, meso_son, TRUE)
verify_ergodic(P_macro_global, "macro global")
verify_ergodic(P_meso_global,  "meso global")
steady_macro  = pi_star(P_macro_global)
steady_meso   = pi_star(P_meso_global)

# EM/SM bootstrap
macro_om = cache_load("macro_om", quote(
  mobility_curve_with_boot(data, macro_pop, macro_son,
                            df_linked = data, df_full = aian_full,
                            ts = 1:4, R = 500, .seed = 123)
), force = force_rerun)
meso_om  = cache_load("meso_om", quote(
  mobility_curve_with_boot(data, meso_pop, meso_son,
                            df_linked = data, df_full = aian_full,
                            ts = 1:4, R = 500, .seed = 123)
), force = force_rerun)

# lo/hi are approximate 95% bands from SE; OM recoverable as EM + SM at caller level.
om_total = bind_rows(macro_om |> mutate(level = 1), meso_om |> mutate(level = 0)) |>
  mutate(
    level = factor(level, labels = c("meso", "macro")),
    lo    = est - 1.96 * se,
    hi    = est + 1.96 * se
  )

################################################################################
############################# 5. PLOTTING ######################################
################################################################################

meso_level_order = meso_order

## Transition matrix heatmaps ----

g_macro      = plot_pmat_heatmap(p_mat_macro, macro_pop, macro_son, title_expr = expression(P))
g0_macro     = plot_pi_column(pi0_vec_macro, title_expr = expression(pi[0]))
g_star_macro = plot_pi_column(steady_macro, title_expr = expression(pi^"*"))
combined_plot_macro = g_macro + g0_macro + g_star_macro + plot_layout(widths = c(6, 1, 1))

g_meso      = plot_pmat_heatmap(p_mat_meso, meso_pop, meso_son,
                                levels = meso_level_order, text_size = 3, title_expr = "P")
g0_meso     = plot_pi_column(pi0_vec_meso, expression(pi[0]), levels = meso_level_order)
g_star_meso = plot_pi_column(steady_meso, expression(pi^"*"), levels = meso_level_order)
combined_plot_meso = g_meso + g0_meso + g_star_meso + plot_layout(widths = c(6, 1, 1))

################################################################################
#################### 7. SUPPLEMENTARY: COHORT STATIONARITY ####################
################################################################################

# Split into three birth cohorts and compare macro transition matrices.
# Tests whether the pooled P matrix is stable across the observation window.
# Sons born 1896-1920 span WWI, the 1920s labour market, and the early Depression.

cohort_labels = c("1896-1905", "1906-1915", "1916-1920")

data_cohorts = data |>
  mutate(cohort_group = cut(birthyr_son,
                            breaks = c(1895, 1905, 1915, 1921),
                            labels = cohort_labels))

cohort_ns = data_cohorts |>
  count(cohort_group) |>
  mutate(ess = map_dbl(cohort_group, function(cg) {
    d = filter(data_cohorts, cohort_group == cg)
    sum(d$w_atc_norm)^2 / sum(d$w_atc_norm^2)
  }))

cat("\n--- Cohort sample sizes ---\n")
print(cohort_ns)

cohort_results = cohort_labels |>
  purrr::map(function(cg) {
    d = filter(data_cohorts, cohort_group == cg) |>
      mutate(w_atc_norm = w_atc_norm / sum(w_atc_norm) * n())
    P   = p_matrix(d, macro_pop, macro_son)
    pi0 = pi_0(d, macro_pop)
    list(
      cohort  = cg,
      P       = P,
      om      = round(om(P, pi0, t = 1), 3),
      d_prime = round(exp(d_prime(d, macro_pop, macro_son, t = 1)), 3)
    )
  })

cohort_tbl = purrr::map_dfr(cohort_results, function(res) {
  tibble(cohort = res$cohort, om_1 = res$om, d_prime_1 = res$d_prime)
})

cat("\n--- Cohort stationarity: OM(1) and d'(1) by birth cohort ---\n")
print(cohort_tbl)

cat("\n--- Macro transition matrices by cohort ---\n")
for (res in cohort_results) {
  cat("\nCohort:", res$cohort, "\n")
  print(round(res$P, 3))
}

################################################################################
####################### 9. REGIONAL SUMMARY TABLE #############################
################################################################################

regional_table = results_region |>
  select(region, n, sm, em, p_manual_fm_farming,
         p_farming_fm_farming, p_nonemp_fm_farming) |>
  arrange(desc(sm))

cat("\n--- Regional summary table ---\n")
print(
  knitr::kable(
    regional_table,
    col.names = c("Region", "N", "SM", "EM",
                  "P(Manual|Farm)", "P(Farm|Farm)", "P(Nonemp|Farm)"),
    digits  = 3,
    caption = "Regional mobility statistics (weighted macro 4x4)"
  )
)

# SM share of total mobility at t=1 — the paper's headline ratio
em_sm_ratio = om_total |>
  filter(t == 1, level == "macro") |>
  select(measure, est) |>
  pivot_wider(names_from = measure, values_from = est) |>
  mutate(sm_share = round(SM / (SM + EM), 3),
         om       = round(SM + EM, 3))

cat("\n--- SM share of OM at t=1 (macro) ---\n")
print(em_sm_ratio)

################################################################################
####################### 10. ROBUSTNESS: AGE 25-44 #############################
################################################################################

stats_main = compute_mobility_stats(data)

data_2544 = data |>
  filter((1940 - birthyr_son) >= 25, (1940 - birthyr_son) <= 44) |>
  mutate(w_atc_norm = w_atc_norm / sum(w_atc_norm) * n())

cat(sprintf("\n--- Age restriction: main n = %d, 25-44 n = %d (%.0f%% loss) ---\n",
            nrow(data), nrow(data_2544),
            (1 - nrow(data_2544) / nrow(data)) * 100))

stats_2544 = compute_mobility_stats(data_2544)

robustness_tbl = bind_rows(
  stats_main |> mutate(sample = "main (20-44)"),
  stats_2544 |> mutate(sample = "restricted (25-44)")
) |> select(sample, everything())

cat("\n--- Robustness: main vs 25-44 restricted sample ---\n")
print(robustness_tbl)

# Transition matrix for visual comparison
P_2544 = p_matrix(data_2544, macro_pop, macro_son)
cat("\nRestricted (25-44) macro transition matrix:\n")
print(round(P_2544, 3))
cat("\nMain macro transition matrix:\n")
print(round(P_macro_global, 3))


################################################################################
############################# 12. SAVE FIGURES ################################
################################################################################

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# Transition matrix heatmaps (wide: P matrix + pi_0 + pi*)
ggsave("output/figures/pmat_macro.png",        combined_plot_macro,      width = 14, height = 6,  dpi = 200)
ggsave("output/figures/pmat_meso.png",         combined_plot_meso,       width = 14, height = 6,  dpi = 200)
