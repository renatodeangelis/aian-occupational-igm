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

data = read_csv("data/aian_weighted.csv") |>
  mutate(w_atc_norm = w_trim_norm)

aian_full = readRDS("data/aian_full.rds")

macro_levels = unique(data$macro_pop)
meso_levels = unique(data$meso_pop)

################################################################################
############################# 2. PLOT HELPERS ##################################
################################################################################

# Weighted proportion table for a single grouping variable
weighted_prop_table = function(data, var) {
  var_sym = ensym(var)
  data |>
    group_by(!!var_sym) |>
    summarise(wsum = sum(w_atc_norm), .groups = "drop") |>
    mutate(prop = wsum / sum(wsum) * 100) |>
    select(-wsum) |>
    arrange(!!var_sym)
}

# Transition matrix heatmap with bootstrap CIs
plot_pmat_heatmap = function(boot_df, dad_var, son_var,
                             levels = NULL, text_size = 4, title_expr = "P") {
  dad_sym = ensym(dad_var)
  son_sym = ensym(son_var)

  plot_df = boot_df
  if (!is.null(levels)) {
    plot_df = plot_df |>
      mutate(!!dad_sym := factor(!!dad_sym, levels = levels),
             !!son_sym := factor(!!son_sym, levels = rev(levels)))
  } else {
    plot_df = plot_df |>
      mutate(!!dad_sym := factor(!!dad_sym, levels = rev(unique(!!dad_sym))))
  }

  ggplot(plot_df, aes(x = !!son_sym, y = !!dad_sym, fill = est)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f\n(%.3f)", est, se)),
              vjust = 0.3, size = text_size) +
    scale_fill_gradient(low = "lightyellow", high = "firebrick") +
    labs(x = "Son's occupation", y = "Father's occupation",
         fill = "Transition Prob.", title = title_expr) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 10),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
}

# Single-column heatmap for pi_0 or pi_star vectors
plot_pi_column = function(vec, title_expr, levels = NULL) {
  df = tibble(father = names(vec), value = as.numeric(vec))
  if (!is.null(levels)) {
    df = df |> mutate(father = factor(father, levels = levels))
  } else {
    df = df |> mutate(father = factor(father, levels = rev(unique(father))))
  }

  ggplot(df, aes(x = 1, y = father, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 4, color = "black") +
    scale_fill_gradient(low = "lightyellow", high = "firebrick") +
    labs(title = title_expr) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
}

# Region choropleth map
plot_region_map = function(states_sf, regions_sf, centroids,
                           fill_expr, label_expr, title) {
  ggplot() +
    geom_sf(data = states_sf, aes(fill = {{ fill_expr }}), color = "white", size = 0) +
    geom_sf(data = regions_sf, fill = NA, color = "black", size = 0.8) +
    geom_sf_text(data = centroids, aes(label = {{ label_expr }}),
                 size = 4, color = "black", fontface = "bold") +
    scale_fill_gradient(low = "lightyellow", high = "firebrick") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 15, face = "bold")) +
    labs(title = title)
}

################################################################################
############################# 3. DESCRIPTIVES ##################################
################################################################################

desc_main = data |>
  mutate(age_pop = 1940 - birthyr_pop,
         age_son = 1940 - birthyr_son) |>
  summarise(age_pop_mean = weighted.mean(age_pop, w_atc_norm),
            age_pop_sd = wtd.sd(age_pop, w_atc_norm),
            age_son_mean = weighted.mean(age_son, w_atc_norm),
            age_son_sd = wtd.sd(age_son, w_atc_norm))

desc_region = weighted_prop_table(data, region) |>
  mutate(prop = round(prop, 1))
desc_macro_pop = weighted_prop_table(data, macro_pop)
desc_macro_pop_alt = weighted_prop_table(data, macro_pop_alt)
desc_macro_son = weighted_prop_table(data, macro_son)
desc_macro_son_alt = weighted_prop_table(data, macro_son_alt)
desc_educd = weighted_prop_table(data, education)

################################################################################
############################# 4. COMPUTATION ###################################
################################################################################

ts = 1:4

# Bootstrap cache — set force_rerun = TRUE to discard cache and re-estimate.
# Delete individual files from cache/ to selectively re-run one bootstrap.
# Note: cache is not auto-invalidated if aian_weighted.csv changes; delete
# cache/ manually after re-running weighting.R.
force_rerun <- FALSE
cache_dir   <- "cache"
dir.create(cache_dir, showWarnings = FALSE)

cache_load <- function(name, expr, force = force_rerun) {
  path <- file.path(cache_dir, paste0(name, ".rds"))
  if (!force && file.exists(path)) {
    message("Loading cached: ", name)
    return(readRDS(path))
  }
  result <- eval(expr, envir = parent.frame())
  saveRDS(result, path)
  result
}

# Transition matrices with bootstrap CIs
p_mat_macro = cache_load("p_mat_macro", quote(
  boot_pmatrix_ci(data, macro_pop, macro_son,
                  df_linked = data, df_full = aian_full,
                  R = 500, .seed = 123)
))
p_mat_meso  = cache_load("p_mat_meso", quote(
  boot_pmatrix_ci(data, meso_pop, meso_son,
                  df_linked = data, df_full = aian_full,
                  R = 500, .seed = 123)
))

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
))
meso_om  = cache_load("meso_om", quote(
  mobility_curve_with_boot(data, meso_pop, meso_son,
                            df_linked = data, df_full = aian_full,
                            ts = 1:4, R = 500, .seed = 123)
))

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

macro_level_order = c("nonemp", "nonmanual", "manual", "farmer")
meso_level_order = c("nonemp", "nonmanual", "crafts", "unskilled", "farmer", "farmworker")

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

## Overall mobility ----

om_plot = ggplot(om_total, aes(x = t, y = est)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = level,
                  group = interaction(level, measure)),
              alpha = 0.2, color = NA) +
  geom_line(aes(color = level, linetype = measure,
                group = interaction(level, measure)), linewidth = 1) +
  geom_point(aes(shape = measure, color = level), size = 2.5) +
  scale_x_continuous(breaks = 0:6) +
  scale_y_continuous(breaks = seq(0, 0.7, 0.1)) +
  scale_linetype_manual(name = "Measure",
                        labels = c(EM = "Exchange Mobility", SM = "Structural Mobility"),
                        values = c(EM = "dashed", SM = "solid")) +
  scale_shape_manual(name = "Measure",
                     labels = c(EM = "Exchange Mobility", SM = "Structural Mobility"),
                     values = c(EM = 17, SM = 16)) +
  scale_fill_manual(name = "Level",
                    labels = c("Meso (6 categories)", "Macro (4 categories)"),
                    values = c("meso" = "#009E73", "macro" = "#D55E00")) +
  scale_color_manual(name = "Level",
                     labels = c("Meso (6 categories)", "Macro (4 categories)"),
                     values = c("meso" = "#009E73", "macro" = "#D55E00")) +
  labs(x = "Generation (t)", y = "Mobility") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(xlim = c(0, 6), ylim = c(0, 0.7))
om_plot

################################################################################
############################# 6. REGIONAL ######################################
################################################################################

compute_mobility_stats = function(df) {
  df = df |> mutate(w_atc_norm = w_atc_norm / sum(w_atc_norm) * n())

  P      = p_matrix(df, macro_pop, macro_son)
  pi0    = pi_0(df, macro_pop)
  pistar = pi_star(P)

  om_val = om(P, pi0, t = 0)
  sm_val = sm(P, pi0, t = 0)
  em_val = om_val - sm_val

  # Farming-family decomposition
  farming_rows = df |> filter(macro_pop == "farming")
  p_manual_given_farming = weighted.mean(
    farming_rows$macro_son == "manual", farming_rows$w_atc_norm, na.rm = TRUE)
  p_farming_given_farming = weighted.mean(
    farming_rows$macro_son == "farming", farming_rows$w_atc_norm, na.rm = TRUE)
  p_nonemp_given_farming = weighted.mean(
    farming_rows$macro_son == "nonemp", farming_rows$w_atc_norm, na.rm = TRUE)

  # Nonemp-family decomposition
  nonemp_rows = df |> filter(macro_pop == "nonemp")
  p_nonemp_given_nonemp = weighted.mean(
    nonemp_rows$macro_son == "nonemp", nonemp_rows$w_atc_norm, na.rm = TRUE)

  # Marginal son distribution
  p_son_nonemp    = weighted.mean(df$macro_son == "nonemp",    df$w_atc_norm, na.rm = TRUE)
  p_son_nonmanual = weighted.mean(df$macro_son == "nonmanual", df$w_atc_norm, na.rm = TRUE)

  tibble(
    sm                   = round(sm_val, 3),
    em                   = round(em_val, 3),
    om                   = round(om_val, 3),
    p_manual_fm_farming  = round(p_manual_given_farming, 3),
    p_farming_fm_farming = round(p_farming_given_farming, 3),
    p_nonemp_fm_farming  = round(p_nonemp_given_farming, 3),
    p_nonemp_fm_nonemp   = round(p_nonemp_given_nonemp, 3),
    p_son_nonemp         = round(p_son_nonemp, 3),
    p_son_nonmanual      = round(p_son_nonmanual, 3)
  )
}

results_region = data |>
  group_by(region) |>
  group_modify(~ compute_mobility_stats(.x)) |>
  ungroup() |>
  left_join(count(data, region), by = "region")

state_regions = tibble(
  state_name = tolower(c(
    "alabama","arizona","arkansas","california","colorado","connecticut","delaware",
    "florida","georgia","idaho","illinois","indiana","iowa","kansas","kentucky",
    "louisiana","maine","maryland","massachusetts","michigan","minnesota",
    "mississippi","missouri","montana","nebraska","nevada","new hampshire",
    "new jersey","new mexico","new york","north carolina","north dakota","ohio",
    "oklahoma","oregon","pennsylvania","rhode island","south carolina","south dakota",
    "tennessee","texas","utah","vermont","virginia","washington","west virginia","wisconsin","wyoming")),
  statefip_1940 = c(
    1,4,5,6,8,9,10,
    12,13,16,17,18,19,20,21,
    22,23,24,25,26,27,28,29,
    30,31,32,33,34,35,36,37,
    38,39,40,41,42,44,45,46,
    47,48,49,50,51,53,54,55,56),
  region = assign_region(statefip_1940))

states_sf = st_as_sf(map("state", plot = FALSE, fill = TRUE)) |>
  rename(state_name = ID) |>
  left_join(state_regions, by = "state_name") |>
  st_make_valid() |>
  st_buffer(dist = 0) |>
  left_join(results_region, by = "region")

sf::sf_use_s2(FALSE)

regions_sf = states_sf |>
  filter(!is.na(region)) |>
  group_by(region) |>
  summarize(geom = st_union(geom), .groups = "drop")

centroids = st_centroid(regions_sf) |>
  left_join(results_region, by = "region")

## Regional maps ----

sm_plot = plot_region_map(states_sf, regions_sf, centroids,
                            sm, sm,
                            "Structural Mobility by Region")

em_plot = plot_region_map(states_sf, regions_sf, centroids,
                                em, em,
                                "Exchange Mobility by Region")

p_manual_plot = plot_region_map(states_sf, regions_sf, centroids,
                                p_manual_fm_farming, p_manual_fm_farming,
                                "P(Manual | Father Farming) by Region")

p_farming_plot = plot_region_map(states_sf, regions_sf, centroids,
                                 p_farming_fm_farming, p_farming_fm_farming,
                                 "P(Farming | Father Farming) by Region")

p_nonemp_plot = plot_region_map(states_sf, regions_sf, centroids,
                                p_nonemp_fm_nonemp, p_nonemp_fm_nonemp,
                                "P(Nonemp | Father Nonemp) by Region")

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
##################### 8. ALT MATRIX (ATTACHMENT-FIRST) ########################
################################################################################

# Men with a valid occ code but no recorded labor market activity are
# reclassified to nonemp, applied symmetrically to fathers and sons.
# data_alt renames the alt columns into the main column slots so all
# existing functions (pi_0, p_matrix, compute_mobility_stats) work unchanged.
data_alt = data |>
  select(-macro_pop, -macro_son, -meso_pop, -meso_son) |>
  rename(
    macro_pop = macro_pop_alt,
    macro_son = macro_son_alt,
    meso_pop  = meso_pop_alt,
    meso_son  = meso_son_alt
  )

## Alt transition matrices ----

p_mat_macro_alt = cache_load("p_mat_macro_alt", quote(
  boot_pmatrix_ci(data_alt, macro_pop, macro_son,
                  df_linked = data_alt, df_full = aian_full,
                  R = 500, .seed = 123)
))

pi0_vec_macro_alt  = pi_0(data_alt, macro_pop)
P_macro_alt_global = p_matrix(data_alt, macro_pop, macro_son)
verify_ergodic(P_macro_alt_global, "macro alt global")
steady_macro_alt   = pi_star(P_macro_alt_global)

## Mobility stats: main vs alt ----

stats_main = compute_mobility_stats(data)
stats_alt  = compute_mobility_stats(data_alt)

comparison_tbl = bind_rows(
  stats_main |> mutate(matrix = "main (occ-first)"),
  stats_alt  |> mutate(matrix = "alt (attachment-first)")
) |> select(matrix, everything())

cat("\n--- Main vs Alt: global mobility statistics ---\n")
print(comparison_tbl)

## Regional comparison ----

results_region_alt = data_alt |>
  group_by(region) |>
  group_modify(~ compute_mobility_stats(.x)) |>
  ungroup() |>
  left_join(count(data_alt, region), by = "region")

comparison_region_tbl = bind_rows(
  results_region     |> mutate(matrix = "main"),
  results_region_alt |> mutate(matrix = "alt")
) |> select(matrix, region, n, sm, em, p_nonemp_fm_farming, p_manual_fm_farming,
            p_farming_fm_farming)

cat("\n--- Main vs Alt: regional mobility statistics ---\n")
print(comparison_region_tbl)

## Alt regional maps ----

states_sf_alt = st_as_sf(map("state", plot = FALSE, fill = TRUE)) |>
  rename(state_name = ID) |>
  left_join(state_regions, by = "state_name") |>
  st_make_valid() |>
  st_buffer(dist = 0) |>
  left_join(results_region_alt, by = "region")

regions_sf_alt = states_sf_alt |>
  filter(!is.na(region)) |>
  group_by(region) |>
  summarize(geom = st_union(geom), .groups = "drop")

centroids_alt = st_centroid(regions_sf_alt) |>
  left_join(results_region_alt, by = "region")

sm_plot_alt = plot_region_map(states_sf_alt, regions_sf_alt, centroids_alt,
                              sm, sm,
                              "Structural Mobility by Region (Alt)")

em_plot_alt = plot_region_map(states_sf_alt, regions_sf_alt, centroids_alt,
                              em, em,
                              "Exchange Mobility by Region (Alt)")

p_manual_plot_alt = plot_region_map(states_sf_alt, regions_sf_alt, centroids_alt,
                                    p_manual_fm_farming, p_manual_fm_farming,
                                    "P(Manual | Father Farming) by Region (Alt)")

p_farming_plot_alt = plot_region_map(states_sf_alt, regions_sf_alt, centroids_alt,
                                     p_farming_fm_farming, p_farming_fm_farming,
                                     "P(Farming | Father Farming) by Region (Alt)")

p_nonemp_plot_alt = plot_region_map(states_sf_alt, regions_sf_alt, centroids_alt,
                                    p_nonemp_fm_nonemp, p_nonemp_fm_nonemp,
                                    "P(Nonemp | Father Nonemp) by Region (Alt)")

## Alt heatmap plots ----

g_macro_alt      = plot_pmat_heatmap(p_mat_macro_alt, macro_pop, macro_son,
                                     title_expr = expression(P^{alt}))
g0_macro_alt     = plot_pi_column(pi0_vec_macro_alt, expression(pi[0]^{alt}))
g_star_macro_alt = plot_pi_column(steady_macro_alt, expression(pi^{"*,alt"}))
combined_plot_macro_alt = g_macro_alt + g0_macro_alt + g_star_macro_alt +
  plot_layout(widths = c(6, 1, 1))

combined_main_vs_alt = combined_plot_macro / combined_plot_macro_alt
combined_main_vs_alt

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
####################### 11. UNIDIFF BY REGION #################################
################################################################################

# Tests whether regional variation in father-son association is a difference
# in degree (phi_k varies, pattern fixed) or kind (pattern changes).
# Macro 4x4 used to avoid sparse cells in smaller regions.
# Weighted cell frequencies passed as response; Poisson family is valid with
# fractional counts — coefficient estimates are consistent.

library(gnm)

region_long = data |>
  count(region, macro_pop, macro_son, wt = w_atc_norm, name = "freq") |>
  mutate(
    region    = factor(region),
    macro_pop = factor(macro_pop, levels = macro_order),
    macro_son = factor(macro_son, levels = macro_order)
  )

# Conditional independence baseline (no association layer)
unidiff_null = gnm(
  freq ~ region + macro_pop + macro_son,
  data = region_long, family = poisson, trace = FALSE
)

# UNIDIFF: single multiplier phi_k scales a common association pattern per region
unidiff_mod = gnm(
    freq ~ region + macro_pop + macro_son +                                  
           Mult(region, macro_pop:macro_son),
    data = region_long, family = poisson, trace = FALSE                      
  )                                                         

# Saturated association (region-specific full interaction — no constraint)
unidiff_sat = gnm(
  freq ~ region + macro_pop + macro_son + region:(macro_pop:macro_son),
  data = region_long, family = poisson, trace = FALSE
)

cat("\n--- UNIDIFF model fit ---\n")
cat(sprintf("Null (indep):  df = %d,  deviance = %.2f\n",
            df.residual(unidiff_null), deviance(unidiff_null)))
cat(sprintf("UNIDIFF:       df = %d,  deviance = %.2f\n",
            df.residual(unidiff_mod),  deviance(unidiff_mod)))
cat(sprintf("Saturated:     df = %d,  deviance = %.2f\n",
            df.residual(unidiff_sat),  deviance(unidiff_sat)))

# Extract phi_k (the UNIDIFF multipliers, one per region)
phi_idx = pickCoef(unidiff_mod, "Exp\\(region\\)")
phi_est = coef(unidiff_mod)[phi_idx]
phi_se  = sqrt(diag(vcov(unidiff_mod)))[phi_idx]

phi_tbl = tibble(
  region = levels(region_long$region),
  phi    = round(exp(phi_est), 3),   # exponentiate: phi > 1 = stronger association
  se     = round(phi_se, 3)
)

cat("\n--- UNIDIFF phi parameters by region (ref = first level) ---\n")
print(phi_tbl)

################################################################################
############################# 12. SAVE FIGURES ################################
################################################################################

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# Transition matrix heatmaps (wide: P matrix + pi_0 + pi*)
ggsave("output/figures/pmat_macro.pdf",        combined_plot_macro,      width = 14, height = 6)
ggsave("output/figures/pmat_meso.pdf",         combined_plot_meso,       width = 14, height = 6)
ggsave("output/figures/pmat_alt_macro.pdf",    combined_plot_macro_alt,  width = 14, height = 6)
ggsave("output/figures/pmat_main_vs_alt.pdf",  combined_main_vs_alt,     width = 14, height = 12)

# EM/SM mobility curves
ggsave("output/figures/om_plot.pdf",           om_plot,         width = 8,  height = 6)

# Regional maps
ggsave("output/figures/map_sm.pdf",            sm_plot,         width = 10, height = 6)
ggsave("output/figures/map_em.pdf",            em_plot,         width = 10, height = 6)
ggsave("output/figures/map_p_manual.pdf",      p_manual_plot,   width = 10, height = 6)
ggsave("output/figures/map_p_farming.pdf",     p_farming_plot,  width = 10, height = 6)
ggsave("output/figures/map_p_nonemp.pdf",      p_nonemp_plot,   width = 10, height = 6)

# Alt regional maps
ggsave("output/figures/map_sm_alt.pdf",        sm_plot_alt,        width = 10, height = 6)
ggsave("output/figures/map_em_alt.pdf",        em_plot_alt,        width = 10, height = 6)
ggsave("output/figures/map_p_manual_alt.pdf",  p_manual_plot_alt,  width = 10, height = 6)
ggsave("output/figures/map_p_farming_alt.pdf", p_farming_plot_alt, width = 10, height = 6)
ggsave("output/figures/map_p_nonemp_alt.pdf",  p_nonemp_plot_alt,  width = 10, height = 6)
