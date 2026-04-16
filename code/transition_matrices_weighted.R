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

source("research/projects/aian-igm/code/utils.R")

data = read_csv("data/aian_weighted.csv") |>
  mutate(w_atc_norm = w_trim_norm)

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
    geom_text(aes(label = sprintf("%.2f\n[%.2f, %.2f]", est, lo, hi)),
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

# d(t), d'(t), AM convergence curves
plot_measure_curves = function(results_df) {
  ggplot(results_df, aes(x = t, y = est, color = measure, fill = measure)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, linetype = 0) +
    geom_line(linetype = 1) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c(d = "#332288", d_prime = "#009E73", AM = "darkred"),
                       labels = c(d = expression(log(d(t))),
                                  d_prime = expression(log(d*"'"*(t))),
                                  AM = expression(log(AM(pi[0], t))))) +
    scale_fill_manual(values = c(d = "#332288", d_prime = "#009E73", AM = "darkred"),
                      labels = c(d = expression(log(d(t))),
                                 d_prime = expression(log(d*"'"*(t))),
                                 AM = expression(log(AM(pi[0], t))))) +
    scale_y_continuous(breaks = c(0, -2, -4, -6, -8)) +
    labs(x = "Generation (t)", y = "log of Measure Value") +
    theme_minimal() +
    theme(legend.position = c(0.5, 0.05),
          legend.justification = c("right", "bottom"),
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          axis.line.x = element_line(linewidth = 0.5),
          axis.line.y = element_line(linewidth = 0.5),
          axis.text = element_text(size = 10),
          axis.ticks = element_line(size = 0.5)) +
    coord_cartesian(ylim = c(-8, 0))
}

# IM curves by origin category
plot_im_curves = function(im_df, color_values) {
  ggplot(im_df, aes(x = t, y = est, color = origin)) +
    geom_line(linetype = 1, linewidth = 1) +
    geom_point(size = 1.5) +
    scale_x_continuous(breaks = 0:5) +
    scale_y_continuous(breaks = c(0, -2, -4, -6, -8)) +
    scale_color_manual(name = "Occ. Origin", values = color_values) +
    labs(x = "Generation (t)", y = expression(log~IM(t,i))) +
    theme_minimal() +
    theme(legend.position = c(0.5, 0.05),
          legend.justification = c("right", "bottom"),
          legend.direction = "vertical",
          legend.title = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 12),
          axis.line.x = element_line(linewidth = 0.5),
          axis.line.y = element_line(linewidth = 0.5),
          axis.text = element_text(size = 10),
          axis.ticks = element_line(size = 0.5)) +
    coord_cartesian(ylim = c(-8, 0))
}

# IM min/max envelope plot
plot_im_minmax = function(im_df) {
  im_df |>
    group_by(t) |>
    summarise(est_max = max(est), lo_max = max(lo), hi_max = max(hi),
              est_min = min(est), lo_min = min(lo), hi_min = min(hi),
              .groups = "drop") |>
    pivot_longer(
      cols = c(est_max, lo_max, hi_max, est_min, lo_min, hi_min),
      names_to = c(".value", "type"),
      names_pattern = "(.*)_(max|min)") |>
    ggplot(aes(x = t, y = est, color = type, fill = type)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = type), alpha = 0.15, linetype = 0) +
    geom_point(size = 1.5) +
    geom_line(linetype = 1, linewidth = 1) +
    scale_x_continuous(breaks = 0:5) +
    scale_y_continuous(breaks = c(0, -2, -4, -6, -8)) +
    labs(x = "Generation (t)", y = expression(log~IM(t,i))) +
    theme_minimal() +
    theme(legend.position = c(0.5, 0.05),
          legend.justification = c("right", "bottom"),
          legend.direction = "vertical",
          legend.title = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 12),
          axis.line.x = element_line(linewidth = 0.5),
          axis.line.y = element_line(linewidth = 0.5),
          axis.text = element_text(size = 10),
          axis.ticks = element_line(size = 0.5)) +
    coord_cartesian(ylim = c(-8, 0))
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

desc_region   = weighted_prop_table(data, region) |>
  mutate(prop = round(prop, 1))
desc_meso_pop = weighted_prop_table(data, meso_pop)
desc_meso_son = weighted_prop_table(data, meso_son)
desc_educd    = weighted_prop_table(data, education)

################################################################################
############################# 4. COMPUTATION ###################################
################################################################################

ts = 0:4

# Transition matrices with bootstrap CIs
p_mat_macro = boot_pmatrix_ci(data, macro_pop, macro_son, R = 1000, conf = 0.95, .seed = 123)
p_mat_meso  = boot_pmatrix_ci(data, meso_pop, meso_son, R = 1000, conf = 0.95, .seed = 123)

# Initial and stationary distributions
pi0_vec_macro = pi_0(data, macro_pop)
pi0_vec_meso  = pi_0(data, meso_pop)
P_macro_global = p_matrix(data, macro_pop, macro_son, TRUE)
P_meso_global  = p_matrix(data, meso_pop, meso_son, TRUE)
verify_ergodic(P_macro_global, "macro global")
verify_ergodic(P_meso_global,  "meso global")
steady_macro  = pi_star(P_macro_global)
steady_meso   = pi_star(P_meso_global)

# d(t), d'(t), AM bootstrap
results_macro = boot_measures_by_t(data, macro_pop, macro_son, ts = ts, R = 1000, .seed = 123)
results_meso  = boot_measures_by_t(data, meso_pop, meso_son, ts = ts, R = 1000, .seed = 123)

# IM bootstrap
im_boot_macro = boot_im_by_t(data, macro_pop, macro_son, ts = 0:4, R = 1000, .seed = 123) |>
  mutate(origin = recode(origin,
                         farming = "Farming", manual = "Manual",
                         nonmanual = "Nonmanual", nilf = "Not working"))

im_boot_meso = boot_im_by_t(data, meso_pop, meso_son, ts = 0:4, R = 1000, .seed = 123) |>
  mutate(origin = recode(origin, farmer = "Farmer", unskilled = "Unskilled",
                         crafts = "Crafts", white_collar = "White Collar",
                         nilf = "Not working", farmworker = "Farmworker"),
         origin = factor(origin, levels = c("Farmer", "Farmworker", "Crafts",
                                            "Unskilled", "White Collar",
                                            "Not working")))

# Overall mobility bootstrap
macro_om = mobility_curve_with_boot(data, macro_pop, macro_son, ts = 1:6)
meso_om  = mobility_curve_with_boot(data, meso_pop, meso_son, ts = 1:6)

om_total = bind_rows(macro_om |> mutate(level = 1), meso_om |> mutate(level = 0)) |>
  mutate(level = factor(level, labels = c("meso", "macro"))) |>
  filter(measure != "SM")

################################################################################
############################# 5. PLOTTING ######################################
################################################################################

meso_level_order = c("nilf", "nonmanual", "crafts", "unskilled", "farmer", "farmworker")

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

## d(t), d'(t), AM convergence curves ----

measure_plot_macro = plot_measure_curves(results_macro)
measure_plot_meso  = plot_measure_curves(results_meso)
combined_measures  = measure_plot_macro + measure_plot_meso

## IM curves ----

im_macro_plot = plot_im_curves(im_boot_macro,
  c("Farming" = "forestgreen", "Manual" = "firebrick",
    "Nonmanual" = "steelblue", "Not working" = "darkorange"))

im_meso_plot = plot_im_curves(im_boot_meso,
  c("Farmer" = "darkgreen", "Farmworker" = "#90EE90",
    "Crafts" = "darkred", "Unskilled" = "pink",
    "White Collar" = "darkblue", "Not working" = "darkorange"))

im_combined = im_macro_plot + im_meso_plot + plot_layout(widths = c(2, 2))

## IM min/max envelopes ----

im_minmax_macro = plot_im_minmax(im_boot_macro)
im_minmax_meso  = plot_im_minmax(im_boot_meso)

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
                        labels = c("Exchange Mobility", "Overall Mobility"),
                        values = c(EM = "dashed", OM = "solid")) +
  scale_shape_manual(name = "Measure",
                     labels = c("Exchange Mobility", "Overall Mobility"),
                     values = c(EM = 17, OM = 16)) +
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
  p_nilf_given_farming = weighted.mean(
    farming_rows$macro_son == "nilf", farming_rows$w_atc_norm, na.rm = TRUE)

  tibble(
    sm                   = round(sm_val, 3),
    em                   = round(em_val, 3),
    om                   = round(om_val, 3),
    p_manual_fm_farming  = round(p_manual_given_farming, 3),
    p_farming_fm_farming = round(p_farming_given_farming, 3),
    p_nilf_fm_farming    = round(p_nilf_given_farming, 3)
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

cat("\n--- Macro transition matrices by cohort ---\n")
for (res in cohort_results) {
  cat("\nCohort:", res$cohort,
      "  OM(1) =", res$om,
      "  d'(1) =", res$d_prime, "\n")
  print(round(res$P, 3))
}
