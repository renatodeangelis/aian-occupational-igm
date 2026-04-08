library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(expm)
library(rlang)
library(purrr)
library(tigris)
library(jtools)
library(maps)
library(sf)
library(weights)

source("code/utils.R")

data = read_csv("data/aian_weighted.csv") |>
  mutate(meso_pop = if_else(meso_pop %in% c("prof", "clerical"), "white_collar", meso_pop),
         meso_son = if_else(meso_son %in% c("prof", "clerical"), "white_collar", meso_son))

## DATA DESCRIPTION

desc_main = data |>
  mutate(age_pop = 1940 - birthyr_pop,
         age_son = 1940 - birthyr_son) |>
  summarise(age_pop_mean = weighted.mean(age_pop, w_atc_norm),
            age_pop_sd = wtd.sd(age_pop, w_atc_norm),
            age_son_mean = weighted.mean(age_son, w_atc_norm),
            age_son_sd = wtd.sd(age_son, w_atc_norm))

desc_region = data |>
  group_by(region) |>
  summarise(wsum = sum(w_atc_norm), .groups = "drop") |>
  mutate(prop = round(wsum / sum(wsum) * 100, 1)) |>
  select(-wsum) |>
  arrange(region)

desc_meso_pop = data |>
  group_by(meso_pop) |>
  summarise(wsum = sum(w_atc_norm), .groups = "drop") |>
  mutate(prop = wsum / sum(wsum) * 100) |>
  select(-wsum) |>
  arrange(meso_pop)

desc_meso_son = data |>
  group_by(meso_son) |>
  summarise(wsum = sum(w_atc_norm), .groups = "drop") |>
  mutate(prop = wsum / sum(wsum) * 100) |>
  select(-wsum) |>
  arrange(meso_son)

desc_educd = data |>
  group_by(education) |>
  summarise(wsum = sum(w_atc_norm), .groups = "drop") |>
  mutate(prop = wsum / sum(wsum) * 100) |>
  select(-wsum) |>
  arrange(education)

################################################################################
############################### BASIC MEASURES #################################
################################################################################

macro_levels = unique(data$macro_pop)
meso_levels = unique(data$meso_pop)

################################################################################
############################### IMPLEMENTATION #################################
################################################################################

## MACRO
p_mat_macro = boot_pmatrix_ci(
  data, macro_pop, macro_son,
  R = 1000, conf = 0.95, .seed = 123)

g_macro = ggplot(
  p_mat_macro |>
    mutate(macro_pop = factor(macro_pop, levels = rev(unique(macro_pop)))),
  aes(x = macro_son, y = macro_pop, fill = est)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n[%.2f, %.2f]", est, lo, hi)),
            vjust = 0.3, size = 4) +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(x = "Son's occupation",
       y = "Father's occupation",
       fill = "Transition Prob.",
       title = expression(P)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 12),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

pi0_vec_macro = pi_0(data, macro_pop)
df_pi0_macro = tibble(
  father = names(pi0_vec_macro),
  pi0 = as.numeric(pi0_vec_macro)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g0_macro = ggplot(df_pi0_macro, aes(x = 1, y = father, fill = pi0)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi0)),
            size = 4, color = "black") +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(title = expression(pi[0])) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

steady_macro = pi_star(p_matrix(data, macro_pop, macro_son, TRUE))
df_pi_star_macro = tibble(
  father = names(steady_macro),
  pi_star = as.numeric(steady_macro)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g_star_macro = ggplot(df_pi_star_macro, aes(x = 1, y = father, fill = pi_star)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi_star)),
            size = 4, color = "black") +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(title = expression(pi^"*")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

combined_plot_macro = g_macro + g0_macro + g_star_macro + plot_layout(widths = c(6, 1, 1))

## MESO
p_mat_meso = boot_pmatrix_ci(
  data, meso_pop, meso_son,
  R = 1000, conf = 0.95, .seed = 123)

g_meso = ggplot(
  p_mat_meso |>
    mutate(meso_pop = factor(meso_pop, levels = c("unemp", "white_collar",
                                                  "crafts", "unskilled",
                                                  "farmer", "farmworker")),
           meso_son = factor(meso_son, levels = rev(c("unemp", "white_collar",
                                                      "crafts", "unskilled",
                                                      "farmer", "farmworker")))),
  aes(x = meso_son, y = meso_pop, fill = est)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n[%.2f, %.2f]", est, lo, hi)),
            vjust = 0.3, size = 3) +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(x = "Son's occupation",
       y = "Father's occupation",
       fill = "Transition Prob.",
       title = "P") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

pi0_vec_meso = pi_0(data, meso_pop)
df_pi0_meso = tibble(
  father = names(pi0_vec_meso),
  pi0 = as.numeric(pi0_vec_meso)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g0_meso = ggplot(df_pi0_meso |>
                   mutate(father = factor(father, levels = c("unemp", "white_collar",
                                                                 "crafts", "unskilled",
                                                                 "farmer", "farmworker"))),
                          aes(x = 1, y = father, fill = pi0)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi0)),
            size = 4, color = "black") +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(title = expression(pi[0])) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

steady_meso = pi_star(p_matrix(data, meso_pop, meso_son, TRUE))
df_pi_star_meso = tibble(
  father = names(steady_meso),
  pi_star = as.numeric(steady_meso)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g_star_meso = ggplot(df_pi_star_meso |>
                       mutate(father = factor(father, levels = c("unemp", "white_collar",
                                                                     "crafts", "unskilled",
                                                                     "farmer", "farmworker"))),
                     aes(x = 1, y = father, fill = pi_star)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi_star)),
            size = 4, color = "black") +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(title = expression(pi^"*")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

combined_plot_meso = g_meso + g0_meso + g_star_meso + plot_layout(widths = c(6, 1, 1))

################################################################################

## d(t), d'(t), AM CURVES

ts = 0:4

results_macro = boot_measures_by_t(data, macro_pop, macro_son, ts = ts, R = 1000, .seed = 123)

measure_plot_macro = ggplot(results_macro, aes(x = t, y = est, color = measure, fill = measure)) +
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

results_meso = boot_measures_by_t(data, meso_pop, meso_son, ts = ts, R = 1000, .seed = 123)

measure_plot_meso = ggplot(results_meso, aes(x = t, y = est, color = measure, fill = measure)) +
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

combined_measures = measure_plot_macro + measure_plot_meso

## IM CURVES

# MACRO
im_boot_macro = boot_im_by_t(data, macro_pop, macro_son, ts = 0:4, R = 1000, .seed = 123) |>
  dplyr::mutate(origin = dplyr::recode(origin,
                                       farming = "Farming", manual = "Manual",
                                       nonmanual = "Nonmanual", unemp = "Not working"))

im_macro_plot = ggplot(im_boot_macro, aes(x = t, y = est, color = origin)) +
#  geom_ribbon(aes(ymin = lo, ymax = hi, fill = origin), alpha = 0.15, linetype = 0) +
  geom_line(linetype = 1, linewidth = 1) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8)) +
  scale_color_manual(name = "Occ. Origin",
                     values = c("Farming" = "forestgreen",
                                "Manual" = "firebrick",
                                "Nonmanual" = "steelblue",
                                "Not working" = "darkorange")) +
  labs(x = "Generation (t)", y = expression(log~IM(t,i))) +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.05),
        legend.justification = c("right","bottom"),
        legend.direction = "vertical",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(ylim = c(-8, 0))

# MESO
im_boot_meso = boot_im_by_t(data, meso_pop, meso_son, ts = 0:4, R = 1000, .seed = 123) |>
  mutate(origin = recode(origin, farmer = "Farmer", unskilled = "Unskilled",
                         crafts = "Crafts", white_collar = "White Collar",
                         unemp = "Not working", farmworker = "Farmworker"),
         origin = factor(origin, levels = c("Farmer", "Farmworker", "Crafts",
                                            "Unskilled", "White Collar",
                                            "Not working")))

im_meso_plot = ggplot(im_boot_meso, aes(x = t, y = est, color = origin)) +
#  geom_ribbon(aes(ymin = lo, ymax = hi, fill = origin), alpha = 0.15, linetype = 0) +
  geom_line(linetype = 1, linewidth = 1) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8)) +
  scale_color_manual(name = "Occ. Origin",
                     values = c("Farmer" = "darkgreen",
                                "Farmworker" = "#90EE90",
                                "Crafts" = "darkred",
                                "Unskilled" = "pink",
                                "White Collar" = "darkblue",
                                "Not working" = "darkorange")) +
  labs(x = "Generation (t)", y = expression(log~IM(t,i))) +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.05),
        legend.justification = c("right","bottom"),
        legend.direction = "vertical",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(ylim = c(-8, 0))

im_combined = im_macro_plot + im_meso_plot +
  plot_layout(widths = c(2, 2))

# MIN AND MAX

im_minmax_macro = im_boot_macro |>
  group_by(t) |>
  summarise(est_max = max(est),
            lo_max = max(lo),
            hi_max = max(hi),
            est_min = min(est),
            lo_min = min(lo),
            hi_min = min(hi),
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
        legend.justification = c("right","bottom"),
        legend.direction = "vertical",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(ylim = c(-8, 0))

im_minmax_meso = im_boot_meso |>
  group_by(t) |>
  summarise(est_max = max(est),
            lo_max = max(lo),
            hi_max = max(hi),
            est_min = min(est),
            lo_min = min(lo),
            hi_min = min(hi),
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
        legend.justification = c("right","bottom"),
        legend.direction = "vertical",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(ylim = c(-8, 0))


## OVERALL MOBILITY

macro_om = mobility_curve_with_boot(data, macro_pop, macro_son, ts = 0:6)
meso_om = mobility_curve_with_boot(data, meso_pop, meso_son, ts = 0:6)

om_total = bind_rows(macro_om |> mutate(level = 1), meso_om |> mutate(level = 0)) |>
  mutate(level = factor(level, labels = c("meso", "macro"))) |>
  filter(measure != "SM")

om_plot = ggplot(om_total, aes(x = t, y = est)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = level,
                  group = interaction(level, measure)),
              alpha = 0.2, color = NA) +
  geom_line(aes(color = level, linetype = measure,
                group = interaction(level, measure)), linewidth = 1) +
  geom_point(aes(shape = measure, color = level), size = 2.5) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6)) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)) +
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
  labs(x = "Generation (t)",
       y = "Mobility") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(xlim = c(0, 6),
                  ylim = c(0, 0.7)); om_plot

## COUNTY ESTIMATION

compute_mobility_stats = function(df) {

  p_upward = df |>
    filter(macro_pop != "nonmanual") |>
    mutate(count = (macro_son == "nonmanual") |
             (macro_pop == "unemp" & macro_son != "unemp")) |>
    summarise(tot = sum(w_atc_norm), up = sum(count * w_atc_norm)) |>
    mutate(prop = up / tot) |>
    select(prop)

  p_downward = df |>
    filter(macro_pop == "nonmanual" | meso_pop == "farmer" | meso_pop == "crafts") |>
    mutate(count = (macro_pop == "nonmanual" & macro_son != "nonmanual") |
             (meso_pop == "farmer" & meso_son %in% c("unskilled", "farmworker", "unemp")) |
             (meso_pop == "crafts" & meso_son %in% c("unskilled", "farmworker", "unemp"))) |>
    summarise(tot = sum(w_atc_norm), down = sum(count * w_atc_norm)) |>
    mutate(prop = down / tot) |>
    select(prop)

  d1 = d_prime(df, macro_pop, macro_son, t = 1)

  P = p_matrix(df, macro_pop, macro_son)
  pi0 = pi_0(df, macro_pop)

  P_unweighted = p_matrix_unweighted(df, macro_pop, macro_son)
  pi0_unweighted = pi_0_unweighted(df, macro_pop)

  om_weighted = om(P, pi0, t = 0)
  om_unweighted = om(P_unweighted, pi0_unweighted, t = 0)

  tibble(p_upward = round(p_upward$prop, 2),
         p_downward = round(p_downward$prop, 2),
         d_prime_1 = round(exp(d1), 2),
         om_1 = round(om_weighted, 2),
         om_1_unweighted = round(om_unweighted, 2))
}

results_region = data |>
  group_by(region) |>
  count() |>
  mutate(stats = purrr::map(region, function(re) {
    df_sub = data |>
      filter(region == re)
    compute_mobility_stats(df_sub)})) |>
  unnest(stats)

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

om_1_plot = ggplot() +
  geom_sf(data = states_sf, aes(fill = om_1_unweighted), color = "white", size = 0) +
  geom_sf(data = regions_sf, fill = NA, color = "black", size = 0.8) +
  geom_sf_text(data = centroids, aes(label = om_1_unweighted),
               size = 4, color = "black", fontface = "bold") +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, face = "bold")) +
  labs(title = "OM(0) by Region")

d_1_prime_plot = ggplot() +
  geom_sf(data = states_sf, aes(fill = d_prime_1), color = "white", size = 0) +
  geom_sf(data = regions_sf, fill = NA, color = "black", size = 0.8) +
  geom_sf_text(data = centroids, aes(label = d_prime_1),
               size = 4, color = "black", fontface = "bold") +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, face = "bold")) +
  labs(title = "d'(1) by Region")

upward_downward_plot = ggplot() +
  geom_sf(data = states_sf, aes(fill = round(p_upward/p_downward, 2)), color = "white", size = 0) +
  geom_sf(data = regions_sf, fill = NA, color = "black", size = 0.8) +
  geom_sf_text(data = centroids, aes(label = round(p_upward/p_downward, 2)),
               size = 4, color = "black", fontface = "bold") +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, face = "bold")) +
  labs(title = "P(upward) / P(downward) by Region")
