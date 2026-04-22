################################################################################
# generate_pmat_figures.R
# Generates the four transition matrix heatmaps (macro, meso, macro-alt,
# meso-alt) and saves them as PNG files for Overleaf.
# Run from the repo root: Rscript code/generate_pmat_figures.R
################################################################################

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(expm)
library(purrr)

source("code/utils.R")

data      <- read_csv("data/aian_weighted.csv") |> mutate(w_atc_norm = w_trim_norm)
aian_full <- readRDS("data/aian_full.rds")

# ── Display labels for occupation categories ───────────────────────────────────
occ_labels <- c(
  farming    = "Farming",
  farmer     = "Farming",      # meso uses "farmer"; macro uses "farming"
  farmworker = "Farmworker",
  nonemp     = "Non-employed",
  nonmanual  = "Non-manual",
  manual     = "Manual",
  crafts     = "Crafts",
  unskilled  = "Unskilled"
)

recode_occ_df <- function(df, ...) {
  vars <- rlang::ensyms(...)
  for (v in vars) {
    df <- df |> mutate(!!v := dplyr::recode(as.character(!!v), !!!occ_labels))
  }
  df
}

recode_occ_vec <- function(vec) {
  setNames(as.numeric(vec), dplyr::recode(names(vec), !!!occ_labels))
}

# ── Level orders (display names) ───────────────────────────────────────────────
macro_level_order <- occ_labels[c("nonemp", "nonmanual", "manual", "farming")]
canonical_meso    <- c("nonemp", "nonmanual", "crafts", "unskilled", "farmworker", "farmer")
meso_level_order  <- occ_labels[canonical_meso]

# ── Plot helpers ───────────────────────────────────────────────────────────────
plot_pmat_heatmap <- function(boot_df, dad_var, son_var,
                              levels = NULL, text_size = 5, title_str = "P") {
  dad_sym <- rlang::ensym(dad_var)
  son_sym <- rlang::ensym(son_var)
  plot_df <- boot_df
  if (!is.null(levels)) {
    plot_df <- plot_df |>
      mutate(!!dad_sym := factor(!!dad_sym, levels = levels),
             !!son_sym := factor(!!son_sym, levels = rev(levels)))
  } else {
    plot_df <- plot_df |>
      mutate(!!dad_sym := factor(!!dad_sym, levels = rev(unique(!!dad_sym))))
  }
  ggplot(plot_df, aes(x = !!son_sym, y = !!dad_sym, fill = est)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f\n(%.3f)", est, se)),
              vjust = 0.3, size = text_size) +
    scale_fill_gradient(low = "lightyellow", high = "firebrick") +
    labs(x = "Son's occupation", y = "Father's occupation",
         fill = "Transition Prob.", title = title_str) +
    theme_minimal() +
    theme(axis.text.x       = element_text(angle = 45, hjust = 1, size = 13),
          axis.text.y       = element_text(angle = 45, hjust = 1, size = 13),
          axis.ticks.x      = element_blank(),
          axis.ticks.y      = element_blank(),
          axis.title        = element_text(size = 13),
          legend.position   = "bottom",
          legend.text       = element_text(size = 11),
          legend.title      = element_text(size = 12),
          plot.title        = element_text(hjust = 0.5, size = 16),
          panel.grid        = element_blank())
}

plot_pi_column <- function(vec, title_str, levels = NULL) {
  df <- tibble(father = names(vec), value = as.numeric(vec))
  if (!is.null(levels)) {
    df <- df |> mutate(father = factor(father, levels = levels))
  } else {
    df <- df |> mutate(father = factor(father, levels = rev(unique(father))))
  }
  ggplot(df, aes(x = 1, y = father, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 5, color = "black") +
    scale_fill_gradient(low = "lightyellow", high = "firebrick") +
    labs(title = title_str) +
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(), axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(), legend.position = "none",
          plot.title   = element_text(hjust = 0.5, size = 16),
          panel.grid   = element_blank())
}

# ── Cache helpers ──────────────────────────────────────────────────────────────
force_rerun <- FALSE
cache_dir   <- "cache"
dir.create(cache_dir, showWarnings = FALSE)

cache_load <- function(name, expr, force = force_rerun) {
  path <- file.path(cache_dir, paste0(name, ".rds"))
  if (!force && file.exists(path)) {
    message("Loading cached: ", name); return(readRDS(path))
  }
  result <- eval(expr, envir = parent.frame())
  saveRDS(result, path)
  result
}

# ── Bootstrap transition matrices ──────────────────────────────────────────────
R_slides <- 50  # fast for presentation; bump to 500 for paper

message("Running macro bootstrap …")
p_mat_macro <- cache_load("p_mat_macro_slides", quote(
  boot_pmatrix_ci(data, macro_pop, macro_son,
                  df_linked = data, df_full = aian_full, R = R_slides, .seed = 123)
))

message("Running meso bootstrap …")
p_mat_meso <- cache_load("p_mat_meso_slides", quote(
  boot_pmatrix_ci(data, meso_pop, meso_son,
                  df_linked = data, df_full = aian_full, R = R_slides, .seed = 123)
))

# ── Alt dataset ────────────────────────────────────────────────────────────────
data_alt <- data |>
  select(-macro_pop, -macro_son, -meso_pop, -meso_son) |>
  rename(macro_pop = macro_pop_alt, macro_son = macro_son_alt,
         meso_pop  = meso_pop_alt,  meso_son  = meso_son_alt)

message("Running macro-alt bootstrap …")
p_mat_macro_alt <- cache_load("p_mat_macro_alt_slides", quote(
  boot_pmatrix_ci(data_alt, macro_pop, macro_son,
                  df_linked = data_alt, df_full = aian_full, R = R_slides, .seed = 123)
))

message("Running meso-alt bootstrap …")
p_mat_meso_alt <- cache_load("p_mat_meso_alt_slides", quote(
  boot_pmatrix_ci(data_alt, meso_pop, meso_son,
                  df_linked = data_alt, df_full = aian_full, R = R_slides, .seed = 123)
))

# ── Stationary distributions ───────────────────────────────────────────────────
pi0_macro     <- recode_occ_vec(pi_0(data,     macro_pop))
pi0_meso      <- recode_occ_vec(pi_0(data,     meso_pop))
pi0_macro_alt <- recode_occ_vec(pi_0(data_alt, macro_pop))
pi0_meso_alt  <- recode_occ_vec(pi_0(data_alt, meso_pop))

steady_macro     <- recode_occ_vec(pi_star(p_matrix(data,     macro_pop, macro_son, TRUE)))
steady_meso      <- recode_occ_vec(pi_star(p_matrix(data,     meso_pop,  meso_son,  TRUE)))
steady_macro_alt <- recode_occ_vec(pi_star(p_matrix(data_alt, macro_pop, macro_son, TRUE)))
steady_meso_alt  <- recode_occ_vec(pi_star(p_matrix(data_alt, meso_pop,  meso_son,  TRUE)))

# ── Recode bootstrap data frames ───────────────────────────────────────────────
p_mat_macro_r     <- recode_occ_df(p_mat_macro,     macro_pop, macro_son)
p_mat_meso_r      <- recode_occ_df(p_mat_meso,      meso_pop,  meso_son)
p_mat_macro_alt_r <- recode_occ_df(p_mat_macro_alt, macro_pop, macro_son)
p_mat_meso_alt_r  <- recode_occ_df(p_mat_meso_alt,  meso_pop,  meso_son)

# ── Assemble combined plots ────────────────────────────────────────────────────
make_combined <- function(boot_df, pi0, steady, dad_var, son_var,
                          levels = NULL, text_size = 4) {
  g      <- plot_pmat_heatmap(boot_df, !!rlang::ensym(dad_var),
                              !!rlang::ensym(son_var),
                              levels = levels, text_size = text_size,
                              title_str = "P")
  g0     <- plot_pi_column(pi0,    title_str = "π₀",  levels = levels)
  g_star <- plot_pi_column(steady, title_str = "π*",       levels = levels)
  g + g0 + g_star + plot_layout(widths = c(6, 1, 1))
}

combined_macro     <- make_combined(p_mat_macro_r,     pi0_macro,     steady_macro,
                                    macro_pop, macro_son,
                                    levels = macro_level_order)
combined_meso      <- make_combined(p_mat_meso_r,      pi0_meso,      steady_meso,
                                    meso_pop,  meso_son,
                                    levels = meso_level_order, text_size = 4)
combined_macro_alt <- make_combined(p_mat_macro_alt_r, pi0_macro_alt, steady_macro_alt,
                                    macro_pop, macro_son,
                                    levels = macro_level_order)
combined_meso_alt  <- make_combined(p_mat_meso_alt_r,  pi0_meso_alt,  steady_meso_alt,
                                    meso_pop,  meso_son,
                                    levels = meso_level_order, text_size = 4)

# ── EM/SM mobility curves: macro vs macro-alt ─────────────────────────────────
message("Running macro EM/SM bootstrap …")
macro_om <- cache_load("macro_om_slides", quote(
  mobility_curve_with_boot(data, macro_pop, macro_son,
                           df_linked = data, df_full = aian_full,
                           ts = 0:4, R = 100, .seed = 123)
))

message("Running macro-alt EM/SM bootstrap …")
macro_alt_om <- cache_load("macro_alt_om_slides", quote(
  mobility_curve_with_boot(data_alt, macro_pop, macro_son,
                           df_linked = data_alt, df_full = aian_full,
                           ts = 0:4, R = 100, .seed = 123)
))

om_combined <- bind_rows(
  macro_om     |> filter(measure != "SM") |> mutate(series = "Main"),
  macro_alt_om |> filter(measure != "SM") |> mutate(series = "Alt")
) |> mutate(
  lo = est - 1.96 * se,
  hi = est + 1.96 * se
)

om_plot <- ggplot(om_combined, aes(x = t, y = est)) +
  geom_ribbon(data = filter(om_combined, measure == "EM"),
              aes(ymin = lo, ymax = hi, fill = series,
                  group = interaction(series, measure)),
              alpha = 0.2, color = NA) +
  geom_line(aes(color = series, linetype = measure,
                group = interaction(series, measure)), linewidth = 1) +
  geom_point(aes(shape = measure, color = series), size = 2.5) +
  scale_x_continuous(breaks = 0:4) +
  scale_y_continuous(breaks = seq(0, 0.7, 0.1)) +
  scale_linetype_manual(name = "Measure",
                        labels = c(EM = "Exchange Mobility", OM = "Overall Mobility"),
                        values = c(EM = "dashed", OM = "solid")) +
  scale_shape_manual(name = "Measure",
                     labels = c(EM = "Exchange Mobility", OM = "Overall Mobility"),
                     values = c(EM = 17, OM = 16)) +
  scale_fill_manual(name = "Sample",
                    values = c("Main" = "#D55E00",
                               "Alt" = "#0072B2")) +
  scale_color_manual(name = "Sample",
                     values = c("Main" = "#D55E00",
                                "Alt" = "#0072B2")) +
  labs(x = "Generation (t)", y = "Probability to move") +
  theme_minimal() +
  theme(legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.text       = element_text(size = 12),
        legend.title      = element_text(size = 13, face = "bold"),
        axis.line.x       = element_line(linewidth = 0.5),
        axis.line.y       = element_line(linewidth = 0.5),
        axis.text         = element_text(size = 12),
        axis.title        = element_text(size = 13),
        axis.ticks        = element_line(linewidth = 0.5),
        panel.grid.minor  = element_blank()) +
  coord_cartesian(xlim = c(0, 4), ylim = c(0, 0.7))

# ── Save ───────────────────────────────────────────────────────────────────────
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

ggsave("output/figures/overall_mobility.png", om_plot,            width = 10, height = 7, dpi = 200)
ggsave("output/figures/trans_macro.png",     combined_macro,     width = 10, height = 7, dpi = 200)
ggsave("output/figures/trans_meso.png",      combined_meso,      width = 10, height = 7, dpi = 200)
ggsave("output/figures/trans_macro_alt.png", combined_macro_alt, width = 10, height = 7, dpi = 200)
ggsave("output/figures/trans_meso_alt.png",  combined_meso_alt,  width = 10, height = 7, dpi = 200)

message("Done. Figures saved to output/figures/")
