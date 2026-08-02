################################################################################
# regional_matrices.R
# Per-region macro transition matrices, pi_0, and pi_star — no bootstrap.
################################################################################

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(expm)

source("code/utils.R")

data = read_csv("data/aian_weighted.csv") |>
  mutate(w_atc_norm = w_trim_norm) |>
  filter(occ_pop <= 970, occ_son <= 970) |>
  mutate(w_atc_norm = w_atc_norm / sum(w_atc_norm) * n())

macro_levels = setdiff(macro_order, "nonemp")  # pi_0() reads this from the caller environment

################################################################################
# PLOT HELPERS
################################################################################

plot_pmat_heatmap_simple = function(P_df, dad_var, son_var,
                                    levels = NULL, text_size = 3.5,
                                    title_expr = "P") {
  dad_sym = ensym(dad_var)
  son_sym = ensym(son_var)

  plot_df = P_df
  if (!is.null(levels)) {
    plot_df = plot_df |>
      mutate(!!dad_sym := factor(!!dad_sym, levels = levels),
             !!son_sym := factor(!!son_sym, levels = rev(levels)))
  }

  ggplot(plot_df, aes(x = !!son_sym, y = !!dad_sym, fill = est)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", est)), size = text_size) +
    scale_fill_gradient(low = "lightyellow", high = "firebrick", limits = c(0, 1)) +
    labs(x = "Son's occupation", y = "Father's occupation",
         title = title_expr) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(angle = 45, hjust = 1, size = 9),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 9),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 11))
}

plot_pi_column_simple = function(vec, title_expr, levels = NULL) {
  df = tibble(occ = names(vec), value = as.numeric(vec))
  if (!is.null(levels)) {
    df = df |> mutate(occ = factor(occ, levels = levels))
  }

  ggplot(df, aes(x = 1, y = occ, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3.5) +
    scale_fill_gradient(low = "lightyellow", high = "firebrick", limits = c(0, 1)) +
    labs(title = title_expr) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 1, size = 9),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 11))
}

################################################################################
# PER-REGION COMPUTATION AND PLOTTING
################################################################################

regions = sort(na.omit(unique(data$region)))

cat(sprintf("Regions found: %s\n", paste(regions, collapse = ", ")))

dir.create("output/figures/regional", recursive = TRUE, showWarnings = FALSE)

for (reg in regions) {
  df_reg = data |>
    filter(region == reg) |>
    mutate(w_atc_norm = w_atc_norm / sum(w_atc_norm) * n())

  cat(sprintf("\n--- %s (n = %d) ---\n", reg, nrow(df_reg)))

  P_mat = p_matrix(df_reg, macro_pop, macro_son, matrix = TRUE)
  lev   = rownames(P_mat)  # original order (sorted); reversed in plot via rev(lev)
  lev_r = rev(lev)         # display order: passed to plot helpers

  pi0    = pi_0(df_reg, macro_pop)[lev]
  pistar = pi_star(P_mat)

  cat("P:\n"); print(round(P_mat, 3))
  cat("pi_0:  "); print(round(pi0, 3))
  cat("pi*:   "); print(round(pistar, 3))

  # Tidy P matrix for ggplot — column-major matches as.vector(P_mat)
  P_df = tibble(
    macro_pop = rep(rownames(P_mat), times = ncol(P_mat)),
    macro_son = rep(colnames(P_mat), each  = nrow(P_mat)),
    est       = as.vector(P_mat)
  )

  g_P      = plot_pmat_heatmap_simple(P_df, macro_pop, macro_son,
                                      levels = lev_r,
                                      title_expr = paste0("P  —  ", reg,
                                                          "  (n = ", nrow(df_reg), ")"))
  g_pi0    = plot_pi_column_simple(pi0,    expression(pi[0]), levels = lev_r)
  g_pistar = plot_pi_column_simple(pistar, expression(pi^"*"), levels = lev_r)

  combined = g_P + g_pi0 + g_pistar + plot_layout(widths = c(5, 1, 1))

  ggsave(
    filename = file.path("output/figures/regional", paste0("pmat_", reg, ".png")),
    plot     = combined,
    width = 10, height = 6, dpi = 300
  )
}

cat("\nDone. PNGs saved to output/figures/regional/\n")
