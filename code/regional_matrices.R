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
meso_levels  = setdiff(meso_order,  "nonemp")

################################################################################
# PLOT HELPERS
################################################################################

plot_pmat_heatmap_simple = function(P_df, dad_var, son_var,
                                    levels = NULL, text_size = 4,
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
    scale_fill_gradient(low = "lightyellow", high = "firebrick") +
    labs(x = "Son's occupation", y = "Father's occupation",
         fill = "Transition Prob.", title = title_expr) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 10),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
}

plot_pi_column_simple = function(vec, title_expr, levels = NULL) {
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

################################################################################
# PER-REGION COMPUTATION AND PLOTTING
################################################################################

macro_level_order = c("nonmanual", "manual", "farming")
meso_level_order  = c("nonmanual", "crafts", "unskilled", "farmworker", "farmer")

regions = sort(na.omit(unique(data$region)))

cat(sprintf("Regions found: %s\n", paste(regions, collapse = ", ")))

dir.create("output/figures/regional", recursive = TRUE, showWarnings = FALSE)

for (reg in regions) {
  df_reg = data |>
    filter(region == reg) |>
    mutate(w_atc_norm = w_atc_norm / sum(w_atc_norm) * n())

  cat(sprintf("\n--- %s (n = %d) ---\n", reg, nrow(df_reg)))

  # --- Macro ---
  P_mat  = p_matrix(df_reg, macro_pop, macro_son, matrix = TRUE)
  pi0    = pi_0(df_reg, macro_pop)
  pistar = pi_star(P_mat)

  cat("P (macro):\n"); print(round(P_mat, 3))
  cat("pi_0 (macro):  "); print(round(pi0, 3))
  cat("pi*  (macro):  "); print(round(pistar, 3))

  P_df = tibble(
    macro_pop = rep(rownames(P_mat), times = ncol(P_mat)),
    macro_son = rep(colnames(P_mat), each  = nrow(P_mat)),
    est       = as.vector(P_mat)
  )

  g_P      = plot_pmat_heatmap_simple(P_df, macro_pop, macro_son,
                                      levels = macro_level_order,
                                      title_expr = paste0("Macro  —  ", reg,
                                                          "  (n = ", nrow(df_reg), ")"))
  g_pi0    = plot_pi_column_simple(pi0,    expression(pi[0]), levels = macro_level_order)
  g_pistar = plot_pi_column_simple(pistar, expression(pi^"*"), levels = macro_level_order)

  combined_macro = g_P + g_pi0 + g_pistar + plot_layout(widths = c(5, 1, 1))

  # --- Meso ---
  P_meso      = p_matrix(df_reg, meso_pop, meso_son, matrix = TRUE)
  pi0_meso    = pi_0(df_reg, meso_pop)
  pistar_meso = pi_star(P_meso)

  cat("P (meso):\n"); print(round(P_meso, 3))
  cat("pi_0 (meso):  "); print(round(pi0_meso, 3))
  cat("pi*  (meso):  "); print(round(pistar_meso, 3))

  P_meso_df = tibble(
    meso_pop = rep(rownames(P_meso), times = ncol(P_meso)),
    meso_son = rep(colnames(P_meso), each  = nrow(P_meso)),
    est      = as.vector(P_meso)
  )

  g_P_m      = plot_pmat_heatmap_simple(P_meso_df, meso_pop, meso_son,
                                        levels = meso_level_order, text_size = 3,
                                        title_expr = paste0("Meso  —  ", reg,
                                                            "  (n = ", nrow(df_reg), ")"))
  g_pi0_m    = plot_pi_column_simple(pi0_meso,    expression(pi[0]), levels = meso_level_order)
  g_pistar_m = plot_pi_column_simple(pistar_meso, expression(pi^"*"), levels = meso_level_order)

  combined_meso = g_P_m + g_pi0_m + g_pistar_m + plot_layout(widths = c(5, 1, 1))

  combined = combined_macro / combined_meso

  ggsave(
    filename = file.path("output/figures/regional", paste0("pmat_", reg, ".png")),
    plot     = combined,
    width = 10, height = 12, dpi = 300
  )
}

cat("\nDone. PNGs saved to output/figures/regional/\n")
