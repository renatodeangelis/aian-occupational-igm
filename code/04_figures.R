################################################################################
# 04_figures.R
# Publication figures from estimated transition matrices.
# Layout: pi_0 | P | pi* (widths 1:6:1).
# Dobrushin worst-row pair highlighted with black border.
#
# Reads:  output/estimates.rds
# Writes: output/presentation/fig1_macro_matrix.png
#         output/presentation/fig2_meso_matrix.png
################################################################################

library(dplyr)
library(ggplot2)
library(patchwork)

source("code/00_utils.R")

est = readRDS("output/estimates.rds")

# Restore objects from estimates list
p_mat_macro  = est$p_mat_macro
p_mat_meso   = est$p_mat_meso
P_macro      = est$P_macro
P_meso       = est$P_meso
pi0_macro    = est$pi0_macro
pi0_meso     = est$pi0_meso
steady_macro = est$steady_macro
steady_meso  = est$steady_meso
dob_mac      = est$dob_mac
dob_mes      = est$dob_mes

# Rebuild macro_levels / meso_levels for pi_0() if called from this script
macro_levels = macro_order
meso_levels  = meso_order

dir.create("output/presentation", recursive = TRUE, showWarnings = FALSE)

################################################################################
# FIGURE 1: MACRO TRANSITION MATRIX (pi_0 | P | pi*)
################################################################################

g_pmac = plot_pmat(p_mat_macro$P, macro_pop, macro_son,
                   levels = macro_order,
                   title_expr = expression(italic(P))) +
  geom_tile(
    data = filter(p_mat_macro$P, macro_pop %in% c(dob_mac$row1, dob_mac$row2)),
    aes(x = macro_son, y = macro_pop),
    fill = NA, color = "black", linewidth = 0.4,
    inherit.aes = FALSE
  )

g0_mac  = plot_pi(pi0_macro,    title_expr = expression(pi[0]),  levels = macro_order)
gst_mac = plot_pi(steady_macro, title_expr = expression(pi^"*"), levels = macro_order)

fig_macro = g0_mac + g_pmac + gst_mac + plot_layout(widths = c(1, 6, 1))

ggsave("output/presentation/fig1_macro_matrix.png", fig_macro,
       width = 10, height = 7, dpi = 300)

message("Wrote fig1_macro_matrix")

################################################################################
# FIGURE 2: MESO TRANSITION MATRIX (pi_0 | P | pi*)
################################################################################

g_pmes = plot_pmat(p_mat_meso$P, meso_pop, meso_son,
                   levels = meso_order, text_size = 4.5,
                   title_expr = expression(italic(P))) +
  geom_tile(
    data = filter(p_mat_meso$P, meso_pop %in% c(dob_mes$row1, dob_mes$row2)),
    aes(x = meso_son, y = meso_pop),
    fill = NA, color = "black", linewidth = 0.4,
    inherit.aes = FALSE
  )

g0_mes  = plot_pi(pi0_meso,    title_expr = expression(pi[0]),  levels = meso_order)
gst_mes = plot_pi(steady_meso, title_expr = expression(pi^"*"), levels = meso_order)

fig_meso = g0_mes + g_pmes + gst_mes + plot_layout(widths = c(1, 6, 1))

ggsave("output/presentation/fig2_meso_matrix.png", fig_meso,
       width = 12, height = 8, dpi = 300)

message("Wrote fig2_meso_matrix")