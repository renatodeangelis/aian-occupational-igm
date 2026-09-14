################################################################################
# 05_panels.R
# Slide 7 panels — farming-exit bar charts.
# Panel 1 (exit_left):  farming fathers → macro_son destinations (horizontal)
# Panel 2 (exit_right): farming→manual sons → empstatd_1940 (vertical)
# Weight: w_trim_norm throughout (present in data after load_global()).
#
# Reads:  data/aian_weighted.rds  (via load_global)
# Writes: output/figures/exit_left.{png,pdf}
#         output/figures/exit_right.{png,pdf}
################################################################################

library(dplyr)
library(ggplot2)
library(scales)

source("code/00_utils.R")
source("code/expected_values.R")

data = load_global()   # sets macro_levels / meso_levels; adds w_atc_norm = w_trim_norm

# Canonical orders — both panels read from here
macro_exit_order = c("farming", "manual", "nonemp", "nonmanual")
emp_order        = c("10", "11", "21", "other")
residual_codes   = c(12, 13, 31, 32, 33, 34)

################################################################################
# PANEL 1 — farming fathers → macro_son
################################################################################

p1_raw = data |>
  filter(macro_pop == "farming") |>
  group_by(macro_son) |>
  summarise(w = sum(w_trim_norm), .groups = "drop") |>
  mutate(prop = w / sum(w))

p1_props = setNames(
  round(p1_raw$prop[match(macro_exit_order, p1_raw$macro_son)], 3),
  macro_exit_order
)

# Assert against expected values (TOL applied to proportions)
if (!is.null(EXPECTED$panel1) && !any(is.na(EXPECTED$panel1))) {
  deltas = abs(p1_props[names(EXPECTED$panel1)] - EXPECTED$panel1)
  if (any(deltas > TOL))
    stop(sprintf("Panel 1 mismatch: max |delta| = %.4f > TOL = %.3f\n  got: %s\n  exp: %s",
                 max(deltas),
                 TOL,
                 paste(round(p1_props, 3), collapse = " / "),
                 paste(EXPECTED$panel1,   collapse = " / ")))
}
cat("Panel 1 proportions:", paste(names(p1_props), round(p1_props, 3), sep = "=", collapse = "  "), "\n")

p1_data = tibble(
  cat    = factor(macro_exit_order, levels = macro_exit_order),
  prop   = p1_props,
  xend   = cumsum(prop),
  xstart = lag(xend, default = 0),
  xmid   = (xstart + xend) / 2
)

p1_fills = c(
  farming   = "#282828",
  manual    = "#696969",
  nonemp    = "#a8a8a8",
  nonmanual = "#d6d6d6"
)

panel1 = ggplot(p1_data) +
  geom_rect(aes(xmin = xstart, xmax = xend, ymin = 0, ymax = 1, fill = cat),
            colour = NA) +
  geom_text(
    data  = filter(p1_data, cat != "nonmanual"),
    aes(x = xmid, y = 0.5, label = sprintf("%.3f", prop)),
    colour = "white", size = 3.8, fontface = "bold"
  ) +
  geom_segment(
    data = filter(p1_data, cat == "nonmanual"),
    aes(x = xend + 0.004, xend = 1.048, y = 0.5, yend = 0.5),
    linewidth = 0.35, colour = "#282828"
  ) +
  geom_text(
    data  = filter(p1_data, cat == "nonmanual"),
    aes(x = 1.053, y = 0.5, label = sprintf("%.3f", prop)),
    colour = "#282828", size = 3.8, fontface = "bold", hjust = 0
  ) +
  scale_fill_manual(values = p1_fills) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     oob = scales::oob_keep) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    legend.position  = "none",
    plot.margin      = margin(0, 42, 0, 0, "pt"),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA)
  )

################################################################################
# PANEL 2 — farming→manual sons → empstatd_1940
################################################################################

p2_raw = data |>
  filter(macro_pop == "farming", macro_son == "manual") |>
  mutate(emp_grp = case_when(
    empstatd_1940 == 10                ~ "10",
    empstatd_1940 == 11                ~ "11",
    empstatd_1940 == 21                ~ "21",
    empstatd_1940 %in% residual_codes  ~ "other",
    TRUE                               ~ NA_character_
  )) |>
  filter(!is.na(emp_grp)) |>
  group_by(emp_grp) |>
  summarise(w = sum(w_trim_norm), .groups = "drop") |>
  mutate(prop = w / sum(w))

p2_props = setNames(
  round(p2_raw$prop[match(emp_order, p2_raw$emp_grp)], 3),
  emp_order
)

if (!is.null(EXPECTED$panel2) && !any(is.na(EXPECTED$panel2))) {
  deltas = abs(p2_props[names(EXPECTED$panel2)] - EXPECTED$panel2)
  if (any(deltas > TOL))
    stop(sprintf("Panel 2 mismatch: max |delta| = %.4f > TOL = %.3f\n  got: %s\n  exp: %s",
                 max(deltas),
                 TOL,
                 paste(round(p2_props, 3), collapse = " / "),
                 paste(EXPECTED$panel2,   collapse = " / ")))
}
cat("Panel 2 proportions:", paste(names(p2_props), round(p2_props, 3), sep = "=", collapse = "  "), "\n")

p2_data = tibble(
  cat    = factor(emp_order, levels = emp_order),
  prop   = p2_props,
  yend   = cumsum(prop),
  ystart = lag(yend, default = 0),
  ymid   = (ystart + yend) / 2,
  lcol   = if_else(cat == "other", "#282828", "white")
)

p2_fills = c(
  "10"    = "#696969",
  "11"    = "#9E2B25",
  "21"    = "#a8a8a8",
  "other" = "#d6d6d6"
)

panel2 = ggplot(p2_data) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = ystart, ymax = yend, fill = cat),
            colour = NA) +
  geom_text(
    aes(x = 0.5, y = ymid, label = sprintf("%.3f", prop), colour = lcol),
    size = 3.8, fontface = "bold"
  ) +
  scale_fill_manual(values = p2_fills) +
  scale_colour_identity() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position  = "none",
    plot.margin      = margin(0, 0, 0, 0),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA)
  )

################################################################################
# SAVE
################################################################################

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

ggsave("output/figures/exit_left.png",  panel1,
       width = 5.0, height = 1.5, dpi = 300, bg = "transparent")
ggsave("output/figures/exit_left.pdf",  panel1,
       width = 5.0, height = 1.5, bg = "transparent")

ggsave("output/figures/exit_right.png", panel2,
       width = 1.9, height = 3.2, dpi = 300, bg = "transparent")
ggsave("output/figures/exit_right.pdf", panel2,
       width = 1.9, height = 3.2, bg = "transparent")

message("Done. Panels saved to output/figures/")
