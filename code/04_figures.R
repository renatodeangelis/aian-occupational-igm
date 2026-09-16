################################################################################
# 04_figures.R
# Publication figures from estimated transition matrices.
# Layout: pi_0 | P | pi* (widths 1:6:1).
# Dobrushin worst-row pair highlighted with black border.
#
# Reads:  output/estimates.rds
# Writes: output/presentation/fig1_macro_matrix.{pdf,png}
#         output/presentation/fig2_meso_matrix.{pdf,png}
#         output/presentation/region_map_check.png
################################################################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(sf)
library(maps)

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
macro_levels = macro_compute_order
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

ggsave("output/presentation/fig1_macro_matrix.pdf", fig_macro,
       width = 10, height = 7, dpi = 300)
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

ggsave("output/presentation/fig2_meso_matrix.pdf", fig_meso,
       width = 12, height = 8, dpi = 300)
ggsave("output/presentation/fig2_meso_matrix.png", fig_meso,
       width = 12, height = 8, dpi = 300)

message("Wrote fig2_meso_matrix")

################################################################################
# FIGURE 3: REGION REFERENCE MAP
# Categorical fills; n per region as two-line label.
# maps::map() IDs include ":suffix" variants (e.g. "michigan:north") — strip
# before joining to state_fips_1940.
################################################################################

region_order = c("sw", "south", "cali", "ok", "plains", "nw", "north")

region_labels = c(
  sw     = "Southwest",
  south  = "South",
  cali   = "California",
  ok     = "Oklahoma",
  plains = "Plains",
  nw     = "Northwest",
  north  = "North"
)

# Count n per region from the global weighted data
data_global = readRDS("data/aian_weighted.rds")

region_n = data_global |>
  count(region, name = "n") |>
  filter(!is.na(region)) |>
  mutate(
    region      = factor(region, levels = region_order),
    region_name = region_labels[as.character(region)],
    map_label   = sprintf("%s\nn = %s", region_name, format(n, big.mark = ","))
  ) |>
  arrange(region)

sf::sf_use_s2(FALSE)

states_sf = sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) |>
  mutate(state_name = sub(":.*$", "", ID)) |>
  left_join(state_fips_1940, by = "state_name") |>
  filter(!is.na(statefip)) |>
  sf::st_make_valid()

stopifnot(!any(is.na(states_sf$region)))

# Albers Equal Area (5070) for polygon union; reproject to WGS84 (4326) for labels
regions_sf = states_sf |>
  sf::st_transform(5070) |>
  group_by(region) |>
  summarize(geometry = sf::st_union(geom), .groups = "drop") |>
  sf::st_transform(4326) |>
  mutate(region = factor(region, levels = region_order)) |>
  left_join(region_n, by = "region")

# Surface points in Albers (accurate), coordinates extracted in WGS84
label_pts_base = regions_sf |>
  sf::st_transform(5070) |>
  sf::st_point_on_surface() |>
  sf::st_transform(4326) |>
  sf::st_coordinates() |>
  as_tibble() |>
  bind_cols(sf::st_drop_geometry(regions_sf)) |>
  rename(x = X, y = Y)

# ── LABEL NUDGES ── edit dx/dy here and re-run this block ──────────────────
nudge_dx = c(cali = -0.5, ok = 0.0, north =  7.0,
             nw   =  0.0, plains =  0.0, south =  2.5, sw = 0.0)
nudge_dy = c(cali = -1.0, ok = 0.0, north = -3.0,
             nw   = -2.0, plains = -2.5, south =  0.5, sw = 0.0)
# ────────────────────────────────────────────────────────────────────────────

nudge = tibble(
  region = factor(names(nudge_dx), levels = region_order),
  dx     = nudge_dx,
  dy     = nudge_dy
)

label_pts = label_pts_base |>
  left_join(nudge, by = "region") |>
  mutate(xlab = x + dx, ylab = y + dy)

region_fills = c(
  sw     = "#EADBC8",
  south  = "#DCE4D2",
  cali   = "#D9DEE8",
  ok     = "#EFE2DA",
  plains = "#E3E0D5",
  nw     = "#D6E0DE",
  north  = "#E6DCE4"
)

region_map = ggplot() +
  geom_sf(data = regions_sf, aes(fill = region),
          color = "grey35", linewidth = 0.25) +
  geom_text(data = label_pts,
            aes(x = xlab, y = ylab, label = map_label),
            size = 3.5, lineheight = 0.95, color = "grey15",
            fontface = "bold") +
  scale_fill_manual(values = region_fills, guide = "none") +
  coord_sf(crs = sf::st_crs(4326), datum = NA, expand = TRUE) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 2))

ggsave("output/presentation/region_map_check.png", region_map,
       width = 6.8, height = 4.4, units = "in", dpi = 200)

message("Done. Figures written to output/presentation/")
