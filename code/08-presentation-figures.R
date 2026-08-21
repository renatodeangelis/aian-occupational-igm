################################################################################
# Presentation figures: macro- and meso-level transition matrices
# Layout: pi_0 | P matrix | pi*
# Bootstrap: R = 50 (SE only, no caching)
# Output: output/presentation/ as PDF + PNG (300 dpi)
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(expm)
library(sf)
library(maps)

source("code/utils.R")

data = readRDS("data/aian_weighted.rds") |>
  mutate(w_atc_norm = w_trim_norm)

aian_full = readRDS("data/aian_full.rds")

# pi_0() looks for these in the parent frame
macro_levels = unique(data$macro_pop)
meso_levels  = unique(data$meso_pop)

################################################################################
# PLOT HELPERS
################################################################################

# Presentation heatmap for the P matrix
plot_pmat_pres = function(boot_df, dad_var, son_var,
                          levels = NULL, text_size = 5.5, title_expr = "P") {
  dad_sym = ensym(dad_var)
  son_sym = ensym(son_var)

  plot_df = boot_df
  if (!is.null(levels)) {
    plot_df = plot_df |>
      mutate(!!dad_sym := factor(!!dad_sym, levels = levels),
             !!son_sym := factor(!!son_sym, levels = rev(levels)))
  }

  ggplot(plot_df, aes(x = !!son_sym, y = !!dad_sym, fill = est)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f\n(%.3f)", est, se)),
              vjust = 0.3, size = text_size) +
    scale_fill_gradient(low = "lightyellow", high = "firebrick",
                        limits = c(0, 1)) +
    guides(fill = guide_colorbar(barwidth = unit(7, "cm"), barheight = unit(0.5, "cm"))) +
    labs(x = "Son's occupation", y = NULL,
         fill = "Prob.", title = title_expr) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 13),
      axis.text.y     = element_text(angle = 45, hjust = 1, size = 13),
      axis.ticks      = element_blank(),
      axis.title.x    = element_text(size = 14),
      legend.position = "bottom",
      legend.text     = element_text(size = 12),
      legend.title    = element_text(size = 13),
      plot.title      = element_text(hjust = 0.5, size = 20, face = "bold"),
      panel.grid      = element_blank()
    )
}

# Presentation single-column tile for pi vectors
plot_pi_pres = function(vec, title_expr, levels = NULL) {
  df = tibble(occ = names(vec), value = as.numeric(vec))
  if (!is.null(levels)) {
    df = df |> mutate(occ = factor(occ, levels = levels))
  } else {
    df = df |> mutate(occ = factor(occ, levels = rev(unique(occ))))
  }

  ggplot(df, aes(x = 1, y = occ, fill = value)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", value)),
              size = 5) +
    scale_fill_gradient(low = "lightyellow", high = "firebrick",
                        limits = c(0, 1)) +
    labs(y = "Father's occupation", title = title_expr) +
    theme_minimal(base_size = 16) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.text.y  = element_text(angle = 45, hjust = 1, size = 13),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.title   = element_text(hjust = 0.5, size = 20, face = "bold"),
      panel.grid   = element_blank()
    )
}

################################################################################
# BOOTSTRAP (R = 50)
################################################################################

p_mat_macro = boot_pmatrix_ci(
  data, macro_pop, macro_son,
  df_linked = data, df_full = aian_full,
  R = 50, .seed = 42
)

p_mat_meso = boot_pmatrix_ci(
  data, meso_pop, meso_son,
  df_linked = data, df_full = aian_full,
  R = 50, .seed = 42
)

################################################################################
# INITIAL AND STATIONARY DISTRIBUTIONS
################################################################################

pi0_macro = pi_0(data, macro_pop)
pi0_meso  = pi_0(data, meso_pop)

P_macro      = p_matrix(data, macro_pop, macro_son, matrix = TRUE)
P_meso       = p_matrix(data, meso_pop,  meso_son,  matrix = TRUE)
steady_macro = pi_star(P_macro)
steady_meso  = pi_star(P_meso)

################################################################################
# FIGURE 1: MACRO TRANSITION MATRIX
################################################################################

g_pmac   = plot_pmat_pres(p_mat_macro, macro_pop, macro_son,
                          levels = macro_order,
                          title_expr = expression(italic(P)))
g0_mac   = plot_pi_pres(pi0_macro, title_expr = expression(pi[0]),
                        levels = macro_order)
gst_mac  = plot_pi_pres(steady_macro, title_expr = expression(pi^"*"),
                        levels = macro_order)

fig_macro = g0_mac + g_pmac + gst_mac +
  plot_layout(widths = c(1, 6, 1))

################################################################################
# FIGURE 2: MESO TRANSITION MATRIX
################################################################################

g_pmes   = plot_pmat_pres(p_mat_meso, meso_pop, meso_son,
                          levels = meso_order,
                          text_size = 4.5,
                          title_expr = expression(italic(P)))
g0_mes   = plot_pi_pres(pi0_meso, title_expr = expression(pi[0]),
                        levels = meso_order)
gst_mes  = plot_pi_pres(steady_meso, title_expr = expression(pi^"*"),
                        levels = meso_order)

fig_meso = g0_mes + g_pmes + gst_mes +
  plot_layout(widths = c(1, 6, 1))

################################################################################
# SAVE
################################################################################

dir.create("output/presentation", recursive = TRUE, showWarnings = FALSE)

ggsave("output/presentation/fig1_macro_matrix.pdf", fig_macro,
       width = 10, height = 7, dpi = 300)
ggsave("output/presentation/fig1_macro_matrix.png", fig_macro,
       width = 10, height = 7, dpi = 300)

ggsave("output/presentation/fig2_meso_matrix.pdf",  fig_meso,
       width = 12, height = 8, dpi = 300)
ggsave("output/presentation/fig2_meso_matrix.png",  fig_meso,
       width = 12, height = 8, dpi = 300)

message("Done. Figures saved to output/presentation/")

################################################################################
# MOBILITY SUMMARY TABLE: GLOBAL + REGIONAL (MACRO LEVEL)
# λ₂, SM(1), EM(1) — no bootstrap
################################################################################

# Compute λ₂, SM(1), EM(1) from raw father-son pairs at the macro level.
# Weights are renormalized within df before P and pi_0 are estimated,
# consistent with the within-region convention in 03_transition_matrices_weighted.R.
compute_macro_stats = function(df, label) {
  df  = mutate(df, w_atc_norm = w_atc_norm / sum(w_atc_norm) * n())
  P   = p_matrix(df, macro_pop, macro_son, matrix = TRUE)
  pi0 = pi_0(df, macro_pop)
  pis = pi_star(P)

  lambda2 = sort(Mod(eigen(P)$values), decreasing = TRUE)[2]

  sm1 = sm(P, pi0, t = 1)
  em1 = om(P, pi0, t = 1) - sm1

  # Safe indexers for categories that may be absent in sparse regional matrices
  cell  = function(r, c) if (r %in% rownames(P) && c %in% colnames(P)) P[r, c] else NA_real_
  pisel = function(k)    if (k %in% names(pis)) pis[k] else NA_real_

  tibble(
    unit            = label,
    n               = nrow(df),
    lambda2         = round(lambda2,                              3),
    sm1             = round(sm1,                                  3),
    em1             = round(em1,                                  3),
    p_farm_farm     = round(cell("farming", "farming"),           3),
    pistar_farming  = round(pisel("farming"),                     3),
    pistar_nonemp   = round(pisel("nonemp"),                      3),
    p_relief = round(weighted.mean(df$empstatd_1940 == 11, df$w_atc_norm, na.rm = TRUE), 3)
  )
}

regional_data = readRDS("data/aian_regional_weighted.rds")

mobility_stats = bind_rows(
  compute_macro_stats(data, "global"),
  lapply(names(regional_data), function(r)
    compute_macro_stats(regional_data[[r]], r))
)

print(mobility_stats, n = Inf)

################################################################################
# PRINT MATRICES AND DISTRIBUTIONS
################################################################################

print_matrix_block = function(label, P, pi0, pistar) {
  sep = strrep("-", nchar(label) + 4)
  cat(sprintf("\n%s\n  %s\n%s\n", sep, label, sep))
  cat("P:\n");    print(round(P,      3))
  cat("pi_0:\n"); print(round(pi0,    3))
  cat("pi*:\n");  print(round(pistar, 3))
}

print_matrix_block("GLOBAL MACRO", P_macro, pi0_macro, steady_macro)
print_matrix_block("GLOBAL MESO",  P_meso,  pi0_meso,  steady_meso)

for (r in names(regional_data)) {
  df    = mutate(regional_data[[r]], w_atc_norm = w_atc_norm / sum(w_atc_norm) * n())
  P_r   = p_matrix(df, macro_pop, macro_son, matrix = TRUE)
  pi0_r = pi_0(df, macro_pop)
  pis_r = pi_star(P_r)
  print_matrix_block(paste("REGION:", toupper(r)), P_r, pi0_r, pis_r)
}

################################################################################
# FIGURE: SEVEN-REGION REFERENCE MAP (slide 4 inset)
# Categorical fills (identity map); n per region as two-line label.
################################################################################

# Ordered by pi*_farming descending, matching the slide 8 table.
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

region_n = data |>
  count(region, name = "n") |>
  filter(!is.na(region)) |>
  mutate(
    region      = factor(region, levels = region_order),
    region_name = region_labels[as.character(region)],
    map_label   = sprintf("%s\nn = %s", region_name, format(n, big.mark = ","))
  ) |>
  arrange(region)

sf::sf_use_s2(FALSE)

state_fips = tibble(
  state_name = tolower(c(
    "alabama","arizona","arkansas","california","colorado","connecticut","delaware",
    "florida","georgia","idaho","illinois","indiana","iowa","kansas","kentucky",
    "louisiana","maine","maryland","massachusetts","michigan","minnesota",
    "mississippi","missouri","montana","nebraska","nevada","new hampshire",
    "new jersey","new mexico","new york","north carolina","north dakota","ohio",
    "oklahoma","oregon","pennsylvania","rhode island","south carolina","south dakota",
    "tennessee","texas","utah","vermont","virginia","washington","west virginia",
    "wisconsin","wyoming")),
  statefip = c(
     1, 4, 5, 6, 8, 9,10,
    12,13,16,17,18,19,20,21,
    22,23,24,25,26,27,28,29,
    30,31,32,33,34,35,36,37,
    38,39,40,41,42,44,45,46,
    47,48,49,50,51,53,54,55,56)
) |>
  mutate(region = assign_region(statefip))

states_sf = st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) |>
  # maps::map() IDs include suffixes like "michigan:north"; strip before joining
  mutate(state_name = sub(":.*$", "", ID)) |>
  left_join(state_fips, by = "state_name") |>
  # DC and any other non-state polygons returned by maps::map() have no statefip
  # match and must be dropped before region assignment
  filter(!is.na(statefip)) |>
  st_make_valid()

stopifnot(!any(is.na(states_sf$region)))

# Project to Albers Equal Area for accurate polygon union and surface-point
# placement, then reproject to WGS84 for degree-based label coordinates
regions_sf = states_sf |>
  st_transform(5070) |>
  group_by(region) |>
  summarize(geometry = st_union(geom), .groups = "drop") |>
  st_transform(4326) |>
  mutate(region = factor(region, levels = region_order)) |>
  left_join(region_n, by = "region")

# Surface points computed in Albers (accurate), coordinates extracted in WGS84.
# Kept as label_pts_base so re-running only the nudge block below never causes
# a column collision (dx/dy already in label_pts would become dx.x/dx.y on the
# second left_join, silently breaking the mutate).
label_pts_base = regions_sf |>
  st_transform(5070) |>
  st_point_on_surface() |>
  st_transform(4326) |>
  st_coordinates() |>
  as_tibble() |>
  bind_cols(st_drop_geometry(regions_sf)) |>
  rename(x = X, y = Y)

# ── LABEL NUDGES ── edit these values, then run this entire block ──────────
# dx = degrees east (+) or west (-); dy = degrees north (+) or south (-)
nudge_dx = c(cali = -0.5, ok = 0.0, north =  7.0,
             nw   =  0.0, plains =  0.0, south =  2.5, sw = 0.0)
nudge_dy = c(cali = -1.0, ok = 0.0, north = -3.0,
             nw   =  -2.0, plains = -2.5, south =  0.5, sw = 0.0)
# ──────────────────────────────────────────────────────────────────────────

nudge = tibble(
  region = factor(names(nudge_dx), levels = region_order),
  dx     = nudge_dx,
  dy     = nudge_dy
)

label_pts = label_pts_base |>
  left_join(nudge, by = "region") |>
  mutate(xlab = x + dx, ylab = y + dy,
         leader = (dx != 0 | dy != 0))

print(label_pts[, c("region", "xlab", "ylab")])

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
  coord_sf(crs = st_crs(4326), datum = NA, expand = TRUE) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 2))

ggsave("output/presentation/region_map_check.png", region_map,
       width = 6.8, height = 4.4, units = "in", dpi = 200)

