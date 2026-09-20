################################################################################
# 05_tables.R
# LaTeX tables for slides 8 and meso summary.
#
# §1  Regional table (shaded LaTeX tabular)
# §2  Meso summary table
# §3  Zero-cell check and Blume-style ε perturbation
# §4  Relief cross-tab (farm origin × employment status)
#
# Reads:  output/estimates.rds
# Writes: output/figures/slide08_regional_table.tex
#         output/figures/slide08_zero_cell_check.tex
#         output/figures/relief_xtab_global.tex
#         output/figures/relief_xtab_regional.tex
#         (meso table printed to console for copy-paste into .tex)
################################################################################

library(dplyr)
library(tidyr)
library(knitr)

source("code/00_utils.R")

est = readRDS("output/estimates.rds")

regional_results = est$regional
reg_cf           = est$reg_cf
reg_cf_boot      = est$reg_cf_boot
regions_list     = est$compare_regions
region_display   = est$region_display
P_meso           = est$P_meso
pi0_meso         = est$pi0_meso
steady_meso      = est$steady_meso

################################################################################
# §1 REGIONAL TABLE
# Sorted by pi*_farming descending; shaded cells use \sh{shade}{val} command.
# Caller must provide \newcommand{\sh}[2]{\cellcolor{black!#1}#2} and load
# booktabs + colortbl.
################################################################################

tbl = do.call(rbind, lapply(regions_list, function(r) {
  g  = regional_results[[r]]
  cf = reg_cf[reg_cf$region == r, ]
  cb = reg_cf_boot[reg_cf_boot$region == r, ]
  data.frame(
    region        = r,
    display       = region_display[r],
    n             = g$n,
    pi0_farming   = g$pi0_farming,
    farm_ret      = g$farm_ret,
    farmwkr_share = g$farmwkr_share_son,
    relief_lower  = g$relief_lower,
    relief_upper  = g$relief_upper,
    d1            = g$d1,
    d1_row1       = g$d1_row1,
    d1_row2       = g$d1_row2,
    d1_lo         = cb$d1_lo,
    d1_hi         = cb$d1_hi,
    ess           = g$ess,
    farmwkr_share_pop = g$farmwkr_share_pop,
    regime        = cf$regime,
    regime_lo     = cb$regime_lo,
    regime_hi     = cb$regime_hi,
    comp          = cf$comp,
    comp_lo       = cb$comp_lo,
    comp_hi       = cb$comp_hi,
    stringsAsFactors = FALSE
  )
}))

tbl = tbl[order(-tbl$pi0_farming), ]

# Shade integer [2, 22] — column range rescaled linearly
shade_int = function(x) {
  rng = range(x, na.rm = TRUE)
  if (diff(rng) < 1e-10) return(rep(2L, length(x)))
  as.integer(round(2 + (x - rng[1]) / diff(rng) * 20))
}

sh_pi0    = shade_int(tbl$pi0_farming)
sh_rel_lo = shade_int(tbl$relief_lower)
sh_rel_hi = shade_int(tbl$relief_upper)
sh_regime = shade_int(tbl$regime)
sh_comp   = shade_int(tbl$comp)

# Strip leading zero
fmt = function(x) sub("^0", "", sprintf("%.3f", x))

cell_sh = function(shade, val, bold = FALSE) {
  v = if (bold) sprintf("\\textbf{%s}", fmt(val)) else fmt(val)
  sprintf("\\sh{%d}{%s}", shade, v)
}

cell = function(val, bold = FALSE) {
  if (bold) sprintf("\\textbf{%s}", fmt(val)) else fmt(val)
}

# Interval cell: "est [lo, hi]" in footnotesize
cell_ci = function(est, lo, hi, shade = NULL) {
  ci_str = sprintf("\\footnotesize[%s,\\,%s]", fmt(lo), fmt(hi))
  v = sprintf("%s %s", fmt(est), ci_str)
  if (!is.null(shade)) sprintf("\\sh{%d}{%s}", shade, v) else v
}

ln = character(0)
ln = c(ln, "\\begin{tabular}{lrrrrrrrrrrrr}")
ln = c(ln, "\\toprule")
ln = c(ln, paste(
  "Region", "$n$",
  "$\\pi_0^{\\text{farm}}$",
  "Farm ret.",
  "Fwkr.$^{\\text{son}}$",
  "Fwkr.$^{\\text{dad}}$",
  "Relief$^-$",
  "Relief$^+$",
  "$\\log\\delta(P)$",
  "Gen. pair",
  "regime$_k$",
  "comp$_k$",
  sep = " & "
), "\\\\")
ln = c(ln, "\\midrule")

for (i in seq_len(nrow(tbl))) {
  row_str = paste(
    tbl$display[i],
    format(tbl$n[i], big.mark = ",", trim = TRUE),
    cell_sh(sh_pi0[i],    tbl$pi0_farming[i]),
    cell(tbl$farm_ret[i]),
    cell(tbl$farmwkr_share[i]),
    cell(tbl$farmwkr_share_pop[i]),
    cell_sh(sh_rel_lo[i], tbl$relief_lower[i]),
    cell_sh(sh_rel_hi[i], tbl$relief_upper[i]),
    cell_ci(tbl$d1[i], tbl$d1_lo[i], tbl$d1_hi[i]),
    sprintf("%s vs %s", tbl$d1_row1[i], tbl$d1_row2[i]),
    cell_ci(tbl$regime[i], tbl$regime_lo[i], tbl$regime_hi[i], sh_regime[i]),
    cell_ci(tbl$comp[i],   tbl$comp_lo[i],   tbl$comp_hi[i],   sh_comp[i]),
    sep = " & "
  )
  ln = c(ln, paste0(row_str, " \\\\"))
}

ln = c(ln, "\\bottomrule")
ln = c(ln, "\\end{tabular}")

tbl_by_n  = tbl[order(-tbl$n), ]
n_entries = sprintf("%s: %s",
                    tbl_by_n$display,
                    format(tbl_by_n$n, big.mark = ",", trim = TRUE))
ln = c(ln, "")
ln = c(ln, "% South excluded: OCC1950 100 conflates owner and tenant farmers.")
ln = c(ln, "")
ln = c(ln, sprintf("{\\footnotesize %s}", paste(n_entries, collapse = ";\\enspace ")))

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
writeLines(ln, "output/figures/slide08_regional_table.tex")
cat("Wrote output/figures/slide08_regional_table.tex\n")

################################################################################
# §2 MESO SUMMARY TABLE
################################################################################

meso_tbl = tibble(
  Category              = meso_order,
  `Retention rate`      = diag(P_meso[meso_order, meso_order]),
  `Inflow to unskilled` = P_meso[meso_order, "unskilled"],
  `$\\pi_0$`            = pi0_meso[meso_order],
  `$\\pi^*$`            = steady_meso[meso_order]
) |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) |>
  arrange(desc(`Retention rate`))

cat("\n--- Meso summary table (LaTeX) ---\n")
cat(knitr::kable(meso_tbl, format = "latex", booktabs = TRUE,
                 escape = FALSE, linesep = ""))
cat("\n")

################################################################################
# §3 ZERO-CELL CHECK AND BLUME-STYLE EPSILON PERTURBATION
################################################################################

eps_vals    = c(.001, .005, .01)
zero_thresh = 1e-8

perturb_P = function(P, eps) {
  Pp = P
  for (i in seq_len(nrow(P))) {
    row = P[i, ]
    low = row < eps
    if (!any(low)) next
    added    = sum(eps - row[low])
    high     = !low
    high_sum = sum(row[high])
    row[low] = eps
    if (high_sum > added)
      row[high] = row[high] * (high_sum - added) / high_sum
    Pp[i, ] = row / sum(row)
  }
  Pp
}

cl = character(0)
cl = c(cl, "% Zero-cell check and Blume-style epsilon perturbation")
cl = c(cl, "% eps in {.001, .005, .01}; replace P[i,j] < eps with eps,")
cl = c(cl, "% rescale remaining row mass proportionally, renormalise.")
cl = c(cl, "")

max_d_nonman  = 0
max_d_lambda2 = 0
regions_zero  = character(0)

for (r in regions_list) {
  P      = regional_results[[r]]$P
  n_zero = sum(P < zero_thresh)
  pis0   = pi_star(P)
  lam0   = sort(Mod(eigen(P, only.values = TRUE)$values), decreasing = TRUE)[2]

  if (n_zero > 0) {
    regions_zero = c(regions_zero, region_display[r])
    cl = c(cl, sprintf("%% %s: %d zero cell(s)", region_display[r], n_zero))
    for (eps in eps_vals) {
      Pp   = perturb_P(P, eps)
      pisp = pi_star(Pp)
      lamp = sort(Mod(eigen(Pp, only.values = TRUE)$values), decreasing = TRUE)[2]
      d_nm  = abs(pisp["nonmanual"] - pis0["nonmanual"])
      d_l2  = abs(lamp - lam0)
      max_d_nonman  = max(max_d_nonman,  d_nm,  na.rm = TRUE)
      max_d_lambda2 = max(max_d_lambda2, d_l2,  na.rm = TRUE)
      cl = c(cl, sprintf(
        "%%   eps=%.3f  Delta pi_nonmanual=%.4f  Delta lambda2=%.4f",
        eps, d_nm, d_l2))
    }
  } else {
    cl = c(cl, sprintf("%% %s: no zero cells", region_display[r]))
  }
}

cl = c(cl, "")
cl = c(cl, sprintf(
  "%% Summary: max Delta pi_nonmanual = %.4f; max Delta lambda2 = %.4f",
  max_d_nonman, max_d_lambda2))
cl = c(cl, "")

oxford = function(x) {
  if (length(x) <= 1) return(x)
  paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
}
cl = c(cl, sprintf(paste(
  "Only %s have a zero cell.",
  "Under Blume-style $\\varepsilon$ perturbation ($\\varepsilon \\in \\{.001,.005,.01\\}$),",
  "$\\pi^*_{\\text{nonmanual}}$ shifts by at most %.3f",
  "and $\\lambda_2$ by at most %.3f."
), oxford(regions_zero), max_d_nonman, max_d_lambda2))

writeLines(cl, "output/figures/slide08_zero_cell_check.tex")
cat("Wrote output/figures/slide08_zero_cell_check.tex\n")

################################################################################
# §4 RELIEF CROSS-TAB
# Weighted cross-tab of sons' macro destination × 1940 employment status,
# conditional on farm origin. Row percentages; n column = unweighted farm-origin
# sons in that destination row.
#
# Column order: at_work | emergency | unemployed | other | na | n
# Rows: macro_order (nonemp, nonmanual, manual, farming).
#
# Caller must load booktabs. Use \input{} to embed in paper.
################################################################################

relief_xtab_global   = est$relief_xtab_global
relief_xtab_regional = est$relief_xtab_regional

col_order  = c("at_work", "emergency", "unemployed", "other", "na")
col_labels = c("At work", "Emergency", "Unemployed", "Other", "Missing")
row_order  = c("nonemp", "nonmanual", "manual", "farming")
row_labels = c("Non-employed", "Non-manual", "Manual", "Farming")

fmt_pct = function(x) sprintf("%.3f", x)

format_relief_tex = function(xtab, title_str = NULL) {
  wide = tidyr::pivot_wider(xtab,
                             id_cols    = c(macro_son, n_row),
                             names_from = emp_cat,
                             values_from = row_pct)
  wide = wide[match(row_order, wide$macro_son), ]

  ln = character(0)
  if (!is.null(title_str))
    ln = c(ln, sprintf("%% %s", title_str))
  ln = c(ln, "\\begin{tabular}{lrrrrrrr}")
  ln = c(ln, "\\toprule")
  ln = c(ln, paste(
    c("Destination", col_labels, "$n$"),
    collapse = " & "
  ), "\\\\")
  ln = c(ln, "\\midrule")

  for (i in seq_len(nrow(wide))) {
    dest = row_labels[match(wide$macro_son[i], row_order)]
    vals = sapply(col_order, function(cc) {
      v = wide[[cc]][i]
      if (is.na(v)) "---" else fmt_pct(v)
    })
    n   = format(wide$n_row[i], big.mark = ",", trim = TRUE)
    ln = c(ln, paste(c(dest, vals, n), collapse = " & "), "\\\\")
  }

  ln = c(ln, "\\bottomrule")
  ln = c(ln, "\\end{tabular}")
  ln
}

# Global table
gl = format_relief_tex(relief_xtab_global, "Relief cross-tab — global")
writeLines(gl, "output/figures/relief_xtab_global.tex")
cat("Wrote output/figures/relief_xtab_global.tex\n")

# Regional table — stacked panels separated by \midrule
reg_lines = character(0)
reg_lines = c(reg_lines, "% Relief cross-tab — per region (farm origin only)")
reg_lines = c(reg_lines, "\\begin{tabular}{llrrrrrrr}")
reg_lines = c(reg_lines, "\\toprule")
reg_lines = c(reg_lines, paste(
  c("Region", "Destination", col_labels, "$n$"),
  collapse = " & "
), "\\\\")
reg_lines = c(reg_lines, "\\midrule")

for (ri in seq_along(compare_regions)) {
  r    = compare_regions[ri]
  disp = region_display[r]
  xtab = relief_xtab_regional[[r]]
  wide = tidyr::pivot_wider(xtab,
                             id_cols     = c(macro_son, n_row),
                             names_from  = emp_cat,
                             values_from = row_pct)
  wide = wide[match(row_order, wide$macro_son), ]

  for (i in seq_len(nrow(wide))) {
    dest = row_labels[match(wide$macro_son[i], row_order)]
    reg_label = if (i == 1) disp else ""
    vals = sapply(col_order, function(cc) {
      v = wide[[cc]][i]
      if (is.na(v)) "---" else fmt_pct(v)
    })
    n = format(wide$n_row[i], big.mark = ",", trim = TRUE)
    reg_lines = c(reg_lines,
                  paste(c(reg_label, dest, vals, n), collapse = " & "), "\\\\")
  }
  if (ri < length(compare_regions))
    reg_lines = c(reg_lines, "\\midrule")
}

reg_lines = c(reg_lines, "\\bottomrule")
reg_lines = c(reg_lines, "\\end{tabular}")

writeLines(reg_lines, "output/figures/relief_xtab_regional.tex")
cat("Wrote output/figures/relief_xtab_regional.tex\n")
