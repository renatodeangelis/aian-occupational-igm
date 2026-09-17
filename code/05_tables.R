################################################################################
# 06_tables.R
# LaTeX tables for slides 8 and meso summary.
#
# §1  Regional table (shaded LaTeX tabular)
# §2  Meso summary table
# §3  Zero-cell check and Blume-style ε perturbation
#
# Reads:  output/estimates.rds
# Writes: output/figures/slide08_regional_table.tex
#         output/figures/slide08_zero_cell_check.tex
#         (meso table printed to console for copy-paste into .tex)
################################################################################

library(dplyr)
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

# South is excluded from regional comparison (OCC1950 100 conflates owner and
# tenant farmers; tenancy dominated in the South).
regions_cf = setdiff(regions_list, "south")

tbl = do.call(rbind, lapply(regions_cf, function(r) {
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
ln = c(ln, "\\begin{tabular}{lrrrrrrrrr}")
ln = c(ln, "\\toprule")
ln = c(ln, paste(
  "Region", "$n$",
  "$\\pi_0^{\\text{farm}}$",
  "Farm ret.",
  "Farmwkr.",
  "Relief$^-$",
  "Relief$^+$",
  "$\\delta(P)$",
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
    cell_sh(sh_rel_lo[i], tbl$relief_lower[i]),
    cell_sh(sh_rel_hi[i], tbl$relief_upper[i]),
    cell(tbl$d1[i]),
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
