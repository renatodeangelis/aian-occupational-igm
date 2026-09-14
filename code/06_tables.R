################################################################################
# 06_tables.R
# LaTeX tables for slides 8 and meso summary.
#
# §1  Re-check assertions from 03_estimate.R (audit trail)
# §2  Regional table (shaded LaTeX tabular)
# §3  Meso summary table
# §4  Zero-cell check and Blume-style ε perturbation
#
# Reads:  data/estimates.rds
# Writes: output/figures/slide08_regional_table.tex
#         output/figures/slide08_zero_cell_check.tex
#         (meso table printed to console for copy-paste into .tex)
################################################################################

library(dplyr)
library(knitr)

source("code/00_utils.R")
source("code/expected_values.R")

est = readRDS("data/estimates.rds")

regional_results = est$regional
regions_list     = est$regions_list
region_display   = est$region_display
P_meso           = est$P_meso
pi0_meso         = est$pi0_meso
steady_meso      = est$steady_meso

################################################################################
# §1 RE-CHECK ASSERTIONS
################################################################################

cat("--- Re-checking regional assertions from estimates.rds ---\n")
num_fields = c("farm_ret", "pi_farming", "pi_manual", "pi_nonman", "lambda2", "relief")

for (r in regions_list) {
  exp = EXPECTED$regional[[r]]
  got = regional_results[[r]]
  if (is.null(exp)) next

  for (fld in num_fields) {
    if (is.na(exp[[fld]])) next
    delta = abs(got[[fld]] - exp[[fld]])
    if (delta > TOL)
      stop(sprintf("[%s] %s: got %.4f, expected %.4f (|delta|=%.4f > TOL=%.3f)",
                   r, fld, got[[fld]], exp[[fld]], delta, TOL))
  }
  if (!is.na(exp$n) && got$n != exp$n)
    stop(sprintf("[%s] n: got %d, expected %d", r, got$n, exp$n))
}
cat("All regional assertions passed.\n\n")

################################################################################
# §2 REGIONAL TABLE
# Sorted by pi*_farming descending; shaded cells use \sh{shade}{val} command.
# Caller must provide \newcommand{\sh}[2]{\cellcolor{black!#1}#2} and load
# booktabs + colortbl.
################################################################################

tbl = do.call(rbind, lapply(regions_list, function(r) {
  g = regional_results[[r]]
  data.frame(
    region     = r,
    display    = region_display[r],
    farm_ret   = g$farm_ret,
    pi_farming = g$pi_farming,
    pi_manual  = g$pi_manual,
    pi_nonman  = g$pi_nonman,
    lambda2    = g$lambda2,
    relief     = g$relief,
    n          = g$n,
    stringsAsFactors = FALSE
  )
}))

tbl = tbl[order(-tbl$pi_farming), ]

# Shade integer [2, 22] — column range rescaled linearly
shade_int = function(x) {
  rng = range(x, na.rm = TRUE)
  if (diff(rng) < 1e-10) return(rep(2L, length(x)))
  as.integer(round(2 + (x - rng[1]) / diff(rng) * 20))
}

sh_farm   = shade_int(tbl$pi_farming)
sh_manual = shade_int(tbl$pi_manual)
sh_relief = shade_int(tbl$relief)

# Strip leading zero
fmt = function(x) sub("^0", "", sprintf("%.3f", x))

cell_sh = function(shade, val, bold = FALSE) {
  v = if (bold) sprintf("\\textbf{%s}", fmt(val)) else fmt(val)
  sprintf("\\sh{%d}{%s}", shade, v)
}

cell = function(val, bold = FALSE) {
  if (bold) sprintf("\\textbf{%s}", fmt(val)) else fmt(val)
}

ln = character(0)
ln = c(ln, "\\begin{tabular}{lrrrrrr}")
ln = c(ln, "\\toprule")
ln = c(ln, paste(
  "Region",
  "Farm ret.",
  "$\\pi^*_{\\text{farming}}$",
  "$\\pi^*_{\\text{manual}}$",
  "$\\pi^*_{\\text{nonmanual}}$",
  "$\\lambda_2$",
  "Relief",
  sep = " & "
), "\\\\")
ln = c(ln, "\\midrule")

for (i in seq_len(nrow(tbl))) {
  r        = tbl$region[i]
  is_ok    = (r == "ok")
  is_south = (r == "south")

  row_str = paste(
    tbl$display[i],
    cell(tbl$farm_ret[i]),
    cell_sh(sh_farm[i],   tbl$pi_farming[i]),
    cell_sh(sh_manual[i], tbl$pi_manual[i]),
    cell(tbl$pi_nonman[i], bold = is_ok),
    cell(tbl$lambda2[i],   bold = is_south),
    cell_sh(sh_relief[i], tbl$relief[i]),
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
ln = c(ln, sprintf("{\\footnotesize %s}", paste(n_entries, collapse = ";\\enspace ")))

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
writeLines(ln, "output/figures/slide08_regional_table.tex")
cat("Wrote output/figures/slide08_regional_table.tex\n")

################################################################################
# §3 MESO SUMMARY TABLE
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
# §4 ZERO-CELL CHECK AND BLUME-STYLE EPSILON PERTURBATION
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
