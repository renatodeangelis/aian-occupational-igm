################################################################################
# 08_regional_table.R
# Regional macro mobility table for slide 8.
# Replaces regional_measures.tex (which did not reproduce from committed data).
#
# Output:
#   output/figures/slide08_regional_table.tex
#   output/figures/slide08_zero_cell_check.tex
#
# Caller supplies \newcommand{\sh}[2]{\cellcolor{black!#1}#2} and loads
# booktabs + colortbl; this file emits only the tabular environment.
################################################################################

library(dplyr)

source("code/utils.R")

# Computation-order override.  utils.R sets macro_order for plot display
# (nonemp, nonmanual, manual, farming).  We redefine here so every named
# access in this script is consistent.  p_matrix() sorts alphabetically —
# farming / manual / nonmanual / nonemp — which matches this order exactly,
# so named indexing is always correct and positional indexing is never used.
macro_order = c("farming", "manual", "nonmanual", "nonemp")

data = readRDS("data/aian_weighted.rds") |>
  mutate(w_atc_norm = w_trim_norm)

################################################################################
# NATIONAL CHECK
# Confirm pi_0 = .659 / .233 / .046 / .062 with w_trim_norm before proceeding.
################################################################################

macro_levels = macro_order          # pi_0() reads this from parent.frame()

pi0_nat  = pi_0(data, macro_pop)[macro_order]
exp_pi0  = c(farming = .659, manual = .233, nonmanual = .046, nonemp = .062)

if (!all(abs(pi0_nat - exp_pi0) < .002)) {
  stop(sprintf(
    "National pi_0 mismatch.\n  got:      %s\n  expected: %s",
    paste(round(pi0_nat, 3), collapse = " / "),
    paste(exp_pi0,           collapse = " / ")
  ))
}
cat("National pi_0 check passed:", paste(round(pi0_nat, 3), collapse = " / "), "\n")

################################################################################
# PER-REGION COMPUTATION
################################################################################

regions_list = c("sw", "south", "cali", "ok", "plains", "nw", "north")

region_display = c(
  sw     = "Southwest",
  south  = "South",
  cali   = "California",
  ok     = "Oklahoma",
  plains = "Plains",
  nw     = "Northwest",
  north  = "North"
)

compute_region = function(df_reg) {
  # p_matrix() uses w_atc_norm, which is already set to w_trim_norm globally.
  # Named indexing throughout — no positional access.
  P   = p_matrix(df_reg, macro_pop, macro_son, matrix = TRUE)
  pis = pi_star(P)

  eig_mods = sort(Mod(eigen(P, only.values = TRUE)$values), decreasing = TRUE)
  lambda2  = eig_mods[2]

  relief = weighted.mean(df_reg$empstatd_1940 == 11, df_reg$w_atc_norm,
                         na.rm = TRUE)
  list(
    P          = P,
    pis        = pis,
    farm_ret   = P["farming", "farming"],
    pi_farming = pis["farming"],
    pi_manual  = pis["manual"],
    pi_nonman  = pis["nonmanual"],
    lambda2    = lambda2,
    relief     = relief,
    n          = nrow(df_reg)
  )
}

results = setNames(
  lapply(regions_list, function(r) compute_region(filter(data, region == r))),
  regions_list
)

################################################################################
# ASSERTIONS — tolerance .002; fail loudly
################################################################################

expected = list(
  sw     = list(farm_ret=.623, pi_farming=.469, pi_manual=.386, pi_nonman=.032, lambda2=.313, relief=.110, n=2090L),
  south  = list(farm_ret=.739, pi_farming=.376, pi_manual=.482, pi_nonman=.072, lambda2=.682, relief=.062, n=1139L),
  cali   = list(farm_ret=.442, pi_farming=.350, pi_manual=.513, pi_nonman=.016, lambda2=.186, relief=.165, n= 645L),
  ok     = list(farm_ret=.459, pi_farming=.276, pi_manual=.405, pi_nonman=.113, lambda2=.377, relief=.136, n=2531L),
  plains = list(farm_ret=.370, pi_farming=.268, pi_manual=.538, pi_nonman=.067, lambda2=.170, relief=.255, n=2656L),
  nw     = list(farm_ret=.436, pi_farming=.243, pi_manual=.595, pi_nonman=.031, lambda2=.295, relief=.174, n=1155L),
  north  = list(farm_ret=.234, pi_farming=.106, pi_manual=.698, pi_nonman=.061, lambda2=.180, relief=.233, n=2030L)
)

tol = .002
num_fields = c("farm_ret", "pi_farming", "pi_manual", "pi_nonman", "lambda2", "relief")

for (r in regions_list) {
  exp = expected[[r]]
  got = results[[r]]
  for (fld in num_fields) {
    delta = abs(got[[fld]] - exp[[fld]])
    if (delta > tol)
      stop(sprintf("[%s] %s: got %.4f, expected %.4f (|delta| = %.4f > tol %.3f)",
                   r, fld, got[[fld]], exp[[fld]], delta, tol))
  }
  if (got$n != exp$n)
    stop(sprintf("[%s] n: got %d, expected %d", r, got$n, exp$n))
}
cat("All regional assertions passed.\n")

################################################################################
# BUILD TABLE DATA FRAME — sorted by pi_farming descending
################################################################################

tbl = do.call(rbind, lapply(regions_list, function(r) {
  g = results[[r]]
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

################################################################################
# SHADING — column range rescaled linearly to integer [2, 22]
# pi_nonmanual is intentionally left unshaded.
################################################################################

shade_int = function(x) {
  rng = range(x, na.rm = TRUE)
  if (diff(rng) < 1e-10) return(rep(2L, length(x)))
  as.integer(round(2 + (x - rng[1]) / diff(rng) * 20))
}

sh_farm   = shade_int(tbl$pi_farming)
sh_manual = shade_int(tbl$pi_manual)
sh_relief = shade_int(tbl$relief)

################################################################################
# LATEX CELL FORMATTERS
################################################################################

# Strip leading zero: "0.469" -> ".469"
fmt = function(x) sub("^0", "", sprintf("%.3f", x))

cell_sh = function(shade, val, bold = FALSE) {
  v = if (bold) sprintf("\\textbf{%s}", fmt(val)) else fmt(val)
  sprintf("\\sh{%d}{%s}", shade, v)
}

cell = function(val, bold = FALSE) {
  if (bold) sprintf("\\textbf{%s}", fmt(val)) else fmt(val)
}

################################################################################
# EMIT TABULAR
################################################################################

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

# N footnote line — sorted by n descending, separate from tabular
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
# ZERO-CELL CHECK AND BLUME-STYLE EPSILON PERTURBATION
################################################################################

eps_vals   = c(.001, .005, .01)
zero_thresh = 1e-8   # treat P[i,j] < thresh as structurally zero

perturb_P = function(P, eps) {
  Pp = P
  for (i in seq_len(nrow(P))) {
    row = P[i, ]
    low = row < eps
    if (!any(low)) next
    # Floor low cells at eps; remove the added mass proportionally from high cells
    added    = sum(eps - row[low])
    high     = !low
    high_sum = sum(row[high])
    row[low] = eps
    if (high_sum > added) {
      row[high] = row[high] * (high_sum - added) / high_sum
    }
    Pp[i, ] = row / sum(row)          # renormalise to machine precision
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
  P       = results[[r]]$P
  n_zero  = sum(P < zero_thresh)
  pis0    = pi_star(P)
  lam0    = sort(Mod(eigen(P, only.values = TRUE)$values), decreasing = TRUE)[2]

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
        eps, d_nm, d_l2
      ))
    }
  } else {
    cl = c(cl, sprintf("%% %s: no zero cells", region_display[r]))
  }
}

cl = c(cl, "")
cl = c(cl, sprintf(
  "%% Summary: max Delta pi_nonmanual = %.4f; max Delta lambda2 = %.4f",
  max_d_nonman, max_d_lambda2
))
cl = c(cl, "")

# Human-readable sentence for slide notes
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
