library(dplyr)
library(tidyr)
library(purrr)

source("code/utils.R")

set.seed(20260803)
B            = 2000
region_order = c("north", "plains", "nw", "cali", "ok", "sw", "south")

dat = readRDS("data/aian_weighted.rds") |>
  mutate(
    w_atc_norm = w_trim_norm,
    relief     = !is.na(empstatd_1940) & empstatd_1940 == 11,
    # farming wins over relief: a relief worker coded to a farm occupation
    # has NOT left agriculture. Flip this if you prefer relief to dominate.
    dest5 = case_when(
      macro_son == "farming"   ~ "farming",
      relief                   ~ "relief",
      macro_son == "manual"    ~ "manual",
      macro_son == "nonmanual" ~ "nonmanual",
      macro_son == "nonemp"    ~ "nonemp"
    )
  )

# how many cases does that precedence choice affect?
n_relief_farm = dat |> filter(relief, macro_son == "farming") |> nrow()
cat("relief workers coded to farming:", n_relief_farm, "\n")

exit_comp = function(df, dest) {
  fo = filter(df, macro_pop == "farming")
  if (nrow(fo) < 30) return(NULL)
  lv = filter(fo, .data[[dest]] != "farming")
  wl = sum(lv$w_atc_norm)
  tibble(
    n_orig    = nrow(fo),
    exit_rate = wl / sum(fo$w_atc_norm),
    to_manual = sum(lv$w_atc_norm[lv[[dest]] == "manual"])    / wl,
    to_nonman = sum(lv$w_atc_norm[lv[[dest]] == "nonmanual"]) / wl,
    to_nonemp = sum(lv$w_atc_norm[lv[[dest]] == "nonemp"])    / wl,
    to_relief = if ("relief" %in% lv[[dest]])
      sum(lv$w_atc_norm[lv[[dest]] == "relief"]) / wl
    else
      NA_real_
  )
}

specs = list(
  baseline = function(df) exit_comp(df, "macro_son"),
  holdout  = function(df) exit_comp(filter(df, !relief), "macro_son"),
  fivecat  = function(df) exit_comp(df, "dest5")
)

point = imap_dfr(specs, \(f, nm)
  dat |>
    group_by(region) |>
    group_modify(~ f(.x) %||% tibble()) |>
    ungroup() |>
    mutate(spec = nm))

boot_one = function(df) {
  s = df[sample(nrow(df), nrow(df), replace = TRUE), ]
  imap_dfr(specs, \(f, nm) {
    r = f(s)
    if (is.null(r)) return(tibble())
    mutate(r, spec = nm)
  })
}

boot = dat |>
  group_by(region) |>
  group_modify(\(.x, .y) map_dfr(seq_len(B), \(b) mutate(boot_one(.x), rep = b))) |>
  ungroup()

ci = boot |>
  group_by(region, spec) |>
  summarise(across(c(exit_rate, to_manual, to_nonemp),
                   list(lo = \(x) quantile(x, .025, na.rm = TRUE),
                        hi = \(x) quantile(x, .975, na.rm = TRUE))),
            .groups = "drop")

# does the ordering survive? rank correlation and its bootstrap distribution
rank_cor = boot |>
  select(region, spec, rep, exit_rate) |>
  pivot_wider(names_from = spec, values_from = exit_rate) |>
  group_by(rep) |>
  summarise(rho = cor(baseline, holdout, method = "spearman"), .groups = "drop")

print(quantile(rank_cor$rho, c(.025, .5, .975)))

rank_cor_ne <- boot %>%
  select(region, spec, rep, to_nonemp) %>%
  pivot_wider(names_from = spec, values_from = to_nonemp) %>%
  group_by(rep) %>%
  summarise(rho = cor(baseline, holdout, method = "spearman"), .groups = "drop")
quantile(rank_cor_ne$rho, c(.025, .5, .975))