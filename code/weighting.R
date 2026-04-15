library(dplyr)
library(readr)
library(tidyr)
library(cobalt)
library(ggplot2)

source("code/utils.R")

aian_merged = read_csv("data/aian_merged.csv") |>
  mutate(region = assign_region(statefip_1940),
         education = classify_education(educd_1940))

aian_full = read_csv(
  file = "https://www.dropbox.com/scl/fi/rdgea171s0uztwasmojkx/usa_00022.csv?rlkey=msbupznr1hq12ejf305vdk2kx&st=s1qmh2in&dl=1") |>
  janitor::clean_names() |>
  filter(age >= 20 & age < 45,
         school == 1) |>
  mutate(birthyr_son = 1940 - age,
         meso_son = classify_meso(occ1950),
         macro_son = classify_macro(meso_son),
         meso_son_alt = classify_meso(occ1950, split_farmer = FALSE),
         macro_son_alt = classify_macro(meso_son_alt),
         region = assign_region(statefip),
         education = classify_education(educd)) |>
  rename(statefip_1940 = statefip, urban_1940 = urban)

# --- Point estimate weights via compute_weights() ---
weights_out = compute_weights(aian_merged, aian_full)
aian_ps     = weights_out$data

# --- Common support diagnostics (reuses p_hat from compute_weights) ---
aian_comb = bind_rows(
  aian_ps   |> mutate(linked = 1),
  aian_full |> mutate(linked = 0, p_hat = weights_out$p_hat_full)
) |>
  mutate(
    cohort    = cut(birthyr_son,
                    breaks = c(1895, 1900, 1905, 1910, 1915, 1921),
                    labels = c("1896-1900", "1901-1905", "1906-1910",
                               "1911-1915", "1916-1920")),
    region    = as.factor(region),
    education = as.factor(education))

p_support = ggplot(aian_comb,
    aes(x = p_hat, fill = factor(linked, labels = c("Unlinked", "Linked")))) +
  geom_density(alpha = 0.5) +
  labs(x = "Propensity score", y = "Density", fill = "Sample") +
  theme_minimal()
ggsave("figures/ps_common_support.png", p_support, width = 7, height = 4)

ps_range_linked   = range(aian_comb$p_hat[aian_comb$linked == 1])
ps_range_unlinked = range(aian_comb$p_hat[aian_comb$linked == 0])
cat("PS range — linked:", round(ps_range_linked, 4),
    " unlinked:", round(ps_range_unlinked, 4), "\n")

cs_lower = max(ps_range_linked[1], ps_range_unlinked[1])
cs_upper = min(ps_range_linked[2], ps_range_unlinked[2])
cat("Common support region:", round(cs_lower, 4), "to",
    round(cs_upper, 4), "\n")

n_outside_cs = sum(aian_ps$p_hat < cs_lower | aian_ps$p_hat > cs_upper)
cat("Linked obs outside common support:", n_outside_cs, "of", nrow(aian_ps),
    sprintf("(%.1f%%)\n", 100 * n_outside_cs / nrow(aian_ps)))

# --- 2.5 fix: Weight trimming and diagnostics ---
cat("\n--- Weight diagnostics (untrimmed) ---\n")
cat("Summary of w_atc_norm:\n")
print(summary(aian_ps$w_atc_norm))

ess_global = sum(aian_ps$w_atc_norm)^2 / sum(aian_ps$w_atc_norm^2)
cat("Effective sample size:", round(ess_global, 1),
    "of", nrow(aian_ps), "observations\n")

trim_threshold = quantile(aian_ps$w_atc, 0.99)
cat("99th percentile trim threshold:", round(trim_threshold, 3), "\n")
cat("Observations trimmed:", sum(aian_ps$w_atc > trim_threshold), "\n")

aian_ps = aian_ps |>
  mutate(w_trim = pmin(w_atc, trim_threshold),
         w_trim_norm = w_trim * n() / sum(w_trim))

ess_trimmed = sum(aian_ps$w_trim_norm)^2 / sum(aian_ps$w_trim_norm^2)
cat("ESS after trimming:", round(ess_trimmed, 1), "\n")

# --- 2.6 fix: Covariate balance diagnostics ---
aian_comb_bal = aian_comb |>
  left_join(aian_ps |> select(histid_1940, w_atc_norm), by = "histid_1940") |>
  mutate(w_atc_norm = if_else(linked == 0, 1, w_atc_norm))

bt = bal.tab(linked ~ cohort + region + education + urban_1940,
             data = aian_comb_bal, weights = "w_atc_norm",
             method = "weighting", estimand = "ATC",
             un = TRUE)
cat("\n--- Covariate balance (substantive variables) ---\n")
print(bt)
dir.create("output", showWarnings = FALSE)
write_csv(as.data.frame(bt$Balance), "output/balance_table.csv")

bt_state = bal.tab(linked ~ statefip_1940,
                   data = aian_comb_bal, weights = "w_atc_norm",
                   method = "weighting", estimand = "ATC",
                   un = TRUE)
state_smds = bt_state$Balance$Diff.Adj
cat("\n--- State balance summary ---\n")
cat("State balance — max |SMD|:", round(max(abs(state_smds)), 3),
    " mean |SMD|:", round(mean(abs(state_smds)), 3), "\n")
if (any(abs(state_smds) > 0.2)) {
  bad_states = rownames(bt_state$Balance)[abs(state_smds) > 0.2]
  cat("WARNING: States with |SMD| > 0.2:",
      paste(bad_states, collapse = ", "), "\n")
}
write_csv(as.data.frame(bt_state$Balance), "output/balance_table_state.csv")

# --- Finalize and write ---
aian_ps = aian_ps |>
  select(-p_hat, -w_atc, -w_trim) |>
  select(where(~ !all(is.na(.)))) |>
  relocate(starts_with("w_parent"), .after = last_col()) |>
  relocate(w_atc_norm, w_trim_norm, .before = starts_with("w_parent"))

write_csv(aian_ps, "data/aian_weighted.csv")
