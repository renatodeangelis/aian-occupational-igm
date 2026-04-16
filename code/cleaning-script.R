library(dplyr)
library(readr)
library(tidyr)

source("code/utils.R")

aian_raw = read_csv(
  file = "https://www.dropbox.com/scl/fi/3gwg0wvb0sj0njjdjxg1p/clp_mlp1850_1940_linked_subsample_300raced_2022-11-5.csv?rlkey=qvub7x5cy6zsbgp74a90l381q&st=kdbmvn65&dl=1",
  col_types = cols(.default = col_integer(),
                   histid_1850 = col_character(),
                   histid_1860 = col_character(),
                   histid_1870 = col_character(),
                   histid_1880 = col_character(),
                   histid_1900 = col_character(),
                   histid_1910 = col_character(),
                   histid_1920 = col_character(),
                   histid_1930 = col_character(),
                   histid_1940 = col_character(),
                   histid_pop_1850 = col_character(),
                   histid_pop_1860 = col_character(),
                   histid_pop_1870 = col_character(),
                   histid_pop_1880 = col_character(),
                     histid_pop_1900 = col_character(),
                   histid_pop_1910 = col_character(),
                   histid_pop_1920 = col_character(),
                   histid_pop_1930 = col_character(),
                   histid_pop_1940 = col_character())) |>
  janitor::clean_names()

cat("Raw linked records:", nrow(aian_raw), "\n")

aian_clean = aian_raw |>
  select(where(~ !all(is.na(.))),
         -starts_with(c("bpld", "birthyr", "gqtyped", "raced", "school_pop",
                        "sex_pop", "sizepl")),
         -ends_with(c("1850", "1860", "1870", "1880")),
         -age_1900, -age_1910, -age_1920, -age_1930, -countyicp_1900,
         -countyicp_1910, -countyicp_1920, -countyicp_1930, -empstatd_1910,
         -empstatd_1930, -histid_1900, -histid_1910, -histid_1920, -histid_1930,
         -labforce_1910, -labforce_1920, -labforce_1930, -occ1950_1900,
         -occ1950_1910, -occ1950_1920, -occ1950_1930, -school_1900, -school_1910,
         -school_1920, -school_1930, -sex_1900, -sex_1910, -sex_1920, -sex_1930,
         -statefip_1900, -statefip_1910, -statefip_1920, -statefip_1930) |>
  rename(occ_son = occ1950_1940) |>
  filter(between(age_1940, 20, 44),
         sex_1940 == 1) |>
  (\(x) { cat("After age/sex filter:", nrow(x), "son-records\n"); x })() |>
  group_by(histid_1940) |>
  summarise(
    # Son's variables (age, occ, educ, etc.) and weights: identical across
    # duplicate linked records, so max is safe
    across(where(is.integer) & !contains("_pop_"), max),
    # Father's histid columns: paste to detect multiple fathers per son
    across(where(is.character),
           ~ paste(unique(na.omit(.x)), collapse = "; ")),
    # Father's integer characteristics: paste to detect conflicts, recover below
    across(where(is.integer) & contains("_pop_"),
           ~ paste(unique(na.omit(.x)), collapse = "; ")),
    .groups = "drop") |>
  mutate(
    # Empty strings from all-NA groups → NA
    across(where(is.character), ~ na_if(.x, "")),
    # Former integer columns were pasted; take first value and restore type
    across(contains("_pop_") & !starts_with("histid_pop"),
           ~ as.integer(sub(";.*", "", .x)))) |>
  filter(! if_any(starts_with("histid_pop_"),
                  ~ grepl(";", .x, fixed = TRUE))) |>
  (\(x) { cat("After multiple-father drop:", nrow(x), "unique sons\n"); x })() |>
  mutate(pid = coalesce(histid_pop_1940, histid_pop_1930, histid_pop_1920,
                        histid_pop_1910, histid_pop_1900))

cat("Unique fathers:", n_distinct(aian_clean$pid), "\n")

aian_age = aian_clean |>
  select(pid, starts_with("age_pop")) |>
  mutate(birthyr_1900 = 1900 - age_pop_1900,
         birthyr_1910 = 1910 - age_pop_1910,
         birthyr_1920 = 1920 - age_pop_1920,
         birthyr_1930 = 1930 - age_pop_1930,
         birthyr_1940 = 1940 - age_pop_1940) |>
  group_by(pid) |>
  summarise(across(starts_with("birthyr"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  rowwise() |>
  mutate(
    birth_median = round(median(c_across(starts_with("birthyr")), na.rm = TRUE)),
    spread = diff(range(c_across(starts_with("birthyr")), na.rm = TRUE)),
    spread_mad = median(abs(c_across(starts_with("birthyr")) - birth_median), na.rm = TRUE)) |>
  ungroup()

modal_occ_pick = pick_modal_occ(aian_clean, aian_age) |>
  rename(occ_pop = occ) |>
  select(-year)

cat("Fathers with recovered occupation:", nrow(modal_occ_pick),
    "of", n_distinct(aian_clean$pid), "unique fathers\n")

modal_occ_pick_attach = pick_modal_occ(aian_clean, aian_age,
                                        prefer_employed   = FALSE,
                                        empstatd_tiebreak = TRUE) |>
  rename(occ_pop_attach = occ, picked_year_attach = year) |>
  select(-birthyr_pop)

aian_merged = aian_clean |>
  left_join(modal_occ_pick, by = "pid") |>
  left_join(modal_occ_pick_attach, by = "pid") |>
  left_join(aian_age |> select(pid, birthyr_spread = spread, spread_mad), by = "pid") |>
  mutate(birthyr_son = 1940 - age_1940) |>
  select(-pid, -starts_with("age"), -starts_with("occ1950_pop")) |>
  filter(!is.na(occ_pop)) |>
  (\(x) { cat("After missing occ_pop drop:", nrow(x), "father-son pairs\n"); x })() |>
  filter(is.na(spread_mad) | spread_mad <= 4) |>
  (\(x) { cat("After spread filter (MAD <= 4):", nrow(x), "father-son pairs\n"); x })() |>
  mutate(spread_flag = !is.na(spread_mad) & spread_mad > 2,
         attachment_level_son = case_when(
           occ_son <= 970 & empstatd_1940 %in% c(21, 22, 31, 32, 33, 34) & wkswork1_1940 == 0  ~ 1L,
           occ_son <= 970 & empstatd_1940 %in% c(21, 22, 31, 32, 33, 34) & wkswork1_1940 <= 13 ~ 2L,
           occ_son <= 970 & empstatd_1940 %in% c(21, 22, 31, 32, 33, 34) & wkswork1_1940 <= 26 ~ 3L,
           occ_son <= 970 & empstatd_1940 %in% c(21, 22, 31, 32, 33, 34)                       ~ 4L,
           .default = NA_integer_),
         attachment_level_pop = {
           ep <- case_when(
             picked_year_attach == 1910 ~ empstatd_pop_1910,
             picked_year_attach == 1930 ~ empstatd_pop_1930,
             picked_year_attach == 1940 ~ empstatd_pop_1940,
             .default = NA_integer_
           )
           lf <- if_else(picked_year_attach == 1920, labforce_pop_1920, NA_integer_)
           case_when(
             occ_pop_attach <= 970 & ep %in% c(21,22,31,32,33,34) ~ 1L,
             occ_pop_attach <= 970 & is.na(ep) & lf == 1          ~ 2L,
             .default = NA_integer_
           )
         },
         meso_pop = classify_meso(occ_pop),
         macro_pop = classify_macro(meso_pop),
         meso_son = classify_meso(occ_son),
         macro_son = classify_macro(meso_son),
         # Attachment-matrix father classification (no occ <= 970 preference)
         meso_pop_attach = classify_meso(occ_pop_attach),
         macro_pop_attach = classify_macro(meso_pop_attach)) |>
  mutate(across(starts_with("macro_"),
                ~ factor(.x, levels = macro_order, ordered = TRUE)),
         across(starts_with("meso_"),
                ~ factor(.x, levels = meso_order, ordered = TRUE))) |>
  mutate(
    lit_son = do.call(pmax, c(pick(starts_with("lit_19")), na.rm = TRUE)),
    lit_pop = do.call(pmax, c(pick(starts_with("lit_pop")), na.rm = TRUE))) |>
  filter(school_1940 == 1) |>
  (\(x) { cat("After school filter:", nrow(x), "father-son pairs\n"); x })() |>
  select(-starts_with("lit_19"), -starts_with("lit_pop_"), -sex_1940, -school_1940) |>
  relocate(statefip_1940, .after = histid_1940) |>
  relocate(starts_with("statefip_pop"), .after = histid_pop_1940) |>
  relocate(birthyr_son, .after = histid_1940) |>
  relocate(birthyr_pop, .after = histid_pop_1940) |>
  relocate(birthyr_spread, .after = birthyr_pop) |>
  relocate(spread_mad, .after = birthyr_spread) |>
  relocate(spread_flag, .after = spread_mad) |>
  relocate(attachment_level_son, .after = spread_flag) |>
  relocate(attachment_level_pop, .after = attachment_level_son) |>
  relocate(occ_pop_attach, picked_year_attach, .after = occ_pop) |>
  relocate(starts_with("macro_son"), .after = occ_son) |>
  relocate(starts_with("meso_son"), .after = macro_son_alt) |>
  relocate(lit_son, .after = educd_1940) |>
  relocate(lit_pop, .after = educd_pop_1940) |>
  relocate(starts_with("w_parent"), .after = last_col())

cat("\nFinal analysis sample:", nrow(aian_merged), "father-son pairs\n")

write_csv(aian_merged, "data/aian_merged.csv")
