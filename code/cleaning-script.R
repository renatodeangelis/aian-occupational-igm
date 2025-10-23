library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)

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

aian_clean = aian_raw |>
  select(where(~ !all(is.na(.))),
         -starts_with(c("bpld", "birthyr", "gqtyped", "raced", "school_pop",
                        "sex_pop", "sizepl", "urban", "wkswork1")),
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
  group_by(histid_1940) |>
  summarise(across(starts_with("w_parent"), ~ max(.x)),
            across(6:18, ~ max(.x)),
            across(19:23, ~ paste(unique(na.omit(.x)), collapse = "; "),
                   .names = "{.col}"),
            across(24:55, ~ paste(unique(na.omit(.x)), collapse = "; "),
                   .names = "{.col}"),
            .groups = "drop") |>
  mutate(across(20:56, ~ na_if(.x, "")),
         across(25:56, as.integer)) |>
  filter(! if_any(starts_with("histid_pop_"), 
                  ~ grepl(";", .x, fixed = TRUE))) |>
  mutate(pid = coalesce(histid_pop_1940, histid_pop_1930, histid_pop_1920,
                        histid_pop_1910, histid_pop_1900))

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
    n_ages = sum(!is.na(c_across(starts_with("birthyr")))),
    birth_median = round(median(c_across(starts_with("birthyr")), na.rm = TRUE)),
    spread = diff(range(c_across(starts_with("birthyr")), na.rm = TRUE))) |>
  ungroup()

modal_occ_pick = aian_clean |>
  mutate(birthyr_son = 1940 - age_1940) |>
  left_join(aian_age |> select(pid, birth_median), by = "pid") |>
  rename(birthyr_pop = birth_median) |>
  select(-starts_with("age")) |>
  filter(birthyr_son > birthyr_pop + 15) |>
  select(pid, starts_with("occ1950_pop_"), birthyr_pop) |>
  pivot_longer(
    cols = starts_with("occ1950_pop_"),
    names_to = "year",
    names_pattern = "occ1950_pop_(\\d{4})",
    names_transform = list(year = as.integer),
    values_to = "occ",
    values_drop_na = TRUE) |>
  group_by(pid, year) |>
  summarise(
    birthyr_pop = dplyr::first(birthyr_pop),
    occ = {
      vals <- unique(na.omit(occ))
      if (length(vals) == 0) NA_integer_
      else if (any(vals <= 970)) vals[vals <= 970][1] else vals[1]
    },
    .groups = "drop") |>
  mutate(implied_age = ifelse(!is.na(birthyr_pop), year - birthyr_pop, NA_real_)) |>
  group_by(pid) |>
  mutate(
    has_pref = any(occ <= 970, na.rm = TRUE),
    occ_used = if_else(has_pref & occ <= 970, occ,
                       if_else(has_pref, NA_integer_, occ))) |>
  filter(!is.na(occ_used)) |>
  add_count(pid, occ_used, name = "freq") |>
  filter(freq == max(freq)) |>
  mutate(age_dist = dplyr::coalesce(abs(implied_age - 40), Inf)) |>
  arrange(pid, age_dist, year) |>
  filter(dplyr::coalesce(implied_age, Inf) <= 65) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(pid,
            occ_pop = occ_used,
            picked_year = year,
            age_at_pick = implied_age,
            birthyr_pop)

classify_meso = function(occ, split_farmer = TRUE) {
  farmer_codes = c(100, 123)
  prof_codes = c(0:99, 200:290)
  crafts_codes = c(762, 773, 781, 782)
  
  case_when(
    occ %in% farmer_codes ~ "farmer",
    split_farmer & occ %in% 810:840 ~ "farmworker",
    occ %in% prof_codes ~ "prof",
    occ %in% 300:490 ~ "clerical",
    occ %in% 500:594 | occ %in% crafts_codes ~ "crafts",
    occ %in% 595:970 & !(occ %in% crafts_codes) & !(split_farmer & occ %in% 810:840) ~ "unskilled",
    occ > 970 ~ "unemp"
  )
}

classify_macro = function(meso) {
  case_when(
    meso %in% c("farmer", "farmworker") ~ "farming",
    meso %in% c("prof", "clerical") ~ "nonmanual",
    meso %in% c("crafts", "unskilled") ~ "manual",
    meso == "unemp" ~ "unemp")
}

macro_order = c("farming", "manual", "nonmanual", "unemp")
meso_order = c("farmworker", "farmer", "unskilled", "crafts", "clerical", "prof", "unemp")
meso_order_alt = c("farmer", "unskilled", "crafts", "clerical", "prof", "unemp")

aian_merged = aian_clean |>
  left_join(modal_occ_pick, by = c("pid")) |>
  mutate(birthyr_son = 1940 - age_1940) |>
  select(-pid, -picked_year, -starts_with("age"), -starts_with("occ1950_pop")) |>
  filter(!is.na(occ_pop)) |>
  mutate(meso_pop = classify_meso(occ_pop),
         macro_pop = classify_macro(meso_pop),
         meso_son = classify_meso(occ_son),
         macro_son = classify_macro(meso_son),
         # Alternative classification of farmworkers
         meso_pop_alt = classify_meso(occ_pop, split_farmer = FALSE),
         macro_pop_alt = classify_macro(meso_pop_alt),
         meso_son_alt = classify_meso(occ_son, split_farmer = FALSE),
         macro_son_alt = classify_macro(meso_son_alt)) |>
  mutate(across(starts_with("macro_"),
                ~ factor(.x, levels = macro_order, ordered = TRUE)),
         across(starts_with("meso_") & !ends_with("_alt"),
                ~ factor(.x, levels = meso_order, ordered = TRUE)),
         across(starts_with("meso_") & ends_with("_alt"),
                ~ factor(.x, levels = meso_order_alt, ordered = TRUE))) |>
  rowwise() |>
  mutate(
    lit_son = {
      v = c_across(starts_with("lit_19"))
      if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)},
    lit_pop = {
      v = c_across(starts_with("lit_pop"))
      if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)}) |>
  ungroup() |>
  select(-starts_with("lit_19"), -starts_with("lit_pop_"), -sex_1940) |>
  relocate(statefip_1940, .after = histid_1940) |>
  relocate(starts_with("statefip_pop"), .after = histid_pop_1940) |>
  relocate(birthyr_son, .after = histid_1940) |>
  relocate(birthyr_pop, .after = histid_pop_1940) |>
  relocate(starts_with("macro_son"), .after = occ_son) |>
  relocate(starts_with("meso_son"), .after = macro_son_alt) |>
  relocate(lit_son, .after = educd_1940) |>
  relocate(school_1940, .after = lit_son) |>
  relocate(lit_pop, .after = educd_pop_1940) |>
  relocate(starts_with("w_parent"), .after = last_col())

#write_csv(aian_merged, "data/aian_merged.csv")
