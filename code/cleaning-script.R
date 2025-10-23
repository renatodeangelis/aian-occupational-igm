library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)

aian_raw = readr::read_csv(
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
  select(where(~ !all(is.na(.)))) |>
  select(-starts_with("bpld"),
         -starts_with("birthyr"),
         -starts_with("gqtyped"),
         -starts_with("raced"),
         -starts_with("school_pop"),
         -starts_with("sex_pop"),
         -starts_with("sizepl"),
         -starts_with("urban"),
         -ends_with("1850"),
         -ends_with("1860"),
         -ends_with("1870"),
         -ends_with("1880"),
         -age_1900,
         -age_1910,
         -age_1920,
         -age_1930,
         -countyicp_1900,
         -countyicp_1910,
         -countyicp_1920,
         -countyicp_1930,
         -empstatd_1910,
         -empstatd_1930,
         -histid_1900,
         -histid_1910,
         -histid_1920,
         -histid_1930,
         -labforce_1910,
         -labforce_1920,
         -labforce_1930,
         -occ1950_1900,
         -occ1950_1910,
         -occ1950_1920,
         -occ1950_1930,
         -school_1900,
         -school_1910,
         -school_1920,
         -school_1930,
         -sex_1900,
         -sex_1910,
         -sex_1920,
         -sex_1930,
         -statefip_1900,
         -statefip_1910,
         -statefip_1920,
         -statefip_1930) |>
  rename(son_occ = occ1950_1940) |>
  filter(age_1940 >= 20 & age_1940 < 45,
         sex_1940 == 1) |>
  group_by(histid_1940) |>
  summarise(across(starts_with("w_parent"), ~ max(.x)),
            across(6:19, ~ max(.x)),
            across(20:24, ~ paste(unique(na.omit(.x)), collapse = "; "),
                   .names = "{.col}"),
            across(25:57, ~ paste(unique(na.omit(.x)), collapse = "; "),
                   .names = "{.col}"),
            .groups = "drop") |>
  mutate(across(21:58, ~ na_if(.x, "")),
         across(26:58, as.integer)) |>
  filter(! if_any(starts_with("histid_pop_"), 
                  ~ grepl(";", .x, fixed = TRUE)))

aian_age = aian_clean |>
  mutate(pid = coalesce(histid_pop_1940, histid_pop_1930, histid_pop_1920,
                        histid_pop_1910, histid_pop_1900)) |>
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
  ungroup() |>
  mutate(quality = case_when(
    is.na(spread) ~ "no data",
    spread <= 2 ~ "good",
    spread <= 5 ~ "ok",
    TRUE ~ "poor"))

aian_ages = aian_clean |>
  mutate(pid = coalesce(histid_pop_1940, histid_pop_1930, histid_pop_1920,
                        histid_pop_1910, histid_pop_1900),
         birthyr_son = 1940 - age_1940) |>
  left_join(aian_age |> select(pid, birth_median), by = "pid") |>
  rename(birthyr_father = birth_median) |>
  select(-starts_with("age")) |>
  filter(birthyr_son > birthyr_father + 15)

df = aian_ages |>
  select(starts_with("histid_pop_"),
         starts_with("occ1950_pop_"),
         starts_with("age_pop_")) |>
  mutate(pid = coalesce(histid_pop_1940, histid_pop_1930, histid_pop_1920,
                        histid_pop_1910, histid_pop_1900))

modal_input = aian_ages |>
  select(pid, starts_with("occ1950_pop_"), birthyr_father) |>
  pivot_longer(
    cols = starts_with("occ1950_pop_"),
    names_to = "year",
    names_pattern = "occ1950_pop_(\\d{4})",
    names_transform = list(year = as.integer),
    values_to = "occ",
    values_drop_na = TRUE) |>
  group_by(pid, year) |>
  summarise(
    birthyr_father = dplyr::first(birthyr_father),
    occ = {
      vals <- unique(na.omit(occ))
      if (length(vals) == 0) NA_integer_
      else if (any(vals <= 970)) vals[vals <= 970][1] else vals[1]
    },
    .groups = "drop") |>
  mutate(implied_age = ifelse(!is.na(birthyr_father), year - birthyr_father, NA_real_)) |>
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
            modal_occ = occ_used,
            picked_year = year,
            age_at_pick = implied_age,
            birthyr_father)

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

aian_merged = aian_ages |>
  left_join(modal_input, by = c("pid", "birthyr_father")) |>
  select(-pid,
         -starts_with("occ1950_pop")) |>
  filter(!is.na(modal_occ)) |>
  mutate(meso_pop = classify_meso(modal_occ),
         macro_pop = classify_macro(meso_pop),
         meso_son = classify_meso(son_occ),
         macro_son = classify_macro(meso_son),
         # Alternative classification of farmworkers
         meso_pop_alt = classify_meso(modal_occ, split_farmer = FALSE),
         macro_pop_alt = classify_macro(meso_pop_alt),
         meso_son_alt = classify_meso(son_occ, split_farmer = FALSE),
         macro_son_alt = classify_macro(meso_son_alt)) |>
  mutate(across(starts_with("macro_"),
                ~ factor(.x, levels = macro_order, ordered = TRUE)),
         across(starts_with("meso_") & !ends_with("_alt"),
                ~ factor(.x, levels = meso_order, ordered = TRUE)),
         across(starts_with("meso_") & ends_with("_alt"),
                ~ factor(.x, levels = meso_order_alt, ordered = TRUE)))

write_csv(aian_merged, "data/aian_merged.csv")
