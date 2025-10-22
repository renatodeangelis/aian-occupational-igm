library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(sf)

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
  select(-starts_with("gqtyped"),
         -starts_with("birthyr"),
         -starts_with("raced"),
         -ends_with("1850"),
         -ends_with("1860"),
         -ends_with("1870"),
         -ends_with("1880"),
         -occ1950_1900,
         -occ1950_1910,
         -occ1950_1920,
         -occ1950_1930,
         -histid_1900,
         -histid_1910,
         -histid_1920,
         -histid_1930,
         -age_1900,
         -age_1910,
         -age_1920,
         -age_1930,
         -bpld_1900,
         -bpld_1910,
         -bpld_1920,
         -bpld_1930,
         -countyicp_1900,
         -countyicp_1910,
         -countyicp_1920,
         -countyicp_1930,
         -empstatd_1910,
         -empstatd_1930,
         -statefip_1900,
         -statefip_1910,
         -statefip_1920,
         -statefip_1930,
         -urban_1900,
         -urban_1910,
         -urban_1920,
         -urban_1930) |>
  rename(son_occ = occ1950_1940) |>
  filter(age_1940 >= 20 & age_1940 < 45,
         sex_1940 == 1,
         (age_1940 < age_pop_1940 + 10) | is.na(age_pop_1940)) |>
  group_by(histid_1940) |>
  summarise(across(starts_with("w_parent"), ~ max(.x)),
            across(6:37, ~ max(.x)),
            across(38:42, ~ paste(unique(na.omit(.x)), collapse = "; "),
                   .names = "{.col}"),
            across(43:100, ~ paste(unique(na.omit(.x)), collapse = "; "),
                   .names = "{.col}"),
            .groups = "drop") |>
  mutate(across(39:101, ~ na_if(.x, "")),
         across(44:101, as.integer)) |>
  filter(! if_any(starts_with("histid_pop_"), 
                  ~ grepl(";", .x, fixed = TRUE)))

df = aian_clean |>
  select(starts_with("histid_pop_"),
         starts_with("occ1950_pop_"),
         starts_with("age_pop_")) |>
  # person id = first non-missing histid across years
  mutate(pid = coalesce(histid_pop_1940, histid_pop_1930, histid_pop_1920,
                        histid_pop_1910, histid_pop_1900))

occ_long = df |>
  pivot_longer(
    cols = starts_with("occ1950_pop_"),
    names_to = "year",
    names_pattern = "occ1950_pop_(\\d{4})",
    names_transform = list(year = as.integer),
    values_to = "occ",
    values_drop_na = TRUE
  ) |>
  select(pid, year, occ) |>
  group_by(pid, year) |>
  summarise(
    occ = {
      vals = unique(na.omit(occ))
      if (length(vals) == 0) NA_integer_
      else if (any(vals <= 970)) vals[vals <= 970][1]
      else vals[1]
    },
    .groups = "drop"
  )

age_long = df |>
  pivot_longer(
    cols = starts_with("age_pop_"),
    names_to = "year",
    names_pattern = "age_pop_(\\d{4})",
    names_transform = list(year = as.integer),
    values_to = "age",
    values_drop_na = TRUE
  ) |>
  select(pid, year, age) |>
  group_by(pid, year) |>
  summarise(age = dplyr::first(na.omit(age), default = NA_integer_), .groups = "drop")

modal_input = occ_long |>
  left_join(age_long, by = c("pid", "year")) |> 
  group_by(pid) |>
  mutate(
    # Is there at least one <= 970 for this person?
    has_pref = any(occ <= 970, na.rm = TRUE),
    occ_used = if_else(has_pref & occ <= 970, occ, if_else(has_pref, NA_integer_, occ))
  ) |>
  # Drop rows that aren't in the eligible pool
  filter(!is.na(occ_used)) |>
  # Mode: count frequency of occ_used within pid
  add_count(pid, occ_used, name = "freq") |>
  filter(freq == max(freq)) |>                 # keep modal ties only
  # Tie-breaker: closest age to 40 (NAs last)
  mutate(age_dist = abs(age - 40)) |>
  arrange(pid, age_dist) |>
  slice_head(n = 1) |>                         # one row per pid
  ungroup() |>
  transmute(pid, modal_occ = occ_used,
            picked_year = year, age_at_pick = age) |>
  filter(age_at_pick <= 65)

aian_merged = aian_clean |>
  mutate(across(starts_with("histid_pop"), as.character),
         pid = coalesce(histid_pop_1940, histid_pop_1930, histid_pop_1920,
                        histid_pop_1910, histid_pop_1900)) |>
  left_join(modal_input, by = "pid") |>
  select(-pid) |>
  filter(!is.na(modal_occ)) |>
  mutate(dad_meso = case_when(
    modal_occ == 100 | modal_occ == 123 ~ "farming",
    modal_occ <= 99 | (modal_occ >= 200 & modal_occ <= 290) ~ "prof",
    modal_occ >= 300 & modal_occ <= 490 ~ "clerical",
    (modal_occ >= 500 & modal_occ <= 594) | modal_occ == 762 |
      modal_occ == 773 | modal_occ == 781 | modal_occ == 782 ~ "crafts",
    modal_occ >= 595 & modal_occ <= 970 & modal_occ != 762 &
      modal_occ != 773 & modal_occ != 781 & modal_occ != 782 ~ "unskilled",
    modal_occ > 970 ~ "unemp")) |>
  mutate(dad_macro = case_when(
    dad_meso == "farming" ~ "farming",
    dad_meso == "prof" | dad_meso == "clerical" ~ "white_col",
    dad_meso == "crafts" | dad_meso == "unskilled" ~ "blue_col",
    dad_meso == "unemp" ~ "unemp")) |>
  mutate(son_meso = case_when(
    son_occ == 100 | son_occ == 123 ~ "farming",
    son_occ <= 99 | (son_occ >= 200 & son_occ <= 290) ~ "prof",
    son_occ >= 300 & son_occ <= 490 ~ "clerical",
    (son_occ >= 500 & son_occ <= 594) | son_occ == 762 |
      son_occ == 773 | son_occ == 781 | son_occ == 782 ~ "crafts",
    son_occ >= 595 & son_occ <= 970 & son_occ != 762 &
      son_occ != 773 & son_occ != 781 & son_occ != 782 ~ "unskilled",
    son_occ > 970 ~ "unemp")) |>
  mutate(son_macro = case_when(
    son_meso == "farming" ~ "farming",
    son_meso == "prof" | son_meso == "clerical" ~ "white_col",
    son_meso == "crafts" | son_meso == "unskilled" ~ "blue_col",
    son_meso == "unemp" ~ "unemp"))

write_csv(aian_merged, "data/aian_merged.csv")
