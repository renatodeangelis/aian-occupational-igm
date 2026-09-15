library(dplyr)
library(readr)
library(tidyr)
library(janitor)

source("code/00_utils.R")

aian = 3
pop_years = c(1900, 1910, 1920, 1930, 1940)
son_years = c(1900, 1910, 1920, 1930, 1940, 1950)

son_multiyear = c("lit", "gq", "gqtype", "school", "relate", "age",
                  "statefip", "speakeng", "occ1950")
son_1940only  = c("educ", "educd", "sex", "countyicp", "urban", "metro",
                  "empstat", "empstatd", "labforce", "classwkr", "ind1950",
                  "farm", "ownershp", "marst", "bpl", "birthyr", "classwkrd",
                  "wkswork1", "hrswork1", "hrswork2", "durunemp",
                  "incwage", "incnonwg", "migrate5", "migrate5d", "migplac5")
pop_multiyear = c("histid", "hik", "age", "birthyr", "occ1950", "ind1950",
                  "classwkr", "labforce", "empstat", "empstatd", "lit",
                  "speakeng", "farm", "ownershp", "gq", "gqtype", "relate",
                  "marst", "bpl", "educd", "statefip", "countyicp", "urban",
                  "wkswork1", "durunemp", "incwage", "classwkrd")

path = "https://www.dropbox.com/scl/fi/3x5hlb10sza5gyoqk1l0k/usa_00028.csv?rlkey=id665yubz25czcttpzrp7f5qv&st=p477m379&dl=1"

raw = read_csv(path, col_types = cols(.default = col_character())) |>
  clean_names() |>
  type_convert(col_types = cols(histid = col_character(),
                                hik = col_character(),
                              .default = col_guess()))

sons = raw |>
  filter(year == 1940, sex == 1, between(age, 20, 49), race == aian,
         !is.na(hik), hik != "") |>
  transmute(son_hik = hik, histid_1940 = histid)

cat("Sons (AIAN m 20-44, 1940, linked):", nrow(sons), "\n")

son_records = raw |> semi_join(sons, by = c("hik" = "son_hik")) |>
  rename(son_hik = hik)

father_records  = son_records |>
  filter(poploc > 0, year %in% pop_years) |>
  select(son_hik, year, serial, son_poploc = poploc, son_age = age) |>
  left_join(raw |> select(year, serial, pernum, everything()) |>
              rename_with(~ paste0("dad_", .x),
                          -c(year, serial, pernum)),
            by = c("year", "serial", "son_poploc" = "pernum"))

unresolved = sum(is.na(father_records$dad_histid))
cat("Unresolved POPLOC pointers:", unresolved, " <- must be 0; nonzero means household members are incomplete\n")

father_records = father_records |>
  filter(!is.na(dad_histid), dad_sex == 1)

multi = father_records |>
  group_by(son_hik) |>
  summarise(n_dad = n_distinct(dad_hik[!is.na(dad_hik) & dad_hik != ""]), .groups = "drop")

father_records = father_records |> semi_join(filter(multi, n_dad <= 1), by = "son_hik")
cat("After multiple-father drop:", n_distinct(father_records$son_hik), "sons\n")

father_records = father_records |>
  group_by(son_hik) |>
  mutate(pid = coalesce(first(na_if(dad_hik, "")), paste0("s_", son_hik))) |>
  ungroup()

# HIK-based lookup: recover father records for years when son was not co-residing.
# dad_hik_map maps each known father hik to the son(s) it belongs to and the pid.
dad_hik_map = father_records |>
  filter(!is.na(dad_hik), dad_hik != "") |>
  distinct(son_hik, dad_hik, pid)

hik_extra = raw |>
  filter(year %in% pop_years,
         hik %in% unique(dad_hik_map$dad_hik),
         sex == 1) |>
  inner_join(dad_hik_map, by = c("hik" = "dad_hik"),
             relationship = "many-to-many") |>
  mutate(.source = "hik") |>
  rename_with(~ paste0("dad_", .x), -c(son_hik, pid, year, .source))

cat(sprintf(
  "HIK lookup: %d father-year records found, %d unique fathers, %d unique sons\n",
  nrow(hik_extra),
  n_distinct(hik_extra$dad_hik),
  n_distinct(hik_extra$son_hik)))

# Bind poploc-derived and hik-derived records.
# arrange(desc(.source)) puts "poploc" before "hik" so widen()'s slice_head()
# keeps the poploc-validated row when both sources find the same father-year.
father_records = father_records |>
  mutate(.source = "poploc") |>
  bind_rows(hik_extra) |>
  arrange(desc(.source)) |>
  select(-.source)

widen = function(df, vars, years, suffix, id) {
  vars = intersect(vars, names(df))
  df |>
    filter(year %in% years) |>
    select(all_of(c(id, "year", vars))) |>
    group_by(across(all_of(c(id, "year")))) |>
    slice_head(n = 1) |>
    ungroup() |>
    pivot_wider(id_cols = all_of(id), names_from = year, values_from = all_of(vars),
                names_glue = paste0("{.value}", suffix, "{year}"))
}

son_wide = son_records |>
  widen(union(son_multiyear, son_1940only), son_years, "_", "son_hik")

pop_wide = father_records |>
  select(son_hik, year, pid, starts_with("dad_")) |>
  rename_with(~ sub("^dad_", "", .x)) |>
  widen(pop_multiyear, pop_years, "_pop_", "son_hik")

aian_clean = sons |>
  inner_join(distinct(father_records, son_hik, pid), by = "son_hik") |>
  left_join(pop_wide, by = "son_hik") |>
  left_join(son_wide, by = "son_hik") |>
  rename(occ_son = occ1950_1940) |>
  select(-son_hik) |>
  select(where(~ !all(is.na(.x))))

cat("\nNon-missing counts for father-side variables by year:\n")
print(aian_clean |>
  summarise(across(matches("_pop_\\d{4}$"), ~ sum(!is.na(.x)))) |>
  pivot_longer(everything(), names_to = "col", values_to = "n") |>
  separate_wider_regex(col, c(var = ".*", "_pop_", yr = "\\d{4}")) |>
  pivot_wider(names_from = yr, values_from = n),
n = 40)

stopifnot(
  "duplicate sons" = !any(duplicated(aian_clean$histid_1940)),
  "pid missing" = !any(is.na(aian_clean$pid)),
  "occ_son missing" = "occ_son" %in% names(aian_clean),
  "no occ1950_pop cols" = any(grepl("^occ1950_pop_\\d{4}$", names(aian_clean))),
  "no age_pop cols" = any(grepl("^age_pop_\\d{4}$", names(aian_clean))))

cat("\nFather-year coverage (drives the single-observation problem):\n")
print(aian_clean |>
        summarise(across(matches("^occ1950_pop_\\d{4}$"), ~ sum(!is.na(.x)))) |>
        pivot_longer(everything(), names_to = "col", values_to = "n"))

n_obs = aian_clean |>
  transmute(k = rowSums(!is.na(pick(matches("^occ1950_pop_\\d{4}$")))))
cat("\nFather observations per son:\n")
print(count(n_obs, k))
cat("Single-observation fathers:", sum(n_obs$k == 1),
    sprintf("(%.1f%%)\n", 100 * mean(n_obs$k == 1)))

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

modal_meso_pop = pick_modal_meso(aian_clean, aian_age, prefer_employed = FALSE, empstatd_tiebreak = FALSE) |>
  rename(meso_pop = meso, picked_year = year)

aian_merged = aian_clean |>
  left_join(modal_meso_pop, by = "pid") |>
  left_join(aian_age |> select(pid, birthyr_spread = spread, spread_mad), by = "pid") |>
  mutate(birthyr_son = 1940 - age_1940) |>
  select(-starts_with("age"), -starts_with("occ1950_pop")) |>
  filter(!is.na(meso_pop)) |>
  (\(x) { cat("After missing meso_pop drop:", nrow(x), "father-son pairs\n"); x })() |>
  filter(is.na(spread_mad) | spread_mad <= 4) |>
  (\(x) { cat("After spread filter (MAD <= 4):", nrow(x), "father-son pairs\n"); x })() |>
  mutate(spread_flag = !is.na(spread_mad) & spread_mad > 2,
         macro_pop = classify_macro(meso_pop),
         meso_son = classify_meso(occ_son),
         macro_son = classify_macro(meso_son)) |>
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
  relocate(meso_pop, .after = picked_year) |>
  relocate(macro_pop, .after = meso_pop) |>
  relocate(occ_pop, .after = macro_pop) |>
  relocate(starts_with("macro_son"), .after = occ_son) |>
  relocate(starts_with("meso_son"), .after = macro_son) |>
  relocate(lit_son, .after = educd_1940) |>
  relocate(lit_pop, .after = educd_pop_1940) |>
  relocate(starts_with("w_parent"), .after = last_col())

cat("\nFinal analysis sample:", nrow(aian_merged), "father-son pairs\n")

cat("meso_pop values:", paste(sort(unique(as.character(aian_merged$meso_pop))), collapse = ", "), "\n")
stopifnot(all(aian_merged$meso_pop %in% meso_order))

saveRDS(aian_merged, "data/aian_merged.rds")
