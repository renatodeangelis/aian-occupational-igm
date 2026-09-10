library(dplyr)
library(readr)
library(tidyr)
library(janitor)

source("code/utils.R")

aian = 3
pop_years = c(1910, 1920, 1930, 1940)
son_years = c(1910, 1920, 1930, 1940, 1950)

son_multiyear = c("lit", "gq", "gqtype", "school", "relate", "age",
                  "statefip", "speakeng", "occ1950")
son_1940only  = c("educ", "educd", "sex", "countyicp", "urban", "metro",
                  "empstat", "empstatd", "labforce", "classwkr", "ind1950",
                  "farm", "ownershp", "marst", "bpl", "birthyr",
                  "wkswork1", "hrswork1", "hrswork2", "durunemp",
                  "incwage", "incnonwg", "migrate5", "migrate5d", "migplac5")
pop_multiyear = c("histid", "hik", "age", "birthyr", "occ1950", "ind1950",
                  "classwkr", "labforce", "empstat", "empstatd", "lit",
                  "speakeng", "farm", "ownershp", "gq", "gqtype", "relate",
                  "marst", "bpl", "educd", "statefip", "countyicp", "urban")

path = "https://www.dropbox.com/scl/fi/q4725zk5qw3ltruiucwgc/usa_00026.csv?rlkey=vzauffbwyj61upzhb4fezbq5c&st=y9467v1h&dl=1"

raw = read_csv(path, col_types = cols(.default = col_character())) |>
  clean_names() |>
  type_convert(col_types = cols(histid = col_character(),
                                hik = col_character(),
                              .default = col_guess()))

sons = raw |>
  filter(year == 1940, sex == 1, between(age, 20, 44), race == AIAN,
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
 
saveRDS(aian_clean, "data/aian_clean.rds")
