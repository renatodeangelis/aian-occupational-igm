library(dplyr)
library(readr)
library(janitor)

setwd("~/Downloads")
aian_raw = readr::read_csv(
  file = "clp_mlp1850_1940_linked_subsample_300raced_2022-11-5.csv",
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
  select_if(~ !all(is.na(.))) |>
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
  rename(son_occ = occ1950_1940)

aian_filtered = aian_clean |>
  filter(age_1940 >= 15 & age_1940 <= 44,
         sex_1940 == 1,
         sex_pop_1940 == 1 | is.na(sex_pop_1940)) |>
  filter((!is.na(occ1950_pop_1900) & occ1950_pop_1900 <= 970) | 
           (!is.na(occ1950_pop_1910) & occ1950_pop_1910 <= 970) | 
           (!is.na(occ1950_pop_1920) & occ1950_pop_1920 <= 970) | 
           (!is.na(occ1950_pop_1930) & occ1950_pop_1930 <= 970) | 
           (!is.na(occ1950_pop_1940) & occ1950_pop_1940 <= 970))

write_csv(aian_filtered, "code/aian_filtered.csv")