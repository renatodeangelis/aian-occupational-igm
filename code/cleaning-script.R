library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(sf)

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
  rename(son_occ = occ1950_1940) |>
  filter(age_1940 >= 20 & age_1940 < 40,
         sex_1940 == 1,
         !(sex_pop_1940 == 2))

modal_occ = aian_clean |>
  pivot_longer(
    cols = starts_with("occ1950_pop_"),
    names_to = "year",
    values_to = "occ",
    values_drop_na = TRUE) |>
  mutate(year_num = as.integer(str_extract(year, "\\d{4}"))) |>
  group_by(histid_pop_1940, occ) |>
  summarise(freq = n(),
            first_year = min(year_num),
            .groups = "drop") |>
  group_by(histid_pop_1940) |>
  filter(freq == max(freq)) |>
  slice_min(first_year, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(histid_pop_1940, modal_occ = occ)

aian_filtered = aian_clean |>
  left_join(modal_occ, by = "histid_pop_1940") |>
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

unzip("nhgis0002_shapefile_tl2000_us_county_1940.zip", exdir = "cty_1940_shp")
unzip("cb_2018_us_aiannh_500k.zip", exdir = "res_1920_shp")

cty_40 = st_read("cty_1940_shp/US_county_1940.shp")
res_20 = st_read("res_1920_shp/cb_2018_us_aiannh_500k.shp")

res = st_transform(res_20, 5070)
cty = st_transform(cty_40, 5070) |>
  mutate(state_icp = as.integer(ICPSRST),
         county_icp = as.integer(ICPSRCTY))

overlap = st_intersection(cty, res) |>
  st_drop_geometry() |>
  distinct(state_icp, county_icp) |>
  arrange(state_icp, county_icp) |>
  mutate(res_cty = 1)

aian_res = aian_filtered |>
  mutate(state_icp = case_when(
    statefip_1940 == 1  ~ 41,    # Alabama
    statefip_1940 == 2  ~ 81,   # Alaska
    statefip_1940 == 4  ~ 61,    # Arizona
    statefip_1940 == 5  ~ 42,    # Arkansas
    statefip_1940 == 6  ~ 71,    # California
    statefip_1940 == 8  ~ 62,    # Colorado
    statefip_1940 == 9  ~ 1,    # Connecticut
    statefip_1940 == 10 ~ 11,    # Delaware
    statefip_1940 == 11 ~ 98,    # DC
    statefip_1940 == 12 ~ 43,    # Florida
    statefip_1940 == 13 ~ 44,   # Georgia
    statefip_1940 == 15 ~ 82,   # Hawaii
    statefip_1940 == 16 ~ 63,   # Idaho
    statefip_1940 == 17 ~ 21,   # Illinois
    statefip_1940 == 18 ~ 22,   # Indiana
    statefip_1940 == 19 ~ 31,   # Iowa
    statefip_1940 == 20 ~ 32,   # Kansas
    statefip_1940 == 21 ~ 51,   # Kentucky
    statefip_1940 == 22 ~ 45,   # Louisiana
    statefip_1940 == 23 ~ 2,   # Maine
    statefip_1940 == 24 ~ 52,   # Maryland
    statefip_1940 == 25 ~ 3,   # Massachusetts
    statefip_1940 == 26 ~ 23,   # Michigan
    statefip_1940 == 27 ~ 33,   # Minnesota
    statefip_1940 == 28 ~ 46,   # Mississippi
    statefip_1940 == 29 ~ 34,   # Missouri
    statefip_1940 == 30 ~ 64,   # Montana
    statefip_1940 == 31 ~ 35,   # Nebraska
    statefip_1940 == 32 ~ 65,   # Nevada
    statefip_1940 == 33 ~ 4,   # New Hampshire
    statefip_1940 == 34 ~ 12,   # New Jersey
    statefip_1940 == 35 ~ 66,   # New Mexico
    statefip_1940 == 36 ~ 13,   # New York
    statefip_1940 == 37 ~ 47,   # North Carolina
    statefip_1940 == 38 ~ 36,   # North Dakota
    statefip_1940 == 39 ~ 24,   # Ohio
    statefip_1940 == 40 ~ 53,   # Oklahoma
    statefip_1940 == 41 ~ 72,   # Oregon
    statefip_1940 == 42 ~ 14,   # Pennsylvania
    statefip_1940 == 44 ~ 5,   # Rhode Island
    statefip_1940 == 45 ~ 48,   # South Carolina
    statefip_1940 == 46 ~ 37,   # South Dakota
    statefip_1940 == 47 ~ 54,   # Tennessee
    statefip_1940 == 48 ~ 49,   # Texas
    statefip_1940 == 49 ~ 67,   # Utah
    statefip_1940 == 50 ~ 6,   # Vermont
    statefip_1940 == 51 ~ 40,   # Virginia
    statefip_1940 == 53 ~ 73,   # Washington
    statefip_1940 == 54 ~ 56,   # West Virginia
    statefip_1940 == 55 ~ 25,   # Wisconsin
    statefip_1940 == 56 ~ 68,   # Wyoming
    TRUE ~ NA_real_)) |>
  rename("county_icp" = "countyicp_1940") |>
  left_join(overlap, by = c("state_icp", "county_icp")) |>
  mutate(res_cty = replace_na(res_cty, 0))

setwd("~/")

write_csv(aian_res, "aian-igm/data/aian_filtered.csv")
write_csv(overlap, "aian-igm/data/res_counties.csv")
