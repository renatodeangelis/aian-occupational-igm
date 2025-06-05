library(dplyr)
library(readr)
library(tidyr)
library(cobalt)

aian_filtered = read_csv("data/aian_filtered.csv")
res_counties = read_csv("data/res_counties.csv")

setwd("~/Downloads")
aian_full = read_csv("usa_00018.csv") |>
  janitor::clean_names() |>
  mutate(linked = 0) |>
  mutate(occ_meso = case_when(
    occ1950 == 100 | occ1950 == 123 ~ "farming",
    occ1950 <= 99 | (occ1950 >= 200 & occ1950 <= 290) ~ "prof",
    occ1950 >= 300 & occ1950 <= 490 ~ "clerical",
    (occ1950 >= 500 & occ1950 <= 594) | occ1950 == 762 |
      occ1950 == 773 | occ1950 == 781 | occ1950 == 782 ~ "crafts",
    occ1950 >= 595 & occ1950 <= 970 & occ1950 != 762 &
      occ1950 != 773 & occ1950 != 781 & occ1950 != 782 ~ "unskilled",
    occ1950 > 970 ~ "unemp")) |>
  mutate(occ_macro = case_when(
    occ_meso == "farming" ~ "farming",
    occ_meso == "prof" | occ_meso == "clerical" ~ "white_col",
    occ_meso == "crafts" | occ_meso == "unskilled" ~ "blue_col",
    occ_meso == "unemp" ~ "unemp")) |>
  mutate(region = case_when(
    stateicp %in% c(62, 63, 65, 67, 68) ~ "basin",
    stateicp == 71 ~ "cali",
    stateicp %in% c(33, 25) ~ "lakes",
    stateicp %in% c(21, 22, 23, 24) ~ "midwest",
    stateicp == 47 ~ "nc",
    stateicp %in% c(1, 11, 2, 52, 3, 4, 12, 13, 14, 5, 6, 98) ~ "ne",
    stateicp %in% c(64, 36, 37) ~ "plains",
    stateicp %in% c(72, 73) ~ "nw",
    stateicp == 53 ~ "ok",
    stateicp %in% c(31, 32, 34, 35) ~ "prairie",
    stateicp %in% c(41, 42, 43, 44, 51, 45, 46, 48, 54, 49, 40, 56) ~ "south",
    stateicp %in% c(61, 66) ~ "sw",
    TRUE ~ NA_character_)) |>
  rename(state_icp = stateicp,
         county_icp = countyicp) |>
  left_join(res_counties, by = c("state_icp", "county_icp")) |>
  mutate(res_cty = replace_na(res_cty, 0))

aian_filtered = aian_filtered |>
  mutate(linked = 1) |>
  rename(occ_meso = son_meso,
         occ_macro = son_macro,
         stateicp = state_icp,
         age = age_1940) |>
  mutate(region = case_when(
    stateicp %in% c(62, 63, 65, 67, 68) ~ "basin",
    stateicp == 71 ~ "cali",
    stateicp %in% c(33, 25) ~ "lakes",
    stateicp %in% c(21, 22, 23, 24) ~ "midwest",
    stateicp == 47 ~ "nc",
    stateicp %in% c(1, 11, 2, 52, 3, 4, 12, 13, 14, 5, 6, 98) ~ "ne",
    stateicp %in% c(64, 36, 37) ~ "plains",
    stateicp %in% c(72, 73) ~ "nw",
    stateicp == 53 ~ "ok",
    stateicp %in% c(31, 32, 34, 35) ~ "prairie",
    stateicp %in% c(41, 42, 43, 44, 51, 45, 46, 48, 54, 49, 40, 56) ~ "south",
    stateicp %in% c(61, 66) ~ "sw",
    TRUE ~ NA_character_))

aian_comb = bind_rows(aian_filtered, aian_full) |>
  mutate(
    age = as.factor(age),
    res_cty = as.factor(res_cty),
    occ_meso = as.factor(occ_meso),
    region = as.factor(region)
  )

model = glm(linked ~ age + occ_meso + region + res_cty,
            data = aian_comb,
            family = binomial)

matched_weighted = aian_comb |>
  mutate(ps = predict(model, type = "response")) |>
  filter(linked == 1) |>
  mutate(ipw = 1 / ps,
         ipw_norm = ipw / mean(ipw))

bal.tab(linked ~ age + occ_meso + region + res_cty,
        data = aian_comb,
        weights = matched_weighted$ipw_norm,
        estimate = "ATE")

write_csv(aian_weighted, "aian-igm/data/aian_weighted.csv")