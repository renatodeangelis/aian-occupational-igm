library(dplyr)
library(readr)
library(tidyr)
library(cobalt)

aian_filtered = read_csv("data/aian_merged.csv")
res_counties = read_csv("data/res_counties.csv")

setwd("~/Downloads")
aian_full = read_csv("usa_00020.csv") |>
  janitor::clean_names() |>
  filter(age < 45) |>
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
         county_icp = countyicp,
         educd_1940 = educd)

aian_filtered = aian_filtered |>
  mutate(linked = 1) |>
  rename(occ_meso = son_meso,
         occ_macro = son_macro,
#         stateicp = state_icp,
         age = age_1940) #|>
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
    age_group = cut(age, breaks = c(20, 25, 30, 35, 40, 45), right = FALSE),
    res_cty = as.factor(res_cty),
    occ_meso = as.factor(occ_meso),
    region = as.factor(region))

sum_linked = sum(aian_ps$linked == 1)

ps_model = glm(linked ~ age_group + occ_meso,
            data = aian_comb,
            family = binomial)

aian_ps = aian_comb %>%
  mutate(p_hat = predict(ps_model, newdata =., type = "response")) |>
  mutate(w_atc = if_else(linked == 1,
                         (1 - p_hat) / p_hat, NA_real_)) |>
  mutate(w_atc_norm = if_else(linked == 1,
                              w_atc * (sum_linked / sum(w_atc[linked == 1])),
                              NA_real_))

comb_for_bal = aian_ps |>
  mutate(w_for_nal = if_else(is.na(w_atc_norm), 1, w_atc_norm))

bal.tab(linked ~ age_group + occ_meso + region + res_cty,
        data = comb_for_bal,
        weights = comb_for_bal$w_for_nal,
        estimand = "ATC",
        un = TRUE)

aian_weighted = aian_ps |>
  filter(linked == 1) |>
  select(-linked, -(year:histid), -p_hat) |>
  rename(weight = w_atc)

setwd("~/aian-igm/data")
write_csv(aian_weighted, "aian_weighted.csv")














