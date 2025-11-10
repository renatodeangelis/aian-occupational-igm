library(dplyr)
library(readr)
library(tidyr)

aian_merged = read_csv("data/aian_merged.csv") |>
  mutate(linked = 1,
         region = case_when(
           statefip_1940 %in% c(8, 16, 32, 49, 56)  ~ "basin",
           statefip_1940 == 6                         ~ "cali",
           statefip_1940 %in% c(27, 55)               ~ "lakes",
           statefip_1940 %in% c(17, 18, 26, 39)       ~ "midwest",
           statefip_1940 == 37                        ~ "nc",
           statefip_1940 %in% c(9, 10, 23, 24, 25, 33, 34, 36, 42, 44, 50, 11) ~ "ne",
           statefip_1940 %in% c(30, 38, 46)           ~ "plains",
           statefip_1940 %in% c(41, 53)               ~ "nw",
           statefip_1940 == 40                        ~ "ok",
           statefip_1940 %in% c(19, 20, 29, 31)       ~ "prairie",
           statefip_1940 %in% c(1, 5, 12, 13, 21, 22, 28, 45, 47, 48, 51, 54) ~ "south",
           statefip_1940 %in% c(4, 35)                ~ "sw",
           TRUE ~ NA_character_),
         education = case_when(
           educd_1940 == 2 ~ "0",
           educd_1940 %in% 14:17 ~ "1-4",
           educd_1940 %in% 22:26 ~ "5-8",
           educd_1940 %in% 30:60 ~ "9-12",
           educd_1940 %in% 70:113 ~ "12+",
           educd_1940 == 999 ~ "missing"))

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

aian_full = read_csv(
  file = "https://www.dropbox.com/scl/fi/imhz0zujc9dfhx3dc6kox/usa_00021.csv?rlkey=1b5xq14od1cky4iyvfbrob413&st=ndwbat7q&dl=1") |>
  janitor::clean_names() |>
  filter(age < 45,
         school == 1) |>
  mutate(linked = 0,
         birthyr_son = 1940 - age,
         meso_son = classify_meso(occ1950),
         macro_son = classify_macro(meso_son),
         meso_son_alt = classify_meso(occ1950, split_farmer = FALSE),
         macro_son_alt = classify_macro(meso_son_alt),
         region = case_when(
           statefip %in% c(8, 16, 32, 49, 56)  ~ "basin",
           statefip == 6                         ~ "cali",
           statefip %in% c(27, 55)               ~ "lakes",
           statefip %in% c(17, 18, 26, 39)       ~ "midwest",
           statefip == 37                        ~ "nc",
           statefip %in% c(9, 10, 23, 24, 25, 33, 34, 36, 42, 44, 50, 11) ~ "ne",
           statefip %in% c(30, 38, 46)           ~ "plains",
           statefip %in% c(41, 53)               ~ "nw",
           statefip == 40                        ~ "ok",
           statefip %in% c(19, 20, 29, 31)       ~ "prairie",
           statefip %in% c(1, 5, 12, 13, 21, 22, 28, 45, 47, 48, 51, 54) ~ "south",
           statefip %in% c(4, 35)                ~ "sw",
           TRUE ~ NA_character_),
         education = case_when(
           educd == 2 ~ "0",
           educd %in% 14:17 ~ "1-4",
           educd %in% 22:26 ~ "5-8",
           educd %in% 30:60 ~ "9-12",
           educd %in% 70:113 ~ "12+",
           educd == 999 ~ "missing")) |>
  rename(statefip_1940 = statefip)

aian_comb = bind_rows(aian_merged, aian_full) |>
  mutate(
    meso_son = as.factor(meso_son),
    meso_son_alt = as.factor(meso_son_alt),
    region = as.factor(region),
    education = as.factor(education),
    birthyr_son = as.factor(birthyr_son))

ps_model = glm(linked ~ birthyr_son + meso_son + region + education,
            data = aian_comb,
            family = binomial)

aian_ps = aian_comb %>%
  mutate(p_hat = predict(ps_model, newdata =., type = "response")) |>
  mutate(w_atc = if_else(linked == 1,
                         (1 - p_hat) / p_hat, NA_real_)) |>
  filter(linked == 1) |>
  mutate(w_atc_norm = w_atc * n() / sum(w_atc))
  select(-linked, -(countyicp:p_hat)) |>
  relocate(starts_with("w_parent"), .after = last_col())

write_csv(aian_ps, "data/aian_weighted.csv")













