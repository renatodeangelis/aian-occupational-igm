library(dplyr)
library(readr)
library(tidyr)
library(cobalt)

source("code/utils.R")

aian_merged = read_csv("data/aian_merged.csv") |>
  mutate(linked = 1,
         region = assign_region(statefip_1940),
         education = case_when(
           educd_1940 == 2 ~ "0",
           educd_1940 %in% 14:17 ~ "1-4",
           educd_1940 %in% 22:26 ~ "5-8",
           educd_1940 %in% 30:60 ~ "9-12",
           educd_1940 %in% 70:113 ~ "12+",
           educd_1940 == 999 ~ "missing"))

aian_full = read_csv(
  file = "https://www.dropbox.com/scl/fi/imhz0zujc9dfhx3dc6kox/usa_00021.csv?rlkey=1b5xq14od1cky4iyvfbrob413&st=ndwbat7q&dl=1") |>
  janitor::clean_names() |>
  filter(age >= 20 & age < 45,
         school == 1) |>
  mutate(linked = 0,
         birthyr_son = 1940 - age,
         meso_son = classify_meso(occ1950),
         macro_son = classify_macro(meso_son),
         meso_son_alt = classify_meso(occ1950, split_farmer = FALSE),
         macro_son_alt = classify_macro(meso_son_alt),
         region = assign_region(statefip),
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
    region = as.factor(region),
    education = as.factor(education),
    birthyr_son = as.factor(birthyr_son))

ps_model = glm(linked ~ birthyr_son * region + education * region + statefip_1940,
            data = aian_comb,
            family = binomial)

aian_ps = aian_comb %>%
  mutate(p_hat = predict(ps_model, newdata =., type = "response")) |>
  mutate(w_atc = if_else(linked == 1,
                         (1 - p_hat) / p_hat, NA_real_)) |>
  filter(linked == 1) |>
  mutate(w_atc_norm = w_atc * n() / sum(w_atc)) |>
  select(-linked, -(countyicp:p_hat)) |>
  relocate(starts_with("w_parent"), .after = last_col())

bal.tab(linked ~ birthyr_son + region + education,
        data = aian_comb, weights = "w_atc_norm",
        method = "weighting", estimand = "ATC")

write_csv(aian_ps, "data/aian_weighted.csv")
