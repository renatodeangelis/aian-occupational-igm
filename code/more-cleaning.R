library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

data = read_csv("data/aian_filtered.csv")

modal_occ = data |>
  pivot_longer(
    cols = starts_with("modal_occ_pop_"),
    names_to = "year",
    values_to = "occ",
    values_drop_na = TRUE) |>
  group_by(histid_pop_1940, occ) |>
  summarise(freq = n(), .groups = "drop") |>
  group_by(histid_pop_1940) |>
  filter(freq == max(freq)) |>
  left_join(
    data |>
      pivot_longer(
        cols = starts_with("modal_occ_pop_"),
        names_to = "year",
        values_to = "occ",
        values_drop_na = TRUE) |>
      mutate(year_num = str_extract(year, "\\d{4}") |> as.integer()),
    by = c("histid_pop_1940", "occ")) |>
  group_by(histid_pop_1940) |>
  slice_min(year_num, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(histid_pop_1940, modal_occ = occ)
           
data_occ = data |>
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

transition_df <- data_occ |>
  count(dad_macro, son_macro, name = "n") |>            # count father→son dyads
  group_by(dad_macro) |>
  mutate(
    n_i      = sum(n),                                  # total dyads for dad_macro = i
    P        = n / n_i,                                 # transition probability P[i,j]
    se       = sqrt(P * (1 - P) / n_i)                  # approx. SE under binomial
  ) |>
  ungroup()

#─── 3. Pivot into matrices ─────────────────────────────────────────────────
# Transition matrix P
P_mat <- transition_df |>
  select(dad_macro, son_macro, P) |>
  pivot_wider(
    names_from  = son_macro,
    values_from = P
  ) |>
  column_to_rownames("dad_macro")

# Standard‐error matrix for P
SE_mat <- transition_df |>
  select(dad_macro, son_macro, se) |>
  pivot_wider(
    names_from  = son_macro,
    values_from = se
  ) |>
  column_to_rownames("dad_macro")

#─── 4. Initial distribution π₀ ──────────────────────────────────────────────
pi0 <- data_occ |>
  count(dad_macro, name = "n") |>
  mutate(pi0 = n / sum(n)) |>
  arrange(dad_macro) |>
  pull(pi0) 
names(pi0) <- levels(data_occ$dad_macro)

#─── 5. Stationary distribution π* ─────────────────────────────────────────
# Solve π* = π* P  ⇒  π* is left eigenvector of P with eigenvalue 1
eig <- eigen(t(as.matrix(P_mat)))
idx <- which.min(abs(eig$values - 1))                # find eigenvalue ≈ 1
v   <- Re(eig$vectors[, idx])
pi_star <- v / sum(v)                                 # normalize to sum to 1
names(pi_star) <- rownames(P_mat)

P_mat
SE_mat
pi0
pi_star

