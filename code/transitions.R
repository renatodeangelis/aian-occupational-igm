library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tibble)

data = read_csv("aian-igm/data/aian_filtered.csv")

transition_df <- data_occ |>
  count(dad_macro, son_macro, name = "n") |>            # count father→son dyads
  group_by(dad_macro) |>
  mutate(
    n_i      = sum(n),                                  # total dyads for dad_macro = i
    P        = n / n_i,                                 # transition probability P[i,j]
    se       = sqrt(P * (1 - P) / n_i)                  # approx. SE under binomial
  ) |>
  ungroup()

# Transition matrix P
P_mat <- transition_df |>
  select(dad_macro, son_macro, P) |>
  pivot_wider(
    names_from  = son_macro,
    values_from = P) |>
  column_to_rownames("dad_macro")

# Standard‐error matrix for P
SE_mat <- transition_df |>
  select(dad_macro, son_macro, se) |>
  pivot_wider(
    names_from  = son_macro,
    values_from = se) |>
  column_to_rownames("dad_macro")

# Initial distribution π₀
pi0 <- data_occ |>
  count(dad_macro, name = "n") |>
  mutate(pi0 = n / sum(n)) |>
  arrange(dad_macro) |>
  pull(pi0) 
names(pi0) <- levels(data_occ$dad_macro)

# Stationary distribution π*
# Solve π* = π* P  ⇒  π* is left eigenvector of P with eigenvalue 1
eig <- eigen(t(as.matrix(P_mat)))
idx <- which.min(abs(eig$values - 1))                # find eigenvalue ≈ 1
v   <- Re(eig$vectors[, idx])
pi_star <- v / sum(v)                                 # normalize to sum to 1
names(pi_star) <- rownames(P_mat)

