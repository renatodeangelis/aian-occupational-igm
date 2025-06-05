library(dplyr)
library(readr)
library(tidyr)
library(tibble)

macro_order = c("farming", "blue_col", "white_col", "unemp")
meso_order = c("farming", "unskilled", "crafts", "clerical", "prof", "unemp")
data = read_csv("aian_weighted.csv") |>
  mutate(
    dad_macro = factor(dad_macro, levels = macro_order),
    occ_macro = factor(occ_macro, levels = macro_order),
    dad_meso = factor(dad_meso, levels = meso_order),
    occ_meso = factor(occ_meso, levels = meso_order))

pi_0 = function(data, level) {
  df = data |> 
    group_by(across(all_of(level))) |>
    summarise(total_w = sum(weight),
                         .groups = "drop") |>
    mutate(pi0 = total_w / sum(total_w)) |>
    arrange(across(all_of(level)))
  
  pi0_vec = df$pi0
  names(pi0_vec) = df[[level]]
  return(pi0_vec)
}

p_matrix = function(data, level_dad, level_son) {
  transition_df = data |>
    group_by({{ level_dad }}, {{ level_son }}) |>
    summarise(
      n_w = sum(weight),
      .groups = "drop"
    ) |>
    group_by({{ level_dad }}) |>
    mutate(
      n_i_w = sum(n_w),
      P     = n_w / n_i_w
    ) |>
    ungroup()
  
  p_mat = transition_df |>
    select({{ level_dad }}, {{ level_son }}, P) |>
    pivot_wider(
      names_from  = {{ level_son }},
      values_from = P
    ) |>
    column_to_rownames(var = rlang::as_string(rlang::enexpr(level_dad)))
  
  as.matrix(p_mat)
}

pi_star = function(p_mat) {
  P = as.matrix(p_mat)
  eig = eigen(t(P))
  idx = which.min(abs(eig$values - 1))
  v = Re(eig$vectors[, idx])
  pi_star = v / sum(v)
  names(pi_star) = rownames(P)
  return(pi_star)
}




