library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)

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
    group_by({{ level }}) |>
    summarise(total_w = sum(weight),
                         .groups = "drop") |>
    mutate(pi0 = total_w / sum(total_w)) |>
    arrange({{ level }})
  
  pi0_vec = df$pi0
  names(pi0_vec) = df |> pull({{ level }})
  return(pi0_vec)
}

p_matrix = function(data, level_dad, level_son, matrix = TRUE) {
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
  
  ifelse(matrix == TRUE, return(transition_df), return(p_mat))
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

p_mat = p_matrix(data, dad_meso, occ_meso)
g = ggplot(
  p_mat |> 
    mutate(dad_meso = factor(dad_meso, levels = rev(unique(dad_meso)))),
  aes(x = occ_meso, y = dad_meso, fill = P)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", P)),
            size = 3, color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = "Son occupation",
       y = "Father occupation",
       fill = "Transition Prob.",
       title = "P") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

pi0_vec = pi_0(data, dad_meso)
df_pi0 = tibble(
  father = names(pi0_vec),
  pi0 = as.numeric(pi0_vec)) |>
  mutate(father = factor(father, levels = rev(unique(father))))

g0 = ggplot(df_pi0, aes(x = 1, y = father, fill = pi0)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi0)),
            size = 3, color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "π_0") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))


steady = pi_star(p_matrix(data, dad_meso, occ_meso, FALSE))
df_pi_star = tibble(
  father = names(steady),
  pi_star = as.numeric(steady)) |>
  mutate(father = factor(father, levels = rev(unique(father))))

g_star = ggplot(df_pi_star, aes(x = 1, y = father, fill = pi_star)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi_star)),
            size = 3, color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "π*") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

combined_plot = g + g0 + g_star + plot_layout(widths = c(6, 1, 1))



