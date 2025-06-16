library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(expm)
library(rlang)

macro_order = c("farming", "blue_col", "white_col", "unemp")
meso_order = c("farming", "unskilled", "crafts", "clerical", "prof", "unemp")
data = read_csv("data/aian_weighted.csv") |>
  mutate(
    dad_macro = factor(dad_macro, levels = macro_order),
    occ_macro = factor(occ_macro, levels = macro_order),
    dad_meso = factor(dad_meso, levels = meso_order),
    occ_meso = factor(occ_meso, levels = meso_order))

################################################################################
############################### BASIC MEASURES #################################
################################################################################

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
  dad_nm = as_name(ensym(level_dad))
  son_nm = as_name(ensym(level_son))
  
  df = data |>
    group_by({{ level_dad }}, {{ level_son }}) |>
    summarise(n_w = sum(weight), .groups = "drop") |>
    group_by({{ level_dad }}) |>
    mutate(
      n_i_w = sum(n_w),
      P = n_w / n_i_w) |>
    ungroup()
  
  wide = df |>
    select(all_of(c(dad_nm, son_nm, "P"))) |>
    pivot_wider(
      names_from = all_of(son_nm),
      values_from = P)
  
  mat = as.matrix(wide[-1])
  rownames(mat) = wide[[dad_nm]]
  
  ifelse(matrix == TRUE, return(mat), return(df))
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


################################################################################
############################## ADVANCED MEASURES ###############################
################################################################################

tv_norm = function(mu, nu){
  0.5 * sum(abs(mu - nu))
}

d_t = function(data, level_dad, level_son, t = 1){
  P_mat = p_matrix(data, {{ level_dad }}, {{ level_son }})
  pi_star = pi_star(P_mat)
  
  P_t = P_mat %^% t
  
  d_i = apply(P_t, 1, function(row_i) tv_norm(row_i, pi_star))
  
  log(max(d_i))
}

d_prime = function(data, level_dad, level_son, t = 1){
  P_mat = p_matrix(data, {{ level_dad }}, {{ level_son }})
  
  P_t = P_mat %^% t
  
  n = nrow(P_t)
  
  pairs = combn(n, 2)
  dvals = apply(pairs, 2, function(idx) {
    i = idx[1]; j = idx[2]
    tv_norm(P_t[i, ], P_t[j, ])
  })
  
  log(max(dvals))
}

am = function(data, level_dad, level_son, t = 1){
  pi_0 = pi_0(data, {{ level_dad }})
  P_mat = p_matrix(data, {{ level_dad }}, {{ level_son }})
  pi_star = pi_star(P_mat)
  
  P_t = P_mat %^% t
  
  pi_t = as.numeric(pi_0 %*% P_t)
  
  am = tv_norm(pi_t, pi_star)
  
  log(am)
}

im = function(data, level_dad, level_son, t = 1){
  P_mat = p_matrix(data, {{ level_dad }}, {{ level_son }})
  pi_star = pi_star(P_mat)
  
  P_t = P_mat %^% t
  
  im_i = apply(P_t, 1, function(row_i) tv_norm(row_i, pi_star))
  
  log(im_i)
}

aim = function(data, level_dad, level_son, t = 1){
  imv = im(data, {{ level_dad }}, {{ level_son}}, t = t)
  pi0 = pi_0(data, {{ level_dad }})
  sum(imv * pi0[names(imv)])
}

sdm_core = function(pi0, P1, mu0, P2 = P1, t = 1){
  P_t_pi = P1 %^% t
  P_t_mu = P2 %^% t
  
  pi0_t = as.numeric(pi0 %*% P_t_pi)
  mu0_t = as.numeric(mu0 %*% P_t_mu)
  
  tv_norm(pi0_t, mu0_t)
}

sdm_data = function(data, level_dad, level_son, subgroup_var, g1, g2, t = 1){
  P_mat = p_matrix(data, {{ level_dad }}, {{ level_son}}, matrix = FALSE)
  
  pi0 = data |>
    filter({{ subgroup_var }} == g1) |>
    pi_0({{ level_dad }})
  
  mu0 = data |>
    filter({{ subgroup_var }} == g2) |>
    pi_0({{ level_dad }})
  
  sdm_core(pi0, P_mat, mu0, P_mat, t = t)
}

################################################################################
############################### IMPLEMENTATION #################################
################################################################################

## OVERALL TRANSITIONS (MACRO)

p_mat = p_matrix(data, dad_macro, occ_macro, matrix = FALSE)
g = ggplot(
  p_mat |> 
    mutate(dad_macro = factor(dad_macro, levels = rev(unique(dad_macro)))),
  aes(x = occ_macro, y = dad_macro, fill = P)) +
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

pi0_vec = pi_0(data, dad_macro)
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

steady = pi_star(p_matrix(data, dad_macro, occ_macro, TRUE))
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

## RESERVATION ONLY

res = data |> filter(res_cty == 1 & statefip_1940 != 40)

p_mat_res = p_matrix(res, dad_macro, occ_macro, matrix = FALSE)
g_res = ggplot(
  p_mat_res |> 
    mutate(dad_macro = factor(dad_macro, levels = rev(unique(dad_macro)))),
  aes(x = occ_macro, y = dad_macro, fill = P)) +
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

pi0_vec_res = pi_0(res, dad_macro)
df_pi0_res = tibble(
  father = names(pi0_vec_res),
  pi0 = as.numeric(pi0_vec_res)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g0_res = ggplot(df_pi0_res, aes(x = 1, y = father, fill = pi0)) +
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

steady_res = pi_star(p_matrix(res, dad_macro, occ_macro, TRUE))
df_pi_star_res = tibble(
  father = names(steady_res),
  pi_star = as.numeric(steady_res)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g_star_res = ggplot(df_pi_star_res, aes(x = 1, y = father, fill = pi_star)) +
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

combined_plot_res = g_res + g0_res + g_star_res + plot_layout(widths = c(6, 1, 1))

## NONRESERVATION ONLY

nonres = data |> filter(res_cty == 0 | statefip_1940 == 40)

p_mat_nonres = p_matrix(nonres, dad_macro, occ_macro, matrix = FALSE)
g_nonres = ggplot(
  p_mat_nonres |> 
    mutate(dad_macro = factor(dad_macro, levels = rev(unique(dad_macro)))),
  aes(x = occ_macro, y = dad_macro, fill = P)) +
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

pi0_vec_nonres = pi_0(nonres, dad_macro)
df_pi0_nonres = tibble(
  father = names(pi0_vec_nonres),
  pi0 = as.numeric(pi0_vec_nonres)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g0_nonres = ggplot(df_pi0_nonres, aes(x = 1, y = father, fill = pi0)) +
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

steady_nonres = pi_star(p_matrix(nonres, dad_macro, occ_macro, TRUE))
df_pi_star_nonres = tibble(
  father = names(steady_nonres),
  pi_star = as.numeric(steady_nonres)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g_star_nonres = ggplot(df_pi_star_nonres, aes(x = 1, y = father, fill = pi_star)) +
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

combined_plot_nonres = g_nonres + g0_nonres + g_star_nonres + 
  plot_layout(widths = c(6, 1, 1))

## OKLAHOMA ONLY

## NONRESERVATION ONLY

okl = data |> filter(statefip_1940 == 40)

p_mat_okl = p_matrix(okl, dad_macro, occ_macro, matrix = FALSE)
g_okl = ggplot(
  p_mat_okl |> 
    mutate(dad_macro = factor(dad_macro, levels = rev(unique(dad_macro)))),
  aes(x = occ_macro, y = dad_macro, fill = P)) +
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

pi0_vec_okl = pi_0(okl, dad_macro)
df_pi0_okl = tibble(
  father = names(pi0_vec_okl),
  pi0 = as.numeric(pi0_vec_okl)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g0_okl = ggplot(df_pi0_okl, aes(x = 1, y = father, fill = pi0)) +
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

steady_okl = pi_star(p_matrix(okl, dad_macro, occ_macro, TRUE))
df_pi_star_okl = tibble(
  father = names(steady_okl),
  pi_star = as.numeric(steady_okl)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g_star_okl = ggplot(df_pi_star_okl, aes(x = 1, y = father, fill = pi_star)) +
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

combined_plot_okl = g_okl + g0_okl + g_star_okl + 
  plot_layout(widths = c(6, 1, 1))

## OVERALL TRANSITIONS (MESO)

p_mat = p_matrix(data, dad_meso, occ_meso, matrix = FALSE)
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

steady = pi_star(p_matrix(data, dad_meso, occ_meso, TRUE))
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

## RESERVATION ONLY

res = data |> filter(res_cty == 1 & statefip_1940 != 40)

p_mat_res = p_matrix(res, dad_meso, occ_meso, matrix = FALSE)
g_res = ggplot(
  p_mat_res |> 
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

pi0_vec_res = pi_0(res, dad_meso)
df_pi0_res = tibble(
  father = names(pi0_vec_res),
  pi0 = as.numeric(pi0_vec_res)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g0_res = ggplot(df_pi0_res, aes(x = 1, y = father, fill = pi0)) +
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

steady_res = pi_star(p_matrix(res, dad_meso, occ_meso, TRUE))
df_pi_star_res = tibble(
  father = names(steady_res),
  pi_star = as.numeric(steady_res)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g_star_res = ggplot(df_pi_star_res, aes(x = 1, y = father, fill = pi_star)) +
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

combined_plot_res = g_res + g0_res + g_star_res + plot_layout(widths = c(6, 1, 1))

## NONRESERVATION ONLY

nonres = data |> filter(res_cty == 0 | statefip_1940 == 40)

p_mat_nonres = p_matrix(nonres, dad_meso, occ_meso, matrix = FALSE)
g_nonres = ggplot(
  p_mat_nonres |> 
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

pi0_vec_nonres = pi_0(nonres, dad_meso)
df_pi0_nonres = tibble(
  father = names(pi0_vec_nonres),
  pi0 = as.numeric(pi0_vec_nonres)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g0_nonres = ggplot(df_pi0_nonres, aes(x = 1, y = father, fill = pi0)) +
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

steady_nonres = pi_star(p_matrix(nonres, dad_meso, occ_meso, TRUE))
df_pi_star_nonres = tibble(
  father = names(steady_nonres),
  pi_star = as.numeric(steady_nonres)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g_star_nonres = ggplot(df_pi_star_nonres, aes(x = 1, y = father, fill = pi_star)) +
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

combined_plot_nonres = g_nonres + g0_nonres + g_star_nonres + 
  plot_layout(widths = c(6, 1, 1))

## OKLAHOMA ONLY

okl = data |> filter(statefip_1940 == 40)

p_mat_okl = p_matrix(okl, dad_meso, occ_meso, matrix = FALSE)
g_okl = ggplot(
  p_mat_okl |> 
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

pi0_vec_okl = pi_0(okl, dad_meso)
df_pi0_okl = tibble(
  father = names(pi0_vec_okl),
  pi0 = as.numeric(pi0_vec_okl)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g0_okl = ggplot(df_pi0_okl, aes(x = 1, y = father, fill = pi0)) +
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

steady_okl = pi_star(p_matrix(okl, dad_meso, occ_meso, TRUE))
df_pi_star_okl = tibble(
  father = names(steady_okl),
  pi_star = as.numeric(steady_okl)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g_star_okl = ggplot(df_pi_star_okl, aes(x = 1, y = father, fill = pi_star)) +
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

combined_plot_okl = g_okl + g0_okl + g_star_okl + 
  plot_layout(widths = c(6, 1, 1))

################################################################################

## d(t), d'(t), AIM CURVES



## SDM CURVES

