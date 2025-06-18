library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(expm)
library(rlang)
library(purrr)
library(sf)
library(tigris)

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
  dad_sym = ensym(level_dad)
  son_sym = ensym(level_son)
  
  dad_nm  = as_string(dad_sym)
  son_nm  = as_string(son_sym)
  
  df = data |>
    group_by( !!dad_sym, !!son_sym ) |>
    summarise(
      n_w = sum(weight),
      .groups = "drop"
    ) |>
    group_by( !!dad_sym ) |>
    mutate(
      n_i_w = sum(n_w),
      P     = ifelse(n_i_w > 0, n_w / n_i_w, 0) 
    ) |>
    ungroup()
  
  if (!matrix) return(df)
  
  wide = df |>
    select( !!dad_sym, !!son_sym, P ) |>
    pivot_wider(
      names_from  = !!son_sym,
      values_from = P,
      values_fill  = list(P = 0)
    )
  
  mat = as.matrix(wide |> select(-!!dad_sym))
  rownames(mat) = pull(wide, !!dad_sym)
  
  mat
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

with_boot = function(data, stat_fn, R = 200, ...){
  args_q = enquos(...)
  N = nrow(data)
  
  call_stat = function(d){
    call_expr = call2(stat_fn, data = d, !!!args_q)
    eval_tidy(call_expr)
  }
  
  est = call_stat(data)
  
  boots = replicate(R, {
    d_b = data |> slice_sample(n = N, replace = TRUE)
    call_stat(d_b)
  }, simplify = FALSE)
  
  boots_arr = simplify2array(boots)
  
  if (is.atomic(est) && length(est) == 1) {
    se = sd(boots_arr, na.rm = TRUE)
  } else if (is.matrix(boots_arr) || (is.array(boots_arr) && length(dim(boots_arr)) == 2)){
    se = apply(boots_arr, 1, sd, na.rm = TRUE)
  } else {
    stop("with_boot(): unsupported output type from stat_fn()")
  }
  
  list(estimate = est,
       se = se)
}

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
  
  log(tv_norm(pi0_t, mu0_t))
}

sdm_mobility = function(data,
                         level_dad, level_son,
                         subgroup_var,
                         res_val, nonres_val,
                         t = 1,
                         variant = c("paired","full","resonly","nonresonly")) {
  variant = match.arg(variant)
  
  P_full = p_matrix(data, {{level_dad}}, {{level_son}}, matrix=TRUE)
  
  d_res = filter(data, {{subgroup_var}} == res_val)
  P_res = p_matrix(d_res, {{level_dad}}, {{level_son}}, matrix=TRUE)
  pi0_res = pi_0(d_res, {{level_dad}})
  
  d_nonres = filter(data, {{subgroup_var}} == nonres_val)
  P_nonres = p_matrix(d_nonres, {{level_dad}}, {{level_son}}, matrix=TRUE)
  pi0_nonres = pi_0(d_nonres,   {{level_dad}})
  
  switch(variant,
         paired = sdm_core(pi0_res, P_res, pi0_nonres, P_nonres, t),
         full = sdm_core(pi0_res, P_full, pi0_nonres, P_full, t),
         resonly = sdm_core(pi0_res, P_res, pi0_nonres, P_res, t),
         nonresonly = sdm_core(pi0_res, P_nonres, pi0_nonres, P_nonres, t)
  )
}

id_within = function(data,
                     level_dad, level_son,
                     subgroup_var, group,
                     j, k,
                     t = 1) {
  sub = data |>
    filter({{ subgroup_var }} == group)
  
  P = p_matrix(sub, {{ level_dad }}, {{ level_son }}, matrix = TRUE)
  
  P_t = if (t == 1) P else P %^% t
  
  row_j = P_t[j, ]
  row_k = P_t[k, ]
  
  tv_norm(row_j, row_k)
}

id_between = function(data,
                      level_dad, level_son,
                      subgroup_var, g1, g2,
                      j,
                      t = 1) {
  P1 = data |>
    filter({{ subgroup_var }} == g1) |>
    p_matrix({{ level_dad }}, {{ level_son }}, matrix = TRUE)
  
  P2 = data |>
    filter({{ subgroup_var }} == g2) |>
    p_matrix({{ level_dad }}, {{ level_son }}, matrix = TRUE)
  
  P1_t = if (t == 1) P1 else P1 %^% t
  P2_t = if (t == 1) P2 else P2 %^% t
  
  tv_norm(P1_t[j, ], P2_t[j, ])
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

data = data |> mutate(res_cty_ok = ifelse(res_cty == 0 | statefip_1940 == 40, 0, 1))

## d(t), d'(t), AIM CURVES

ts = 0:5
measures = tibble(
  measure = c("d", "d_prime", "AM"),
  fn = list(d_t, d_prime, am)
)

results = expand_grid(
  t = ts,
  measures) |>
  mutate(boot = pmap(list(fn, t),
                     ~ with_boot(
                       data,
                       ..1,
                       R = 1000,
                       level_dad = dad_macro,
                       level_son = occ_macro,
                       t = ..2))) |>
  mutate(est = map_dbl(boot, "estimate"),
         se = map_dbl(boot, "se")) |>
  select(measure, t, est, se) |>
  mutate(dt_t = est*t)

write_csv(results, "data/advanced_measures.csv")

ggplot(results, aes(x = t, y = est, color = measure, fill = measure)) +
  geom_ribbon(aes(ymin = est - 1.96*se,
                  ymax = est + 1.96*se),
              alpha = 0.2,
              linetype = 0) +
  geom_line(linetype = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("d" = "black",
                                "d_prime" = "orange",
                                "AM" = "blue"),
                     labels = c(d = expression(log(d(t))),
                                d_prime = expression(log(d * "'"~(t))),
                                AM = expression(log(AM(pi[0], t))))) +
  scale_fill_manual(values = c("d" = "black",
                               "d_prime" = "orange",
                               "AM" = "blue"),
                    labels = c(d = expression(log(d(t))),
                               d_prime = expression(log(d * "'"~(t))),
                               AM = expression(log(AM(pi[0], t))))) +
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8)) +
  labs(x = "Generation (t)",
       y = "log of Measure Value") +
  theme_minimal() +
  theme(legend.position = c(0.4, 0.05),
        legend.justification = c("right", "bottom"),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5))


results_res = expand_grid(
  t = ts,
  measures) |>
  mutate(boot = pmap(list(fn, t),
                     ~ with_boot(
                       res,
                       ..1,
                       R = 1000,
                       level_dad = dad_macro,
                       level_son = occ_macro,
                       t = ..2))) |>
  mutate(est = map_dbl(boot, "estimate"),
         se = map_dbl(boot, "se")) |>
  select(measure, t, est, se) |>
  mutate(dt_t = est*t)

write_csv(results_res, "data/advanced_measures_res.csv")

results_res = read_csv("data/advanced_measures_res.csv")

plot_res = ggplot(results_res, aes(x = t, y = est, color = measure, fill = measure)) +
  geom_ribbon(aes(ymin = est - 1.96*se,
                  ymax = est + 1.96*se),
              alpha = 0.2,
              linetype = 0) +
  geom_line(linetype = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("d" = "black",
                                "d_prime" = "orange",
                                "AM" = "blue"),
                     labels = c(d = expression(log(d(t))),
                                d_prime = expression(log(d * "'"~(t))),
                                AM = expression(log(AM(pi[0], t))))) +
  scale_fill_manual(values = c("d" = "black",
                               "d_prime" = "orange",
                               "AM" = "blue"),
                    labels = c(d = expression(log(d(t))),
                               d_prime = expression(log(d * "'"~(t))),
                               AM = expression(log(AM(pi[0], t))))) +
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8, -10)) +
  labs(x = "Generation (t)",
       y = "log of Measure Value") +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.05),
        legend.justification = c("right", "bottom"),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(ylim = c(-10, 0))

results_nonres = expand_grid(
  t = ts,
  measures) |>
  mutate(boot = pmap(list(fn, t),
                     ~ with_boot(
                       nonres,
                       ..1,
                       R = 1000,
                       level_dad = dad_macro,
                       level_son = occ_macro,
                       t = ..2))) |>
  mutate(est = map_dbl(boot, "estimate"),
         se = map_dbl(boot, "se")) |>
  select(measure, t, est, se) |>
  mutate(dt_t = est*t)

write_csv(results_nonres, "data/advanced_measures_nonres.csv")

results_nonres = read_csv("data/advanced_measures_nonres.csv")

plot_nonres = ggplot(results_nonres, aes(x = t, y = est, color = measure, fill = measure)) +
  geom_ribbon(aes(ymin = est - 1.96*se,
                  ymax = est + 1.96*se),
              alpha = 0.2,
              linetype = 0) +
  geom_line(linetype = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("d" = "black",
                                "d_prime" = "orange",
                                "AM" = "blue"),
                     labels = c(d = expression(log(d(t))),
                                d_prime = expression(log(d * "'"~(t))),
                                AM = expression(log(AM(pi[0], t))))) +
  scale_fill_manual(values = c("d" = "black",
                               "d_prime" = "orange",
                               "AM" = "blue"),
                    labels = c(d = expression(log(d(t))),
                               d_prime = expression(log(d * "'"~(t))),
                               AM = expression(log(AM(pi[0], t))))) +
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8, -10)) +
  labs(x = "Generation (t)",
       y = "log of Measure Value") +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.05),
        legend.justification = c("right", "bottom"),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(ylim = c(-10, 0))

combined_measures = plot_res + plot_nonres

## IM CURVES

# RESERVATION

im_df_res = map_dfr(0:4, function(tt) {
  # call your existing function im()
  vals = im(res, dad_macro, occ_macro, t = tt)
  # vals is a named numeric vector: names are the father‐origins
  tibble(
    t      = tt,
    origin = names(vals),
    logIM  = vals
  )
}) |>
  mutate(origin = recode(origin,
                         farming   = "Farming",
                         blue_col  = "Manual",
                         white_col = "Nonmanual",
                         unemp     = "Not working"))

im_res = ggplot(im_df_res, aes(x = t, y = logIM, color = origin)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_color_brewer("Father origin", palette = "Dark2") +
  scale_x_continuous(breaks = 0:5) +
  labs(
    x     = "Generation (t)",
    y     = expression(log~IM(t,i)),
  ) +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.05),
        legend.justification = c("right", "bottom"),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(ylim = c(-8, 0))

# NONRESERVATION
im_df_nonres = map_dfr(0:4, function(tt) {
  # call your existing function im()
  vals = im(nonres, dad_macro, occ_macro, t = tt)
  # vals is a named numeric vector: names are the father‐origins
  tibble(
    t      = tt,
    origin = names(vals),
    logIM  = vals
  )
}) |>
  mutate(origin = recode(origin,
                         farming   = "Farming",
                         blue_col  = "Manual",
                         white_col = "Nonmanual",
                         unemp     = "Not working"))

im_nonres = ggplot(im_df_nonres, aes(x = t, y = logIM, color = origin)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_color_brewer("Father origin", palette = "Dark2") +
  scale_x_continuous(breaks = 0:5) +
  labs(
    x     = "Generation (t)",
    y     = expression(log~IM(t,i)),
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        legend.text = element_text(size = 13),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(ylim = c(-8, 0))

im_combined = im_res + im_nonres + plot_layout(widths = c(2, 2))

## SDM CURVES

variants = c("paired", "full", "resonly", "nonresonly")

results_sdm = expand_grid(
  t       = 0:4,
  variant = variants) |>
  mutate(
    boot = pmap(
      list(variant, t),
      ~ with_boot(
        data         = data,
        stat_fn      = sdm_mobility,
        R            = 1000,
        level_dad    = dad_macro,
        level_son    = occ_macro,
        subgroup_var = res_cty_ok,
        res_val      = 1,
        nonres_val   = 0,
        t            = ..2,
        variant      = ..1
      )
    ),
    estimate = map_dbl(boot, "estimate"),
    se       = map_dbl(boot, "se")
  ) |>
  select(-boot)

write_csv(results_sdm, "data/sdm_results.csv")

results_sdm = read_csv("data/sdm_results.csv") |>
  mutate(est = log(est))

sdm_df = expand_grid(
  t       = ts,
  variant = variants) |>
  mutate(
    logSDM = map2_dbl(
      variant, t,
      ~ sdm_mobility(
        data         = data,
        level_dad    = dad_macro,
        level_son    = occ_macro,
        subgroup_var = res_cty_ok,
        res_val      = 1,
        nonres_val   = 0,
        t            = .y,
        variant      = .x)))

sdm_curve = ggplot(sdm_df, aes(x = t, y = logSDM, color = variant, linetype = variant)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(
    name   = NULL,
    values = c(
      paired     = "black",
      full       = "firebrick",
      resonly    = "steelblue",
      nonresonly = "darkgreen"
    ),
    labels = c(
      paired     = expression(log~SDM(t, pi[0]^R, P[R],    pi[0]^N, P[N])),
      full       = expression(log~SDM(t, pi[0]^R, P,       pi[0]^N, P)),
      resonly    = expression(log~SDM(t, pi[0]^R, P[R],    pi[0]^N, P[R])),
      nonresonly = expression(log~SDM(t, pi[0]^R, P[N],    pi[0]^N, P[N]))
    )
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c(
      paired     = "dashed",
      full       = "solid",
      resonly    = "solid",
      nonresonly = "solid"
    ),
    labels = c(
      paired     = expression(log~SDM(t, pi[0]^R, P[R],    pi[0]^N, P[N])),
      full       = expression(log~SDM(t, pi[0]^R, P,       pi[0]^N, P)),
      resonly    = expression(log~SDM(t, pi[0]^R, P[R],    pi[0]^N, P[R])),
      nonresonly = expression(log~SDM(t, pi[0]^R, P[N],    pi[0]^N, P[N]))
    )
  ) +
  scale_x_continuous(breaks = ts) +
  labs(
    x     = "Generations (t)",
    y     = "log SDM(t)",
    title = "(Reservation vs. Non‐reservation Memory Curves)"
  ) +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.05),
        legend.justification = c("right", "bottom"),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(ylim = c(-8, 0))

s################################################################################

wtd_sd = function(x, w) {
  m = weighted.mean(x, w, na.rm = TRUE)
  sqrt( sum(w * (x - m)^2, na.rm = TRUE) / sum(w, na.rm = TRUE) )
}

demo = data |>
  summarise(age = weighted.mean(age, weight),
            res_cty_ok = weighted.mean(res_cty_ok, weight))

occ = data |>
  group_by(occ_macro) |>
  summarise(wsum = sum(weight), .groups = "drop") |>
  mutate(prop = wsum / sum(wsum)) |>
  arrange(occ_macro)

mobility_macro_by_region = data |>
  # 1) define upward vs. downward
  mutate(upward = case_when(
    dad_macro == "unemp" & occ_macro != "umemp" ~ 1,
    dad_macro == "farming" & occ_macro == "white_col" ~ 1,
    dad_macro == "blue_col"   & occ_macro == "white_col" ~ 1,
    TRUE ~ 0),
    downward = case_when(
      dad_macro == "white_col" & occ_macro != "white_col" ~ 1,
      dad_macro == "blue_col" & occ_macro == "unemp" ~ 1,
      dad_macro == "farming" & occ_macro == "umemp" ~ 1,
      TRUE ~ 0)) |>
  group_by(region) |>
  summarise(upward_rate   = weighted.mean(upward,   weight, na.rm = TRUE),
            downward_rate = weighted.mean(downward, weight, na.rm = TRUE),
            .groups = "drop")

mobility_meso_by_region = data |>
  # 1) define upward vs. downward
  mutate(
    upward = case_when(
      dad_meso == "farming"     & occ_meso %in% c("clerical", "prof")           ~ 1,
      dad_meso == "unskilled"   & occ_meso %in% c("crafts", "clerical", "prof") ~ 1,
      dad_meso == "crafts"      & occ_meso %in% c("clerical", "prof")           ~ 1,
      dad_meso == "clerical"    & occ_meso == "prof" ~ 1,
      dad_meso == "unemp" & occ_meso !=  "unemp" ~ 1,
      TRUE  ~ 0
    ),
    downward = case_when(
      # any non‐stay that isn’t classified as upward
      dad_meso != occ_meso & upward == 0 ~ 1,
      TRUE ~ 0
    )
  ) |>
  # 2) aggregate by region, using your sampling‐weight `weight`
  group_by(region) |>
  summarise(
    upward_rate   = weighted.mean(upward,   weight, na.rm = TRUE),
    downward_rate = weighted.mean(downward, weight, na.rm = TRUE),
    .groups = "drop"
  )

mobility_by_region = data |>
  mutate(
    mobility = case_when(dad_macro != occ_macro ~ 1,
                       TRUE ~ 0)) |>
  group_by(region) |>
  summarise(
    mobility_rate = weighted.mean(mobility, weight, na.rm = TRUE),
    .groups = "drop")

states_sf = states(cb = TRUE) |>
  filter(!STUSPS %in% c("AK","HI","PR","GU","VI","MP", "AS")) |>
  st_transform(crs = 5070) |>
  mutate(statefip_1940 = as.integer(STATEFP),
         region = case_when(
           statefip_1940 %in% c(8, 16, 32, 49, 56) ~ "basin",
           statefip_1940 == 6 ~ "cali",
           statefip_1940 %in% c(27, 55) ~ "lakes",
           statefip_1940 %in% c(17, 18, 26, 39) ~ "midwest",
           statefip_1940 == 37 ~ "nc",
           statefip_1940 %in% c(9,10,23,24,25,33,34,36,42,44,50,11) ~ "ne",
           statefip_1940 %in% c(30, 38, 46) ~ "plains",
           statefip_1940 %in% c(41, 53) ~ "nw",
           statefip_1940 == 40 ~ "ok",
           statefip_1940 %in% c(19, 20, 29, 31) ~ "prairie",
           statefip_1940 %in% c(1,5,12,13,21,22,28,45,47,48,51,54) ~ "south",
           statefip_1940 %in% c(4, 35) ~ "sw",
           TRUE ~ NA_character_)) |>
  left_join(mobility_macro_by_region, by = "region") |>
  left_join(mobility_by_region, by = "region")

ggplot(states_sf) +
  geom_sf(aes(fill = upward_rate), color = "grey40") +
  scale_fill_gradient(name = "Upward Mobility",
                      low = "yellow", high = "red",
                      labels = scales::percent_format(accuracy = 1)) +
  theme_void()

ggplot(states_sf) +
  geom_sf(aes(fill = downward_rate), color = "grey40") +
  scale_fill_gradient(name = "Downward Mobility",
                      low = "lavender", high = "purple",
                      labels = scales::percent_format(accuracy = 1)) +
  theme_void()

ggplot(states_sf) +
  geom_sf(aes(fill = mobility_rate), color = "grey40") +
  scale_fill_gradient(name = "Total Mobility",
                      low = "lightblue", high = "darkblue",
                      labels = scales::percent_format(accuracy = 1)) +
  theme_void()



