library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(expm)
library(rlang)
library(purrr)
library(tigris)


data = read_csv("data/aian_weighted.csv")

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
  
  return(mat)
}

pi_star = function(p_mat) {
  P = as.matrix(p_mat)
  eig = eigen(t(P))
  idx = which.min(abs(eig$values - 1))
  v = Re(eig$vectors[, idx])
  if (any(v < 0)) v = abs(v)
  pi_star = v / sum(v)
  names(pi_star) = rownames(P)
  return(pi_star)
}

p_upward = function(data,
                    level = c("macro", "meso"),
                    farming = TRUE,
                    unemp = FALSE,
                    weight = weight) {
  
  level = match.arg(level)
  
  data |>
    mutate(upward =
             (unemp & dad_macro == "unemp" & occ_macro != "unemp") |
             (farming & dad_macro == "farming" & occ_macro == "white_col") |
             if (level == "macro")
               (dad_macro == "blue_col" & occ_macro == "white_col")
           else 
             (dad_meso == "unskilled" & (occ_macro == "white_col" | occ_meso == "crafts")) |
             (dad_meso == "clerical" & occ_meso == "prof")) |>
    summarise(upward = weighted.mean(as.numeric(upward), w = {{ weight }})) |>
    pull(upward)
}

upward = function(p_mat,
                  level = c("macro", "meso"),
                  farming = TRUE,
                  unemp = FALSE)

p_downward = function(data,
                      level = c("macro", "meso"),
                      farming = TRUE,
                      unemp = FALSE,
                      weight = weight) {
  
  level = match.arg(level)
  
  data |>
    mutate(downward = 
             (unemp & dad_macro != "unemp" & occ_macro == "unemp") |
             (farming & dad_macro == "white_col" & occ_macro == "farming") |
             if (level == "macro")
               (dad_macro == "white_col" & occ_macro == "blue_col")
           else
             ((dad_macro == "white_col" | dad_meso == "crafts") & occ_meso == "unskilled") |
             (dad_meso == "prof" & occ_meso == "clerical")) |>
    summarise(downward = weighted.mean(as.numeric(downward), w = {{ weight }})) |>
    pull(downward)
}

shorrocks = function(data, level_dad, level_son, weight = weight) {
  data |>
    mutate(.same = {{ level_dad }} != {{ level_son }}) |>
    summarise(prop = weighted.mean(as.numeric(.same), w = {{ weight }}, na.rm = TRUE)) |>
    pull(prop)
}

first_passage = function(p_mat) {
  P = as.matrix(p_mat)
}

################################################################################
############################## ADVANCED MEASURES ###############################
################################################################################

with_boot = function(data, stat_fn, R = 200, ...) {
  args_q = rlang::enquos(...)
  N = nrow(data)
  stat_call = function(d) rlang::eval_tidy(rlang::call2(stat_fn, data = d, !!!args_q))
  
  est = stat_call(data)
  
  boots = replicate(R, {
    d_b = data[sample.int(N, N, replace = TRUE), , drop = FALSE]
    stat_call(d_b)
  }, simplify = TRUE)  # for these stats, boots is a numeric vector (scalar per draw)
  
  se = sd(boots, na.rm = TRUE)
  
  list(estimate = est, se = se, draws = boots)
}


boot_pmatrix_ci = function(
    data, level_dad, level_son,
    R = 1000, conf = 0.95, .seed = NULL
){
  if (!is.null(.seed)) set.seed(.seed)
  
  dad_sym = rlang::ensym(level_dad)
  son_sym = rlang::ensym(level_son)
  
  # point estimate
  P_hat = p_matrix(data, !!dad_sym, !!son_sym, matrix = TRUE)
  rnames = rownames(P_hat); cnames = colnames(P_hat)
  nR = nrow(P_hat); nC = ncol(P_hat)
  
  # one classical bootstrap draw: resample indices
  boot_once = function() {
    idx = sample.int(nrow(data), nrow(data), replace = TRUE)
    d_b = data[idx, , drop = FALSE]
    p_matrix(d_b, !!dad_sym, !!son_sym, matrix = TRUE)
  }
  
  boots = replicate(R, boot_once(), simplify = FALSE)
  arr   = simplify2array(boots)  # [rows x cols x R]
  
  alpha = 1 - conf
  se_mat = apply(arr, c(1, 2), sd, na.rm = TRUE)
  q_lo   = apply(arr, c(1, 2), quantile, probs = alpha/2, na.rm = TRUE, names = FALSE)
  q_hi   = apply(arr, c(1, 2), quantile, probs = 1 - alpha/2, na.rm = TRUE, names = FALSE)
  
  tibble::tibble(
    !!dad_sym := rep(rnames, times = nC),
    !!son_sym := rep(cnames, each  = nR),
    est = as.vector(P_hat),
    se  = as.vector(se_mat),
    lo  = as.vector(q_lo),
    hi  = as.vector(q_hi)
  )
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

## ---- helpers to find generators of d(t) and d'(t) ----------------------------

d_generator = function(P_t, pi_star) {
  # returns the class(es) i with max || P_t(i,·) - π* ||_TV and the value
  scores = apply(P_t, 1, function(r) tv_norm(r, pi_star))
  mx = max(scores)
  i_star = which(abs(scores - mx) < 1e-12)
  list(classes = rownames(P_t)[i_star], value = mx)
}

dprime_generator = function(P_t) {
  # returns the pair(s) (i,j) with max || P_t(i,·) - P_t(j,·) ||_TV and the value
  n = nrow(P_t)
  best = -Inf
  keep = list()
  for (i in 1:(n-1)) for (j in (i+1):n) {
    v = tv_norm(P_t[i, ], P_t[j, ])
    if (v > best + 1e-12) {
      best = v
      keep = list(c(i, j))
    } else if (abs(v - best) <= 1e-12) {
      keep = append(keep, list(c(i, j)))
    }
  }
  pairs_named = lapply(keep, \(idx) rownames(P_t)[idx])
  list(pairs = pairs_named, value = best)
}

identify_generators = function(data, level_dad, level_son, ts = 0:5) {
  P = p_matrix(data, {{ level_dad }}, {{ level_son }})
  piS = pi_star(P)
  
  purrr::map_dfr(ts, function(tt) {
    P_t = if (tt == 0) diag(nrow(P)) else P %^% tt
    rownames(P_t) = rownames(P)
    
    dgen  = d_generator(P_t, piS)
    dpgen = dprime_generator(P_t)
    
    tibble::tibble(
      t               = tt,
      d_value         = dgen$value,
      d_classes       = paste(dgen$classes, collapse = " | "),
      dprime_value    = dpgen$value,
      dprime_pairs    = paste(
        vapply(dpgen$pairs, function(p) paste(p, collapse = " vs "), character(1)),
        collapse = "  |  "
      )
    )
  })
}

################################################################################
############################### IMPLEMENTATION #################################
################################################################################

## MACRO
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

## MESO
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

## MID

## RESERVATION ONLY

ci_res = boot_pmatrix_ci(
  res, dad_mid, occ_mid,
  R = 1000, conf = 0.95, .seed = 123
) |>
  mutate(
    dad_mid = factor(dad_mid, levels = rev(mid_order)),
    occ_mid = factor(occ_mid, levels = mid_order)
  )

g_mid_res = ggplot(ci_res, aes(x = occ_mid, y = dad_mid, fill = est)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n[%.2f, %.2f]", est, lo, hi)),
            vjust = 0.3, size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(
    x = "Son occupation",
    y = "Father occupation",
    fill = "Transition Prob.",
    title = "P"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

pi0_vec_mid_res = pi_0(res, dad_mid)
df_pi0_mid_res = tibble(
  father = names(pi0_vec_mid_res),
  pi0    = as.numeric(pi0_vec_mid_res)
) |>
  mutate(father = factor(father, levels = rev(mid_order)))

g0_mid_res = ggplot(df_pi0_mid_res, aes(x = 1, y = father, fill = pi0)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi0)), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(title = expression(pi[0])) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

steady_mid_res = pi_star(p_matrix(res, dad_mid, occ_mid, TRUE))
df_pi_star_mid_res = tibble(
  father  = names(steady_mid_res),
  pi_star = as.numeric(steady_mid_res)
) |>
  mutate(father = factor(father, levels = rev(mid_order)))

g_star_mid_res = ggplot(df_pi_star_mid_res, aes(x = 1, y = father, fill = pi_star)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi_star)), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(title = expression(pi^"*")) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

combined_plot_mid_res = g_mid_res + g0_mid_res + g_star_mid_res +
  plot_layout(widths = c(6, 1, 1))


## NONRESERVATION ONLY

ci_nonres = boot_pmatrix_ci(
  nonres, dad_mid, occ_mid,
  R = 1000, conf = 0.95, .seed = 123
) |>
  mutate(
    dad_mid = factor(dad_mid, levels = rev(mid_order)),
    occ_mid = factor(occ_mid, levels = mid_order)
  )

g_mid_nonres = ggplot(ci_nonres, aes(x = occ_mid, y = dad_mid, fill = est)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n[%.2f, %.2f]", est, lo, hi)),
            vjust = 0.3, size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(
    x = "Son occupation",
    y = "Father occupation",
    fill = "Transition Prob.",
    title = "P"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

pi0_vec_mid_nonres = pi_0(nonres, dad_mid)
df_pi0_mid_nonres = tibble(
  father = names(pi0_vec_mid_nonres),
  pi0    = as.numeric(pi0_vec_mid_nonres)
) |>
  mutate(father = factor(father, levels = rev(mid_order)))

g0_mid_nonres = ggplot(df_pi0_mid_nonres, aes(x = 1, y = father, fill = pi0)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi0)), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(title = expression(pi[0])) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

steady_mid_nonres = pi_star(p_matrix(nonres, dad_mid, occ_mid, TRUE))
df_pi_star_mid_nonres = tibble(
  father  = names(steady_mid_nonres),
  pi_star = as.numeric(steady_mid_nonres)
) |>
  mutate(father = factor(father, levels = rev(mid_order)))

g_star_mid_nonres = ggplot(df_pi_star_mid_nonres, aes(x = 1, y = father, fill = pi_star)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi_star)), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(title = expression(pi^"*")) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

combined_plot_mid_nonres = g_mid_nonres + g0_mid_nonres + g_star_mid_nonres +
  plot_layout(widths = c(6, 1, 1))


################################################################################

data = data |> mutate(res_cty_ok = ifelse(res_cty == 0 | statefip_1940 == 40, 0, 1))

## d(t), d'(t), AIM CURVES

ts = 0:5
measures = tibble(
  measure = c("d", "d_prime", "AM"),
  fn = list(d_t, d_prime, am)
)

results_res = expand_grid(t = ts, measures) |>
  mutate(boot = pmap(list(fn, t),
                     ~ with_boot(res, ..1, R = 1000,
                                 level_dad = dad_mid,
                                 level_son = occ_mid,
                                 t = ..2))) |>
  mutate(
    est   = map_dbl(boot, "estimate"),
    se    = map_dbl(boot, "se"),
    lo    = map_dbl(boot, ~ quantile(.x$draws, 0.025, na.rm = TRUE)),
    hi    = map_dbl(boot, ~ quantile(.x$draws, 0.975, na.rm = TRUE))
  ) |>
  select(measure, t, est, se, lo, hi) |>
  mutate(dt_t = est * t)

write_csv(results_res, "data/advanced_measures_res.csv")

results_res = read_csv("data/advanced_measures_res.csv")

plot_res = ggplot(results_res, aes(x = t, y = est, color = measure, fill = measure)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, linetype = 0) +
  geom_line(linetype = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c(d = "#332288", d_prime = "#009E73", AM = "darkred"),
                     labels = c(d = expression(log(d(t))),
                                d_prime = expression(log(d*"'"*(t))),
                                AM = expression(log(AM(pi[0], t))))) +
  scale_fill_manual(values = c(d = "#332288", d_prime = "#009E73", AM = "darkred"),
                    labels = c(d = expression(log(d(t))),
                               d_prime = expression(log(d*"'"*(t))),
                               AM = expression(log(AM(pi[0], t))))) +
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8, -10)) +
  labs(x = "Generation (t)", y = "log of Measure Value") +
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

results_nonres = expand_grid(t = ts, measures) |>
  mutate(boot = pmap(list(fn, t),
                     ~ with_boot(nonres, ..1, R = 1000,
                                 level_dad = dad_mid,
                                 level_son = occ_mid,
                                 t = ..2))) |>
  mutate(
    est   = map_dbl(boot, "estimate"),
    se    = map_dbl(boot, "se"),
    lo    = map_dbl(boot, ~ quantile(.x$draws, 0.025, na.rm = TRUE)),
    hi    = map_dbl(boot, ~ quantile(.x$draws, 0.975, na.rm = TRUE))
  ) |>
  select(measure, t, est, se, lo, hi) |>
  mutate(dt_t = est * t)

write_csv(results_nonres, "data/advanced_measures_nonres.csv")

results_nonres = read_csv("data/advanced_measures_nonres.csv")

plot_nonres = ggplot(results_nonres, aes(x = t, y = est, color = measure, fill = measure)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, linetype = 0) +
  geom_line(linetype = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c(d = "#332288", d_prime = "#009E73", AM = "darkred"),
                     labels = c(d = expression(log(d(t))),
                                d_prime = expression(log(d*"'"*(t))),
                                AM = expression(log(AM(pi[0], t))))) +
  scale_fill_manual(values = c(d = "#332288", d_prime = "#009E73", AM = "darkred"),
                    labels = c(d = expression(log(d(t))),
                               d_prime = expression(log(d*"'"*(t))),
                               AM = expression(log(AM(pi[0], t))))) +
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8, -10)) +
  labs(x = "Generation (t)", y = "log of Measure Value") +
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

boot_im_by_t <- function(data, level_dad, level_son, ts = 0:5, R = 1000, .seed = NULL) {
  if (!is.null(.seed)) set.seed(.seed)
  
  dad_sym <- rlang::ensym(level_dad)
  son_sym <- rlang::ensym(level_son)
  
  # point estimates first
  P0   <- p_matrix(data, !!dad_sym, !!son_sym, matrix = TRUE)
  piS0 <- pi_star(P0)
  rowlabs <- rownames(P0)
  
  point_by_t <- purrr::map(ts, function(tt) {
    Pt <- if (tt == 0) diag(nrow(P0)) else P0 %^% tt
    im_vals <- apply(Pt, 1, function(r) tv_norm(r, piS0))
    tibble::tibble(t = tt, origin = rowlabs, est = log(im_vals))
  }) |> dplyr::bind_rows()
  
  # bootstrap draws (classical): resample indices, recompute P, pi*, and IM(t,i)
  N <- nrow(data)
  boot_once <- function() {
    idx <- sample.int(N, N, replace = TRUE)
    db  <- data[idx, , drop = FALSE]
    P   <- p_matrix(db, !!dad_sym, !!son_sym, matrix = TRUE)
    piS <- pi_star(P)
    
    # for efficiency, walk powers iteratively
    Pt <- diag(nrow(P))
    out_list <- vector("list", length(ts))
    for (k in seq_along(ts)) {
      tt <- ts[k]
      if (tt > 0) {
        # incrementally multiply from previous Pt (safer than P%^%tt in a loop)
        if (k == 1) {  # if first t > 0
          Pt <- P %^% tt
        } else {
          # if ts is strictly increasing by 1, you can do Pt <- Pt %*% P
          # but to be robust to arbitrary ts, recompute only when needed:
          Pt <- P %^% tt
        }
      } else {
        Pt <- diag(nrow(P))
      }
      im_vals <- apply(Pt, 1, function(r) tv_norm(r, piS))
      out_list[[k]] <- log(im_vals)
    }
    # returns a matrix [origins x length(ts)]
    do.call(cbind, out_list)
  }
  
  boots <- replicate(R, boot_once(), simplify = FALSE)  # list of matrices
  # stack to array: [origin x t x R]
  arr <- simplify2array(boots)
  # names
  dimnames(arr) <- list(origin = rowlabs, t = as.character(ts), rep = NULL)
  
  # summarize: SE + percentile CI per (origin, t)
  alpha <- 0.05
  summ <- lapply(ts, function(tt) {
    a2 <- arr[, as.character(tt), , drop = FALSE]  # [origin x 1 x R]
    draws_mat <- drop(a2)                           # [origin x R]
    se <- apply(draws_mat, 1, sd, na.rm = TRUE)
    lo <- apply(draws_mat, 1, quantile, probs = alpha/2,     na.rm = TRUE, names = FALSE)
    hi <- apply(draws_mat, 1, quantile, probs = 1 - alpha/2, na.rm = TRUE, names = FALSE)
    tibble::tibble(t = tt, origin = rowlabs, se = unname(se), lo = unname(lo), hi = unname(hi))
  }) |> dplyr::bind_rows()
  
  dplyr::left_join(point_by_t, summ, by = c("t", "origin"))
}

tol_bright <- c(
  "#4477AA", # blue
  "#EE6677", # red
  "#228833", # green
  "#CCBB44", # yellow
  "#66CCEE"  # cyan
)

# RESERVATION
im_boot_res_mid <- boot_im_by_t(res, dad_mid, occ_mid, ts = 0:5, R = 1000, .seed = 123) |>
  dplyr::mutate(origin = dplyr::recode(origin,
                                       farming="Farming", unskilled="Unskilled", crafts="Crafts",
                                       nonmanual="Nonmanual", unemp="Not working"
  ))

im_res_mid_plot <- ggplot(im_boot_res_mid, aes(x = t, y = est, color = origin)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = origin), alpha = 0.15, linetype = 0, color = NA) +
  geom_line(linetype = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = tol_bright) +
  scale_fill_manual(values = tol_bright) +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8, -10)) +
  labs(x = "Generation (t)", y = expression(log~IM(t,i))) +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.05),
        legend.justification = c("right","bottom"),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(ylim = c(-10, 0))

# NONRESERVATION
im_boot_nonres_mid <- boot_im_by_t(nonres, dad_mid, occ_mid, ts = 0:5, R = 1000, .seed = 123) |>
  dplyr::mutate(origin = dplyr::recode(origin,
                                       farming="Farming", unskilled="Unskilled", crafts="Crafts",
                                       nonmanual="Nonmanual", unemp="Not working"
  ))

im_nonres_mid_plot <- ggplot(im_boot_nonres_mid, aes(x = t, y = est, color = origin)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = origin), alpha = 0.15, linetype = 0, color = NA) +
  geom_line(linetype = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = tol_bright) +
  scale_fill_manual(values = tol_bright) +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8, -10)) +
  labs(x = "Generation (t)", y = expression(log~IM(t,i))) +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.05),
        legend.justification = c("right","bottom"),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5)) +
  coord_cartesian(ylim = c(-10, 0))

im_combined <- im_res_mid_plot + im_nonres_mid_plot + 
  plot_layout(widths = c(2, 2))

################################################################################

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

