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
  if (any(v < 0)) v = abs(v)
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

with_boot_bb = function(data, stat_fn, R = 1000, ..., prior = c("plain", "weight"), .parallel = FALSE, .seed = NULL) {
  if (!is.null(.seed)) set.seed(.seed)
  prior = match.arg(prior)
  N = nrow(data)
  args_q = rlang::enquos(...)
  
  stat_call = function(d) rlang::eval_tidy(rlang::call2(stat_fn, data = d, !!!args_q))
  
  # One BB draw: multiply existing weights by independent Gammas
  bb_once = function() {
    if (prior == "plain") {
      g = rexp(N, rate = 1)
    } else {
      # prior mass proportional to sampling weights (shift tiny eps to avoid zeros)
      a = pmax(data$weight, 1e-12)
      g = rgamma(N, shape = a, rate = 1)
    }
    d_b = data
    d_b$weight = data$weight * g
    stat_call(d_b)
  }
  
  boots = if (.parallel) {
    parallel::mclapply(seq_len(R), function(i) bb_once(), mc.cores = max(1L, parallel::detectCores() - 1L))
  } else {
    replicate(R, bb_once(), simplify = FALSE)
  }
  
  est = stat_call(data)
  boots_arr = simplify2array(boots)
  
  se =
    if (is.atomic(est) && length(est) == 1L) {
      sd(boots_arr, na.rm = TRUE)
    } else if (is.matrix(boots_arr) || (is.array(boots_arr) && length(dim(boots_arr)) == 2L)) {
      apply(boots_arr, 1L, sd, na.rm = TRUE)
    } else {
      stop("with_boot_bb(): unsupported output type from stat_fn()")
    }
  
  list(estimate = est, se = se)
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

results = read_csv("data/advanced_measures.csv")

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
  scale_color_manual(values = c("d" = "#332288",
                                "d_prime" = "#009E73",
                                "AM" = "darkred"),
                     labels = c(d = expression(log(d(t))),
                                d_prime = expression(log(d*"'"*(t))),
                                AM = expression(log(AM(pi[0], t))))) +
  scale_fill_manual(values = c("d" = "#332288",
                               "d_prime" = "#009E73",
                               "AM" = "darkred"),
                    labels = c(d = expression(log(d(t))),
                               d_prime = expression(log(d*"'"*(t))),
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
  scale_color_manual(values = c("d" = "#332288",
                                "d_prime" = "#009E73",
                                "AM" = "darkred"),
                     labels = c(d = expression(log(d(t))),
                                d_prime = expression(log(d*"'"*(t))),
                                AM = expression(log(AM(pi[0], t))))) +
  scale_fill_manual(values = c("d" = "#332288",
                               "d_prime" = "#009E73",
                               "AM" = "darkred"),
                    labels = c(d = expression(log(d(t))),
                               d_prime = expression(log(d*"'"*(t))),
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

