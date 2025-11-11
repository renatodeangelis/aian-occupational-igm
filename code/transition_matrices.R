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
library(jtools)

data = read_csv("data/aian_weighted.csv")

## DATA DESCRIPTION

desc_main = data |>
  mutate(age_pop = 1940 - birthyr_pop,
         age_son = 1940 - birthyr_son) |>
  summarise(age_pop_mean = weighted.mean(age_pop, w_atc_norm),
            age_pop_sd = wtd.sd(age_pop, w_atc_norm),
            age_son_mean = weighted.mean(age_son, w_atc_norm),
            age_son_sd = wtd.sd(age_son, w_atc_norm))

desc_region = data |>
  group_by(region) |>
  summarise(wsum = sum(w_atc_norm), .groups = "drop") |>
  mutate(prop = round(wsum / sum(wsum) * 100, 1)) |>
  select(-wsum) |>
  arrange(region)

desc_meso_pop = data |>
  group_by(meso_pop) |>
  summarise(wsum = sum(w_atc_norm), .groups = "drop") |>
  mutate(prop = wsum / sum(wsum) * 100) |>
  select(-wsum) |>
  arrange(meso_pop)

desc_meso_son = data |>
  group_by(meso_son) |>
  summarise(wsum = sum(w_atc_norm), .groups = "drop") |>
  mutate(prop = wsum / sum(wsum) * 100) |>
  select(-wsum) |>
  arrange(meso_son)

desc_educd = data |>
  group_by(education) |>
  summarise(wsum = sum(w_atc_norm), .groups = "drop") |>
  mutate(prop = wsum / sum(wsum) * 100) |>
  select(-wsum) |>
  arrange(education)

################################################################################
############################### BASIC MEASURES #################################
################################################################################

pi_0 = function(data, level) {
  df = data |> 
    group_by({{ level }}) |>
    summarise(total_w = sum(w_atc_norm),
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
      n_w = sum(w_atc_norm),
      .groups = "drop") |>
    group_by( !!dad_sym ) |>
    mutate(
      n_i_w = sum(n_w),
      P     = ifelse(n_i_w > 0, n_w / n_i_w, 0)) |>
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

################################################################################
############################## ADVANCED MEASURES ###############################
################################################################################

with_boot = function(data, stat_fn, R = 1000, .seed = NULL, ...) {
  if (!is.null(.seed)) set.seed(.seed)
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
    R = 1000, conf = 0.95, .seed = NULL){
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

mu_t = function(pi0, P, t = 0) {
  P = as.matrix(P)
  if (t == 0) return(as.numeric(pi0))
  as.numeric(pi0 %*% (P %^% t))
}

om = function(P, pi0, t) {
  P = as.matrix(P)
  mu_t = mu_t(pi0, P, t)
  return(1 - sum(mu_t * diag(P)))
}

sm = function(P, pi0, t) {
  P = as.matrix(P)
  mu_t = mu_t(pi0, P, t)
  mu_t1 = as.numeric(mu_t %*% P)
  return(tv_norm(mu_t, mu_t1))
}

em = function(P, pi0, t) {
  OM = om(P, pi0, t)
  SM = sm(P, pi0, t)
  
  return(OM - SM)
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

identify_generators = function(data, level_dad, level_son, ts = 0:4) {
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
p_mat_macro = boot_pmatrix_ci(
  data, macro_pop, macro_son,
  R = 1000, conf = 0.95, .seed = 123)

g_macro = ggplot(
  p_mat_macro |> 
    mutate(macro_pop = factor(macro_pop, levels = rev(unique(macro_pop)))),
  aes(x = macro_son, y = macro_pop, fill = est)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n[%.2f, %.2f]", est, lo, hi)),
            vjust = 0.3, size = 4) +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(x = "Son's occupation",
       y = "Father's occupation",
       fill = "Transition Prob.",
       title = expression(P)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 12),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

pi0_vec_macro = pi_0(data, macro_pop)
df_pi0_macro = tibble(
  father = names(pi0_vec_macro),
  pi0 = as.numeric(pi0_vec_macro)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g0_macro = ggplot(df_pi0_macro, aes(x = 1, y = father, fill = pi0)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi0)),
            size = 4, color = "black") +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(title = expression(pi[0])) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

steady_macro = pi_star(p_matrix(data, macro_pop, macro_son, TRUE))
df_pi_star_macro = tibble(
  father = names(steady_macro),
  pi_star = as.numeric(steady_macro)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g_star_macro = ggplot(df_pi_star_macro, aes(x = 1, y = father, fill = pi_star)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi_star)),
            size = 4, color = "black") +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(title = expression(pi^"*")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

combined_plot_macro = g_macro + g0_macro + g_star_macro + plot_layout(widths = c(6, 1, 1))

## MESO
p_mat_meso = boot_pmatrix_ci(
  data, meso_pop, meso_son,
  R = 1000, conf = 0.95, .seed = 123)

g_meso = ggplot(
  p_mat_meso |> 
    mutate(meso_pop = factor(meso_pop, levels = c("unem", "prof",
                                                  "clerical", "crafts",
                                                  "unskilled", "farmer",
                                                  "farmworker")),
           meso_son = factor(meso_son, levels = rev(c("unemp", "prof",
                                                      "clerical", "crafts",
                                                      "unskilled", "farmer",
                                                  "farmworker")))),
  aes(x = meso_son, y = meso_pop, fill = est)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n[%.2f, %.2f]", est, lo, hi)),
            vjust = 0.3, size = 3) +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(x = "Son's occupation",
       y = "Father's occupation",
       fill = "Transition Prob.",
       title = "P") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

pi0_vec_meso = pi_0(data, meso_pop)
df_pi0_meso = tibble(
  father = names(pi0_vec_meso),
  pi0 = as.numeric(pi0_vec_meso)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g0_meso = ggplot(df_pi0_meso |>
                   mutate(father = factor(father, levels = c("unemp", "prof",
                                                                 "clerical", "crafts",
                                                                 "unskilled", "farmer",
                                                                 "farmworker"))),
                          aes(x = 1, y = father, fill = pi0)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi0)),
            size = 4, color = "black") +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(title = expression(pi[0])) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

steady_meso = pi_star(p_matrix(data, meso_pop, meso_son, TRUE))
df_pi_star_meso = tibble(
  father = names(steady_meso),
  pi_star = as.numeric(steady_meso)) |>
  mutate(father = factor(father, levels = rev(unique(father))))
g_star_meso = ggplot(df_pi_star_meso |>
                       mutate(father = factor(father, levels = c("unemp", "prof",
                                                                     "clerical", "crafts",
                                                                     "unskilled", "farmer",
                                                                     "farmworker"))),
                     aes(x = 1, y = father, fill = pi_star)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi_star)),
            size = 4, color = "black") +
  scale_fill_gradient(low = "lightyellow", high = "firebrick") +
  labs(title = expression(pi^"*")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

combined_plot_meso = g_meso + g0_meso + g_star_meso + plot_layout(widths = c(6, 1, 1))

################################################################################

## d(t), d'(t), AM CURVES

ts = 0:4
measures = tibble(
  measure = c("d", "d_prime", "AM"),
  fn = list(d_t, d_prime, am))

results_macro = expand_grid(t = ts, measures) |>
  mutate(boot = pmap(list(fn, t),
                     ~ with_boot(data, ..1, R = 1000,
                                 .seed = 123,
                                 level_dad = macro_pop,
                                 level_son = macro_son,
                                 t = ..2))) |>
  mutate(
    est   = map_dbl(boot, "estimate"),
    se    = map_dbl(boot, "se"),
    lo    = map_dbl(boot, ~ quantile(.x$draws, 0.025, na.rm = TRUE)),
    hi    = map_dbl(boot, ~ quantile(.x$draws, 0.975, na.rm = TRUE))) |>
  select(measure, t, est, se, lo, hi) |>
  mutate(dt_t = est * t)

measure_plot_macro = ggplot(results_macro, aes(x = t, y = est, color = measure, fill = measure)) +
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
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8)) +
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
  coord_cartesian(ylim = c(-8, 0))

results_meso = expand_grid(t = ts, measures) |>
  mutate(boot = pmap(list(fn, t),
                     ~ with_boot(data, ..1, R = 1000,
                                 .seed = 123,
                                 level_dad = meso_pop,
                                 level_son = meso_son,
                                 t = ..2))) |>
  mutate(
    est   = map_dbl(boot, "estimate"),
    se    = map_dbl(boot, "se"),
    lo    = map_dbl(boot, ~ quantile(.x$draws, 0.025, na.rm = TRUE)),
    hi    = map_dbl(boot, ~ quantile(.x$draws, 0.975, na.rm = TRUE))) |>
  select(measure, t, est, se, lo, hi) |>
  mutate(dt_t = est * t)

measure_plot_meso = ggplot(results_meso, aes(x = t, y = est, color = measure, fill = measure)) +
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
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8)) +
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
  coord_cartesian(ylim = c(-8, 0))

combined_measures = measure_plot_macro + measure_plot_meso

## IM CURVES

boot_im_by_t = function(data, level_dad, level_son, ts = 0:4, R = 1000, .seed = NULL) {
  if (!is.null(.seed)) set.seed(.seed)
  
  dad_sym = rlang::ensym(level_dad)
  son_sym = rlang::ensym(level_son)
  
  # point estimates first
  P0   = p_matrix(data, !!dad_sym, !!son_sym, matrix = TRUE)
  piS0 = pi_star(P0)
  rowlabs = rownames(P0)
  
  point_by_t = purrr::map(ts, function(tt) {
    Pt = if (tt == 0) diag(nrow(P0)) else P0 %^% tt
    im_vals = apply(Pt, 1, function(r) tv_norm(r, piS0))
    tibble::tibble(t = tt, origin = rowlabs, est = log(im_vals))
  }) |> dplyr::bind_rows()
  
  # bootstrap draws (classical): resample indices, recompute P, pi*, and IM(t,i)
  N = nrow(data)
  boot_once = function() {
    idx = sample.int(N, N, replace = TRUE)
    db  = data[idx, , drop = FALSE]
    P   = p_matrix(db, !!dad_sym, !!son_sym, matrix = TRUE)
    piS = pi_star(P)
    
    # for efficiency, walk powers iteratively
    Pt = diag(nrow(P))
    out_list = vector("list", length(ts))
    for (k in seq_along(ts)) {
      tt = ts[k]
      if (tt > 0) {
        # incrementally multiply from previous Pt (safer than P%^%tt in a loop)
        if (k == 1) {  # if first t > 0
          Pt = P %^% tt
        } else {
          # if ts is strictly increasing by 1, you can do Pt = Pt %*% P
          # but to be robust to arbitrary ts, recompute only when needed:
          Pt = P %^% tt
        }
      } else {
        Pt = diag(nrow(P))
      }
      im_vals = apply(Pt, 1, function(r) tv_norm(r, piS))
      out_list[[k]] = log(im_vals)
    }
    # returns a matrix [origins x length(ts)]
    do.call(cbind, out_list)
  }
  
  boots = replicate(R, boot_once(), simplify = FALSE)  # list of matrices
  # stack to array: [origin x t x R]
  arr = simplify2array(boots)
  # names
  dimnames(arr) = list(origin = rowlabs, t = as.character(ts), rep = NULL)
  
  # summarize: SE + percentile CI per (origin, t)
  alpha = 0.05
  summ = lapply(ts, function(tt) {
    a2 = arr[, as.character(tt), , drop = FALSE]  # [origin x 1 x R]
    draws_mat = drop(a2)                           # [origin x R]
    se = apply(draws_mat, 1, sd, na.rm = TRUE)
    lo = apply(draws_mat, 1, quantile, probs = alpha/2,     na.rm = TRUE, names = FALSE)
    hi = apply(draws_mat, 1, quantile, probs = 1 - alpha/2, na.rm = TRUE, names = FALSE)
    tibble::tibble(t = tt, origin = rowlabs, se = unname(se), lo = unname(lo), hi = unname(hi))
  }) |> dplyr::bind_rows()
  
  dplyr::left_join(point_by_t, summ, by = c("t", "origin"))
}

# MACRO
im_boot_macro = boot_im_by_t(data, macro_pop, macro_son, ts = 0:4, R = 1000, .seed = 123) |>
  dplyr::mutate(origin = dplyr::recode(origin,
                                       farming = "Farming", manual = "Manual",
                                       nonmanual = "Nonmanual", unemp = "Not working"))

im_macro_plot = ggplot(im_boot_macro, aes(x = t, y = est, color = origin)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = origin), alpha = 0.15, linetype = 0) +
  geom_line(linetype = 1) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8)) +
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
  coord_cartesian(ylim = c(-8, 0))

# MESO
im_boot_meso = boot_im_by_t(data, meso_pop, meso_son, ts = 0:4, R = 1000, .seed = 123) |>
  dplyr::mutate(origin = dplyr::recode(origin,
                                       farmer = "Farmer", unskilled = "Unskilled",
                                       crafts = "Crafts", prof = "Professional",
                                       clerical = "Clerical", unemp = "Not working",
                                       farmworker = "Farmworker"))

im_meso_plot = ggplot(im_boot_meso, aes(x = t, y = est, color = origin)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = origin), alpha = 0.15, linetype = 0) +
  geom_line(linetype = 1) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(breaks = c(0, -2, -4, -6, -8)) +
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
  coord_cartesian(ylim = c(-8, 0))

im_combined = im_macro_plot + im_meso_plot + 
  plot_layout(widths = c(2, 2))

## COUNTY ESTIMATION

p_upward = function(data,
                    level = c("macro", "meso"),
                    farming = TRUE,
                    unemp = FALSE,
                    w_atc_norm = w_atc_norm) {
  
  level = match.arg(level)
  
  data |>
    mutate(upward =
             (unemp & macro_pop == "unemp" & macro_son != "unemp") |
             (farming & macro_pop == "farming" & macro_son == "white_col") |
             if (level == "macro")
               (macro_pop == "blue_col" & macro_son == "white_col")
           else 
             (dad_meso == "unskilled" & (macro_son == "white_col" | occ_meso == "crafts")) |
             (dad_meso == "clerical" & occ_meso == "prof")) |>
    summarise(upward = w_atc_normed.mean(as.numeric(upward), w = {{ w_atc_norm }})) |>
    pull(upward)
}

county_list = data |>
  group_by(statefip_1940, countyicp_1940) |>
  count() |>
  filter(n >= 100)

compute_mobility_stats = function(df) {
  
  p_upward = df |>
    mutate(count = (macro_pop == "unemp" & macro_pop != "unemp") |
             (macro_pop != "nonmanual" & macro_son == "nonmanual") |
             (meso_pop == "farmworker" & meso_son == "farmer") |
             (meso_pop == "unskilled" & meso_son == "crafts")) |>
    summarise(tot = sum(w_atc_norm), up = sum(count)) |>
    mutate(prop = up / tot) |>
    select(prop)
  
  p_downward = df |>
    mutate(count = (macro_pop != "unemp" & macro_son == "unemp") |
             (macro_pop == "nonmanual" & macro_son != "nonmanual") |
             (meso_pop == "farmer" & meso_son == "farmworker") |
             (meso_pop == "crafts" & meso_son == "unskilled")) |>
    summarise(tot = sum(w_atc_norm), down = sum(count)) |>
    mutate(prop = down / tot) |>
    select(prop)
  
  #d1 = d_prime(df, macro_pop, macro_son, t = 1)
  
  tibble(p_upward = p_upward$prop,
         p_downward = p_downward$prop)
}

results = county_list |>
  mutate(stats = map2(statefip_1940, countyicp_1940, function(st, co) {
    df_sub = data |>
      filter(statefip_1940 == st, countyicp_1940 == co)
    compute_mobility_stats(df_sub)})) |>
  unnest(stats)







  

