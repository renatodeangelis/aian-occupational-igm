# Shared classification functions, region mapping, mobility functions,
# and weight estimation.
# Sourced by cleaning-script.R, weighting.R, and transition_matrices_weighted.R

# --- Occupation classification ---

classify_meso = function(occ, split_farmer = TRUE) {
  farmer_codes = c(100, 123)
  nonman_codes = c(1:99, 200:290, 300:490)
  crafts_codes = c(762, 773, 781, 782)

  case_when(
    occ %in% farmer_codes ~ "farmer",
    split_farmer & occ %in% 810:840 ~ "farmworker",
    occ %in% nonman_codes ~ "nonmanual",
    occ %in% 500:594 | occ %in% crafts_codes ~ "crafts",
    occ %in% 595:970 & !(occ %in% crafts_codes) & !(split_farmer & occ %in% 810:840) ~ "unskilled",
    occ > 970 ~ "nonemp"
  )
}

classify_macro = function(meso) {
  case_when(
    meso %in% c("farmer", "farmworker") ~ "farming",
    meso == "nonmanual" ~ "nonmanual",
    meso %in% c("crafts", "unskilled") ~ "manual",
    meso == "nonemp" ~ "nonemp")
}

macro_order = c("farming", "manual", "nonmanual", "nonemp")
meso_order = c("farmworker", "farmer", "unskilled", "crafts", "nonmanual", "nonemp")

# --- Modal occupation picker ---

pick_modal_occ = function(df, aian_age, prefer_employed = TRUE, empstatd_tiebreak = FALSE) {
  out = df |>
    mutate(birthyr_son = 1940 - age_1940) |>
    left_join(select(aian_age, pid, birth_median), by = "pid") |>
    rename(birthyr_pop = birth_median) |>
    select(-starts_with("age")) |>
    filter(birthyr_son > birthyr_pop + 20) |>
    select(pid, starts_with("occ1950_pop_"), birthyr_pop, birthyr_son) |>
    pivot_longer(
      cols = starts_with("occ1950_pop_"),
      names_to = "year",
      names_pattern = "occ1950_pop_(\\d{4})",
      names_transform = list(year = as.integer),
      values_to = "occ",
      values_drop_na = TRUE) |>
    group_by(pid, year) |>
    summarise(
      # Constant within (pid, year) — carried through, not aggregated
      birthyr_pop = first(birthyr_pop),
      birthyr_son = first(birthyr_son),
      occ = {
        pool <- if (prefer_employed && any(occ <= 970)) occ[occ <= 970] else occ
        as.integer(names(which.max(table(pool))))
      },
      .groups = "drop") |>
    mutate(implied_age    = ifelse(!is.na(birthyr_pop), year - birthyr_pop, NA_real_),
           son_age_at_obs = year - birthyr_son) |>
    group_by(pid) |>
    mutate(
      has_pref = prefer_employed & any(occ <= 970, na.rm = TRUE),
      occ_used = if_else(has_pref & occ <= 970, occ,
                         if_else(has_pref, NA_integer_, occ))) |>
    filter(!is.na(occ_used),
           is.na(implied_age) | implied_age <= 65) |>
    add_count(pid, occ_used, name = "freq") |>
    filter(freq == max(freq)) |>
    mutate(has_empstatd = year %in% c(1910, 1930, 1940),
           age_dist     = coalesce(abs(son_age_at_obs - 10), Inf))

  if (empstatd_tiebreak)
    out = arrange(out, pid, desc(has_empstatd), age_dist, year)
  else
    out = arrange(out, pid, age_dist, year)

  out |>
    slice_head(n = 1) |>
    ungroup() |>
    transmute(pid, occ = occ_used, year, birthyr_pop)
}

# --- Region mapping ---

assign_region = function(statefip) {
  case_when(
    statefip == 6 ~ "cali",
    statefip %in% c(27, 55, 17, 18, 26, 39, 9, 10, 23, 24, 25, 33, 34, 36, 42, 44, 50, 11) ~ "north",
    statefip %in% c(8, 16, 32, 49, 56, 41, 53)  ~ "nw",
    statefip == 40 ~ "ok",
    statefip %in% c(19, 20, 29, 31, 30, 38, 46) ~ "plains",
    statefip %in% c(1, 5, 12, 13, 21, 22, 28, 37, 45, 47, 48, 51, 54) ~ "south",
    statefip %in% c(4, 35) ~ "sw",
    TRUE ~ NA_character_)
}

# --- Education classification ---

classify_education = function(educd) {
  case_when(
    educd == 2 ~ "0",
    educd %in% 14:17 ~ "1-4",
    educd %in% 22:26 ~ "5-8",
    educd %in% 30:60 ~ "9-12",
    educd %in% 70:113 ~ "12+",
    educd == 999 ~ "missing")
}

# --- Propensity score weight estimation ---

compute_weights = function(df_linked, df_full) {
  # df_linked: linked father-son pairs (possibly a bootstrap resample)
  # df_full:   full AIAN comparison sample (held fixed)
  # Returns a list:
  #   $data       — df_linked with p_hat, w_atc, w_atc_norm added
  #   $p_hat_full — PS predictions for df_full from the same model

  comb = dplyr::bind_rows(
    df_linked |> dplyr::mutate(linked = 1),
    df_full   |> dplyr::mutate(linked = 0)
  ) |>
    dplyr::mutate(
      cohort    = cut(birthyr_son,
                      breaks = c(1895, 1900, 1905, 1910, 1915, 1921),
                      labels = c("1896-1900", "1901-1905", "1906-1910",
                                 "1911-1915", "1916-1920")),
      region    = as.factor(region),
      education = as.factor(education))

  model = glm(linked ~ cohort * region + education * region + statefip_1940 + urban_1940,
              data = comb, family = binomial)

  comb_linked = dplyr::filter(comb, linked == 1)
  comb_full   = dplyr::filter(comb, linked == 0)

  list(
    data = comb_linked |>
      dplyr::select(-linked) |>
      dplyr::mutate(
        p_hat      = predict(model, newdata = comb_linked, type = "response"),
        w_atc      = (1 - p_hat) / p_hat,
        w_atc_norm = w_atc * dplyr::n() / sum(w_atc)
      ),
    p_hat_full = predict(model, newdata = comb_full, type = "response")
  )
}

# --- Top-1% weight trimmer ---
#
# WHY THIS EXISTS:
#   compute_weights() returns raw ATC weights (w_atc) and a naive normalisation
#   (w_atc_norm). The main analysis pipeline in weighting.R additionally trims
#   extreme weights at the 99th percentile before renormalising. Any bootstrap
#   that re-estimates weights on each resample must apply the same trim, or the
#   distribution of bootstrap draws will be broader than the estimator it is
#   approximating — inflating SEs.
#
# DESIGN CHOICES:
#   - Trim is applied to w_atc (raw, pre-normalisation), matching weighting.R
#     lines 69-75. Trimming w_atc_norm instead would shift the threshold each
#     time the sample size changes, making resamples non-comparable.
#   - The 99th-percentile threshold is recomputed per call (i.e., per resample)
#     rather than being fixed to the full-sample threshold. This is correct:
#     each draw has its own weight distribution, and the trim should be
#     calibrated to that draw's distribution.
#   - The function overwrites w_atc_norm in place so downstream functions
#     (p_matrix, pi_0) pick up trimmed weights without any argument changes.
#
# ASSUMPTION:
#   df has columns w_atc and w_atc_norm — i.e., it is compute_weights()$data.
#   Calling this on any other data frame will silently produce wrong results.

trim_weights_top1 = function(df) {
  thresh = quantile(df$w_atc, 0.99, na.rm = TRUE)
  dplyr::mutate(df,
    w_atc      = pmin(w_atc, thresh),    # cap extremes; do not remove rows
    w_atc_norm = w_atc * dplyr::n() / sum(w_atc)   # renormalise after cap
  )
}

# --- Transition matrix and distribution functions ---

pi_0 = function(data, level) {
  level_sym = rlang::ensym(level)
  level_nm  = rlang::as_string(level_sym)

  if (level_nm == "macro_pop") {
    if (exists("macro_levels", where = parent.frame(), inherits = TRUE)) {
      all_levels = get("macro_levels", envir = parent.frame())
    } else {
      all_levels = unique(data[[level_nm]])
      warning("macro_levels not found; falling back to unique(data$macro_pop)")
    }
  } else if (level_nm == "meso_pop") {
    if (exists("meso_levels", where = parent.frame(), inherits = TRUE)) {
      all_levels = get("meso_levels", envir = parent.frame())
    } else {
      all_levels = unique(data[[level_nm]])
      warning("meso_levels not found; falling back to unique(data$meso_pop)")
    }
  } else {
    stop("`level` must be either `macro_pop` or `meso_pop`.")
  }

  df = data |>
    dplyr::group_by(!!level_sym) |>
    dplyr::summarise(total_w = sum(w_atc_norm), .groups = "drop") |>
    tidyr::complete(!!level_sym := all_levels, fill = list(total_w = 0)) |>
    dplyr::mutate(pi0 = total_w / sum(total_w)) |>
    dplyr::arrange(factor(!!level_sym, levels = all_levels))

  pi0_vec = df$pi0
  names(pi0_vec) = df[[level_nm]]
  return(pi0_vec)
}

pi_0_unweighted = function(data, level) {
  level_sym = rlang::ensym(level)
  level_nm  = rlang::as_string(level_sym)

  if (level_nm == "macro_pop") {
    if (exists("macro_levels", where = parent.frame(), inherits = TRUE)) {
      all_levels = get("macro_levels", envir = parent.frame())
    } else {
      all_levels = unique(data[[level_nm]])
      warning("macro_levels not found; falling back to unique(data$macro_pop)")
    }
  } else if (level_nm == "meso_pop") {
    if (exists("meso_levels", where = parent.frame(), inherits = TRUE)) {
      all_levels = get("meso_levels", envir = parent.frame())
    } else {
      all_levels = unique(data[[level_nm]])
      warning("meso_levels not found; falling back to unique(data$meso_pop)")
    }
  } else {
    stop("`level` must be either `macro_pop` or `meso_pop`.")
  }

  df = data |>
    dplyr::group_by(!!level_sym) |>
    dplyr::summarise(total_n = dplyr::n(), .groups = "drop") |>
    tidyr::complete(!!level_sym := all_levels, fill = list(total_n = 0)) |>
    dplyr::mutate(pi0 = total_n / sum(total_n)) |>
    dplyr::arrange(factor(!!level_sym, levels = all_levels))

  pi0_vec = df$pi0
  names(pi0_vec) = df[[level_nm]]
  return(pi0_vec)
}

p_matrix = function(data, level_dad, level_son, matrix = TRUE) {
  dad_nm = rlang::as_string(rlang::ensym(level_dad))
  son_nm = rlang::as_string(rlang::ensym(level_son))

  all_levels = sort(union(unique(data[[dad_nm]]), unique(data[[son_nm]])))
  dad_f = factor(data[[dad_nm]], levels = all_levels)
  son_f = factor(data[[son_nm]], levels = all_levels)

  tab = xtabs(data$w_atc_norm ~ dad_f + son_f)
  rs  = rowSums(tab)
  P   = sweep(tab, 1, ifelse(rs > 0, rs, 1), "/")

  mat = base::matrix(as.numeric(P), nrow(P), ncol(P),
                     dimnames = list(all_levels, all_levels))

  if (!matrix) {
    return(
      expand.grid(setNames(list(all_levels, all_levels), c(dad_nm, son_nm)),
                  stringsAsFactors = FALSE) |>
        dplyr::mutate(P = as.vector(t(mat)))
    )
  }

  mat
}

p_matrix_unweighted = function(data, level_dad, level_son, matrix = TRUE) {
  dad_sym = rlang::ensym(level_dad)
  son_sym = rlang::ensym(level_son)

  dad_nm  = rlang::as_string(dad_sym)
  son_nm  = rlang::as_string(son_sym)

  all_levels = sort(union(unique(data[[dad_nm]]), unique(data[[son_nm]])))

  complete_df = tidyr::expand_grid(
    !!dad_sym := all_levels,
    !!son_sym := all_levels
  )

  df = data |>
    dplyr::group_by(!!dad_sym, !!son_sym) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  df = dplyr::right_join(complete_df, df, by = c(dad_nm, son_nm)) |>
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) |>
    dplyr::group_by(!!dad_sym) |>
    dplyr::mutate(
      n_total = sum(n),
      P = ifelse(n_total > 0, n / n_total, 0)
    ) |>
    dplyr::ungroup()

  if (!matrix) return(df)

  wide = df |>
    dplyr::select(!!dad_sym, !!son_sym, P) |>
    tidyr::pivot_wider(
      names_from  = !!son_sym,
      values_from = P,
      values_fill = list(P = 0)
    )

  missing_rows = setdiff(all_levels, wide[[dad_nm]])
  if (length(missing_rows) > 0) {
    extra = data.frame(matrix(0, nrow = length(missing_rows), ncol = ncol(wide)))
    colnames(extra) = colnames(wide)
    extra[[dad_nm]] = missing_rows
    wide = rbind(wide, extra)
  }

  wide = wide[match(all_levels, wide[[dad_nm]]), ]
  mat = as.matrix(wide |> dplyr::select(-!!dad_sym))
  rownames(mat) = dplyr::pull(wide, !!dad_sym)

  return(mat)
}

verify_ergodic = function(P, label = NULL) {
  eigs = abs(Re(eigen(t(as.matrix(P)))$values))
  n_unit = sum(abs(eigs - 1) < 1e-8)
  if (n_unit != 1) {
    msg = sprintf("P has %d unit eigenvalues; stationary distribution is not unique.", n_unit)
    if (!is.null(label)) msg = paste0("[", label, "] ", msg)
    warning(msg)
  }
  invisible(n_unit == 1)
}

pi_star = function(p_mat) {
  P = as.matrix(p_mat)
  eig = eigen(t(P))
  idx = which.min(abs(eig$values - 1))
  v = Re(eig$vectors[, idx])
  if (any(v < 0)) v = abs(v)
  pi_s = v / sum(v)
  names(pi_s) = rownames(P)
  return(pi_s)
}

tv_norm = function(mu, nu) {
  0.5 * sum(abs(mu - nu))
}

# --- Mobility measures ---

d_t = function(data, level_dad, level_son, t = 1) {
  P_mat = p_matrix(data, {{ level_dad }}, {{ level_son }})
  pi_s  = pi_star(P_mat)
  P_t = P_mat %^% t
  d_i = apply(P_t, 1, function(row_i) tv_norm(row_i, pi_s))
  log(max(d_i))
}

d_prime = function(data, level_dad, level_son, t = 1) {
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

am = function(data, level_dad, level_son, t = 1) {
  pi_init = pi_0(data, {{ level_dad }})
  P_mat   = p_matrix(data, {{ level_dad }}, {{ level_son }})
  pi_s    = pi_star(P_mat)
  P_t  = P_mat %^% t
  pi_t = as.numeric(pi_init %*% P_t)
  log(tv_norm(pi_t, pi_s))
}

im = function(data, level_dad, level_son, t = 1) {
  P_mat = p_matrix(data, {{ level_dad }}, {{ level_son }})
  pi_s  = pi_star(P_mat)
  P_t  = P_mat %^% t
  im_i = apply(P_t, 1, function(row_i) tv_norm(row_i, pi_s))
  log(im_i)
}

mu_t = function(pi0, P, t = 0) {
  P = as.matrix(P)
  if (t == 0) return(as.numeric(pi0))
  as.numeric(pi0 %*% (P %^% t))
}

om = function(P, pi0, t) {
  mu = mu_t(pi0, P, t)
  1 - sum(mu * diag(P))
}

sm = function(P, pi0, t) {
  mu  = mu_t(pi0, P, t)
  mu1 = as.numeric(mu %*% P)
  tv_norm(mu, mu1)
}

# --- Generator identification helpers ---

d_generator = function(P_t, pi_star) {
  scores = apply(P_t, 1, function(r) tv_norm(r, pi_star))
  mx = max(scores)
  i_star = which(abs(scores - mx) < 1e-12)
  list(classes = rownames(P_t)[i_star], value = mx)
}

dprime_generator = function(P_t) {
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

# --- Bootstrap functions ---

boot_measures_by_t = function(data, level_dad, level_son,
                               ts = 0:4, R = 1000, .seed = NULL) {
  if (!is.null(.seed)) set.seed(.seed)
  dad_sym = rlang::ensym(level_dad)
  son_sym = rlang::ensym(level_son)

  compute_all = function(P, pi0, pi_s) {
    Pt = diag(nrow(P))
    purrr::map_dfr(ts, function(tt) {
      if (tt > 0) Pt <<- Pt %*% P
      pi_t   = as.numeric(pi0 %*% Pt)
      n      = nrow(Pt)
      pairs  = combn(n, 2)
      tibble::tibble(
        t       = tt,
        d       = log(max(apply(Pt, 1, function(r) tv_norm(r, pi_s)))),
        d_prime = log(max(apply(pairs, 2, function(i) tv_norm(Pt[i[1],], Pt[i[2],])))),
        AM      = log(tv_norm(pi_t, pi_s))
      )
    })
  }

  # point estimates on the full data
  P_hat   = p_matrix(data, !!dad_sym, !!son_sym, matrix = TRUE)
  pi0_hat = as.numeric(pi_0(data, !!dad_sym))
  piS_hat = pi_star(P_hat)
  point_df = compute_all(P_hat, pi0_hat, piS_hat)

  # bootstrap
  N = nrow(data)
  boots = replicate(R, {
    idx   = sample.int(N, N, replace = TRUE)
    P_b   = p_matrix(data[idx, ], !!dad_sym, !!son_sym, matrix = TRUE)
    pi0_b = as.numeric(pi_0(data[idx, ], !!dad_sym))
    piS_b = pi_star(P_b)
    as.matrix(compute_all(P_b, pi0_b, piS_b)[, c("d", "d_prime", "AM")])
  }, simplify = FALSE)

  arr = simplify2array(boots)

  measure_nms = c("d", "d_prime", "AM")
  purrr::map_dfr(seq_along(ts), function(i_t) {
    purrr::map_dfr(seq_along(measure_nms), function(m_idx) {
      v = arr[i_t, m_idx, ]
      tibble::tibble(
        measure = measure_nms[m_idx],
        t       = ts[i_t],
        est     = point_df[[measure_nms[m_idx]]][i_t],
        se      = sd(v, na.rm = TRUE),
        lo      = quantile(v, 0.025, na.rm = TRUE, names = FALSE),
        hi      = quantile(v, 0.975, na.rm = TRUE, names = FALSE)
      )
    })
  }) |> dplyr::mutate(dt_t = est * t)
}

# --- Bootstrap SE for transition matrix cells ---
#
# WHY THIS REPLACES THE OLD boot_pmatrix_ci:
#   The previous version resampled the pre-weighted data object and reused fixed
#   weights. This understates SEs because propensity score estimation uncertainty
#   is not propagated. The fix is to re-run compute_weights() on every resample,
#   then trim, so each draw reflects the full estimator including weight
#   uncertainty.
#
# SCOPE:
#   Returns standard errors only — not confidence intervals. SEs are the only
#   bootstrap output used for the transition matrix tables. Percentile CIs would
#   require R ≥ 1000 and are not needed here; sd() across draws converges faster.
#
# INPUTS:
#   data      — the loaded analysis dataset (aian_weighted.csv). Used only for
#               the occupational classification columns (macro_pop, macro_son,
#               etc.); existing weight columns are overwritten by each call to
#               compute_weights().
#   df_linked — same object as data is fine. compute_weights() only reads
#               birthyr_son, region, education, statefip_1940, urban_1940 for
#               the PS model; all other columns pass through unchanged.
#   df_full   — the full AIAN extract (aian_full.rds). Held fixed across all
#               draws. Resampling df_full as well would be defensible but is
#               not standard practice for ATC weighting where df_full represents
#               a (near-)population target.
#
# ASSUMPTIONS:
#   1. df_linked has region and education columns (added by weighting.R before
#      compute_weights() was called). Passing aian_weighted.csv satisfies this.
#   2. Occasional bootstrap resamples may produce sparse cohort×region cells,
#      causing the GLM to fail or return extreme predictions. These draws are
#      not guarded against here; consider wrapping boot_once() in tryCatch()
#      if convergence warnings appear in practice.
#   3. R = 500 is sufficient for stable SE estimation. For the final paper,
#      bump to 1000 and verify SEs change by < 5%.

boot_pmatrix_ci = function(
    data, level_dad, level_son,
    df_linked, df_full,
    R = 500, .seed = NULL) {

  if (!is.null(.seed)) set.seed(.seed)
  dad_sym = rlang::ensym(level_dad)
  son_sym = rlang::ensym(level_son)
  N = nrow(df_linked)

  # Point estimate: run the full pipeline (weight → trim → P) on the complete
  # linked sample, so the point estimate is on the same pipeline as each draw.
  w_full = compute_weights(df_linked, df_full)
  d_full = trim_weights_top1(w_full$data)
  P_hat  = p_matrix(d_full, !!dad_sym, !!son_sym, matrix = TRUE)
  rnames = rownames(P_hat); cnames = colnames(P_hat)
  nR = nrow(P_hat); nC = ncol(P_hat)

  boot_once = function() {
    idx = sample.int(N, N, replace = TRUE)
    w_b = compute_weights(df_linked[idx, ], df_full)   # re-estimate PS on draw
    d_b = trim_weights_top1(w_b$data)                  # trim this draw's weights
    p_matrix(d_b, !!dad_sym, !!son_sym, matrix = TRUE)
  }

  boots = replicate(R, boot_once(), simplify = FALSE)
  arr   = simplify2array(boots)    # nR × nC × R array

  se_mat = apply(arr, c(1, 2), sd, na.rm = TRUE)

  tibble::tibble(
    !!dad_sym := rep(rnames, times = nC),
    !!son_sym := rep(cnames, each  = nR),
    est = as.vector(P_hat),
    se  = as.vector(se_mat)
  )
}

boot_im_by_t = function(data, level_dad, level_son, ts = 0:4, R = 1000, .seed = NULL) {
  if (!is.null(.seed)) set.seed(.seed)

  dad_sym = rlang::ensym(level_dad)
  son_sym = rlang::ensym(level_son)

  P0   = p_matrix(data, !!dad_sym, !!son_sym, matrix = TRUE)
  piS0 = pi_star(P0)
  rowlabs = rownames(P0)

  Pt0 = diag(nrow(P0))
  point_by_t = purrr::map(ts, function(tt) {
    if (tt > 0) Pt0 <<- Pt0 %*% P0
    im_vals = apply(Pt0, 1, function(r) tv_norm(r, piS0))
    tibble::tibble(t = tt, origin = rowlabs, est = log(im_vals))
  }) |> dplyr::bind_rows()

  N = nrow(data)
  boot_once = function() {
    idx = sample.int(N, N, replace = TRUE)
    db  = data[idx, , drop = FALSE]
    P   = p_matrix(db, !!dad_sym, !!son_sym, matrix = TRUE)
    piS = pi_star(P)

    Pt = diag(nrow(P))
    out_list = vector("list", length(ts))
    for (k in seq_along(ts)) {
      if (ts[k] > 0) Pt = Pt %*% P
      im_vals = apply(Pt, 1, function(r) tv_norm(r, piS))
      out_list[[k]] = log(im_vals)
    }
    do.call(cbind, out_list)
  }

  boots = replicate(R, boot_once(), simplify = FALSE)
  arr = simplify2array(boots)
  dimnames(arr) = list(origin = rowlabs, t = as.character(ts), rep = NULL)

  alpha = 0.05
  summ = lapply(ts, function(tt) {
    a2 = arr[, as.character(tt), , drop = FALSE]
    draws_mat = drop(a2)
    se = apply(draws_mat, 1, sd, na.rm = TRUE)
    lo = apply(draws_mat, 1, quantile, probs = alpha/2, na.rm = TRUE, names = FALSE)
    hi = apply(draws_mat, 1, quantile, probs = 1 - alpha/2, na.rm = TRUE, names = FALSE)
    tibble::tibble(t = tt, origin = rowlabs, se = unname(se), lo = unname(lo), hi = unname(hi))
  }) |> dplyr::bind_rows()

  dplyr::left_join(point_by_t, summ, by = c("t", "origin"))
}

# --- Bootstrap SE for EM/SM mobility curves ---
#
# WHY THIS REPLACES THE OLD mobility_curve_with_boot:
#   Same root problem as boot_pmatrix_ci: the previous version resampled the
#   pre-weighted dataset without re-estimating propensity scores, understating
#   SEs — especially for EM and SM, which depend on the marginal distributions
#   pi_0 and pi_star and are therefore sensitive to weight variance.
#
# SCOPE:
#   Returns SEs for EM and SM only. OM is excluded: the paper's analytical
#   contribution is the EM/SM decomposition, and OM is recoverable as EM + SM
#   at the caller level if needed. Dropping OM halves the output width and
#   avoids implying OM has independent inferential content here.
#   No CIs are returned — same reasoning as boot_pmatrix_ci.
#
# INPUTS:
#   data, df_linked, df_full — same roles as in boot_pmatrix_ci (see above).
#
# DESIGN CHOICES:
#   - compute_em_sm() is a closure that captures dad_sym/son_sym. It is defined
#     inside the outer function so the bootstrap loop can call it without
#     passing symbols explicitly.
#   - pi_0() uses <<- to look up macro_levels/meso_levels in the calling
#     environment. The closure keeps those lookups in scope.
#   - The pre-generated index matrix (inds_mat) from the previous version is
#     replaced with a simpler replicate() loop. The index matrix was a
#     micro-optimisation that added complexity without measurable speed gain
#     at R = 500.
#
# ASSUMPTIONS (shared with boot_pmatrix_ci):
#   1. df_linked already has region and education columns.
#   2. Sparse resample GLM failures are not guarded against; monitor for
#      convergence warnings on the first run.
#   3. R = 500 is sufficient for SE estimation.

mobility_curve_with_boot = function(
    data, level_dad, level_son,
    df_linked, df_full,
    ts = 1:5, R = 500, .seed = NULL) {

  if (!is.null(.seed)) set.seed(.seed)
  dad_sym = rlang::ensym(level_dad)
  son_sym = rlang::ensym(level_son)
  N = nrow(df_linked)

  # Inner stat function: given a weighted data frame, return EM and SM for
  # each t. Called once for the point estimate and once per bootstrap draw.
  compute_em_sm = function(df) {
    P   = p_matrix(df, !!dad_sym, !!son_sym, matrix = TRUE)
    pi0 = as.numeric(pi_0(df, !!dad_sym))

    purrr::map_dfr(ts, function(tt) {
      # mu is the marginal distribution of fathers' occupations at generation t.
      # At t = 0 it equals pi_0; for t > 0 it advances one step through P.
      mu  = as.numeric(pi0 %*% (if (tt == 0) diag(nrow(P)) else P %^% tt))
      mu1 = as.numeric(mu %*% P)           # one step forward from mu
      om_v = 1 - sum(mu * diag(P))         # overall mobility at t
      sm_v = tv_norm(mu, mu1)              # structural component
      tibble::tibble(t = tt, EM = om_v - sm_v, SM = sm_v)
    })
  }

  # Point estimate: run full pipeline on complete linked sample.
  w_full   = compute_weights(df_linked, df_full)
  d_full   = trim_weights_top1(w_full$data)
  point_df = compute_em_sm(d_full)

  boot_once = function() {
    idx = sample.int(N, N, replace = TRUE)
    w_b = compute_weights(df_linked[idx, ], df_full)
    d_b = trim_weights_top1(w_b$data)
    as.matrix(compute_em_sm(d_b)[, c("EM", "SM")])   # nT × 2
  }

  boots = replicate(R, boot_once(), simplify = FALSE)
  arr   = simplify2array(boots)   # nT × 2 × R

  purrr::map_dfr(seq_along(ts), function(i_t) {
    # draws is a 2 × R matrix; rows correspond to EM and SM respectively.
    draws = arr[i_t, , , drop = FALSE]
    draws = matrix(draws, nrow = 2, ncol = R)
    tibble::tibble(
      t       = ts[i_t],
      measure = c("EM", "SM"),
      est     = as.numeric(point_df[i_t, c("EM", "SM")]),
      se      = apply(draws, 1, sd, na.rm = TRUE)
    )
  })
}
