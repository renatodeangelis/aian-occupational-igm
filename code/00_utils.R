# Shared utilities: occupation classification, region mapping, weight estimation,
# transition matrix functions, mobility measures, and plot helpers.
# Sourced by all analysis scripts; never run directly.

# --- Occupation classification ---

classify_meso = function(occ) {
  farmer_codes   = c(100, 123, 830)
  farmwork_codes = c(810, 820, 840)
  nonman_codes   = c(0:99, 200:290, 300:490)
  crafts_codes   = c(762, 773, 781, 782)

  base = case_when(
    occ == 979        ~ NA_character_,
    occ %in% farmer_codes   ~ "farmer",
    occ %in% farmwork_codes ~ "farmworker",
    occ %in% nonman_codes   ~ "nonmanual",
    # crafts_codes (762, 773, 781, 782) fall in 595:970; this branch fires first so they land in crafts, not unskilled.
    occ %in% 500:594 | occ %in% crafts_codes ~ "crafts",
    occ %in% 595:970 & !(occ %in% crafts_codes) & !(occ %in% farmwork_codes) ~ "unskilled",
    occ > 970 ~ "nonemp"
  )
  base
}

classify_macro = function(meso) {
  case_when(
    meso %in% c("farmer", "farmworker") ~ "farming",
    meso == "nonmanual" ~ "nonmanual",
    meso %in% c("crafts", "unskilled") ~ "manual",
    meso == "nonemp" ~ "nonemp")
}

# Canonical display orderings for plots (bottom → top on y-axis).
macro_order         = c("nonemp", "nonmanual", "manual", "farming")
meso_order          = c("nonemp", "nonmanual", "crafts", "unskilled", "farmworker", "farmer")


# --- Modal occupation picker ---

pick_modal_meso = function(df, aian_age, prefer_employed = FALSE, empstatd_tiebreak = FALSE) {
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
    mutate(meso = classify_meso(occ)) |>
    group_by(pid, year) |>
    # Groups by pid (father), not pid × son. Where brothers share a pid, birthyr_son
    # is from the first brother (arbitrary); the picked year and meso apply to all siblings.
    summarise(
      birthyr_pop = first(birthyr_pop),
      birthyr_son = first(birthyr_son),
      res = list({
        keep = if (prefer_employed && any(meso != "nonemp", na.rm = TRUE))
          which(!is.na(meso) & meso != "nonemp")
        else
          which(!is.na(meso))
        if (length(keep) == 0) {
          list(meso = NA_character_, occ_pop = NA_integer_)
        } else {
          m = names(which.max(table(meso[keep])))
          list(meso = m, occ_pop = occ[keep][meso[keep] == m][1])
        }
      }),
      .groups = "drop") |>
    tidyr::unnest_wider(res) |>
    mutate(implied_age    = ifelse(!is.na(birthyr_pop), year - birthyr_pop, NA_real_),
           son_age_at_obs = year - birthyr_son) |>
    group_by(pid) |>
    mutate(
      has_pref  = prefer_employed & any(meso != "nonemp", na.rm = TRUE),
      meso_used = if_else(has_pref & meso != "nonemp", meso,
                          if_else(has_pref, NA_character_, meso))) |>
    ungroup() |>
    filter(!is.na(meso_used),
           is.na(implied_age) | implied_age <= 65) |>
    add_count(pid, meso_used, name = "freq") |>
    group_by(pid) |>
    filter(freq == max(freq)) |>
    mutate(has_empstatd = year %in% c(1910, 1930, 1940),
           age_dist     = coalesce(abs(son_age_at_obs - 10), Inf))

  if (empstatd_tiebreak)
    out = arrange(out, pid, desc(has_empstatd), age_dist, year)
  else
    out = arrange(out, pid, age_dist, year)

  result = out |>
    slice_head(n = 1) |>
    ungroup() |>
    transmute(pid, meso = meso_used, year, birthyr_pop, occ_pop)

  stopifnot(all(is.na(result$occ_pop) | classify_meso(result$occ_pop) == result$meso))
  result
}

# --- Region mapping ---

assign_region = function(statefip) {
  case_when(
    statefip %in% c(4, 35)         ~ "sw",
    statefip %in% c(30, 38, 46)    ~ "nplains",
    statefip %in% c(27, 55)        ~ "glakes",
    statefip %in% c(41, 53)        ~ "nw",
    statefip == 40                 ~ "ok",
    statefip == 6                  ~ "cali",
    statefip == 37                 ~ "nc",
    statefip %in% c(8, 16, 32, 49, 56)                        ~ "basin",
    statefip %in% c(19, 20, 29, 31)                           ~ "prairie",
    statefip %in% c(17, 18, 39, 26)                           ~ "midwest",
    statefip %in% c(9, 10, 11, 23, 24, 25, 33, 34, 36, 42, 44, 50) ~ "northeast",
    statefip %in% c(1, 5, 12, 13, 21, 22, 28, 45, 47, 48, 51, 54)  ~ "south",
    TRUE ~ NA_character_)
}

# Regions large enough to estimate: n >= 800 and ESS >= 500 (fixed ex ante).
compare_regions = c("sw", "nplains", "ok", "glakes", "nw", "cali", "nc", "basin")

# Regions excluded from both comparison and benchmark on size grounds.
small_regions   = c("midwest", "prairie", "northeast", "south")

# --- Education classification ---

classify_education = function(educd) {
  case_when(
    educd == 2 ~ "0",
    educd %in% 14:17 ~ "1-4",
    educd %in% 22:26 ~ "5-8",
    educd %in% 30:60 ~ "9-12",
    educd %in% 70:113 ~ "13+",
    educd == 999 ~ "missing")
}

# --- Propensity score weight estimation ---

compute_weights = function(df_linked, df_full,
                           ps_formula = linked ~ cohort * region + education * region + as.factor(urban_1940)) {
  comb = dplyr::bind_rows(
    df_linked |> dplyr::mutate(linked = 1),
    df_full   |> dplyr::mutate(linked = 0)
  ) |>
    dplyr::mutate(
      cohort    = cut(birthyr_son,
                      breaks = c(1890, 1895, 1900, 1905, 1910, 1915, 1921),
                      labels = c("1891-1895", "1896-1900", "1901-1905", "1906-1910",
                                 "1911-1915", "1916-1920")),
      region    = as.factor(region),
      education = as.factor(education))

  model = speedglm::speedglm(ps_formula,
              data = comb, family = binomial())

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
# Kept for appendix robustness checks; the main analysis uses untrimmed weights.
# Trims w_atc at its 99th percentile and renormalises w_atc_norm in place.
#
# ASSUMPTION:
#   df has columns w_atc and w_atc_norm — i.e., it is compute_weights()$data.

trim_weights_top1 = function(df) {
  thresh = quantile(df$w_atc, 0.99, na.rm = TRUE)
  dplyr::mutate(df,
    w_atc      = pmin(w_atc, thresh),
    w_atc_norm = w_atc * dplyr::n() / sum(w_atc)
  )
}

# --- Data loading helpers ---

# Renormalise w_atc_norm within df so weights average to 1.
renorm = function(df) {
  dplyr::mutate(df, w_atc_norm = w_atc_norm / sum(w_atc_norm) * dplyr::n())
}

# Pool a set of regional weighted frames into a single benchmark frame.
# pop_n: optional named vector of target-population counts by region
#        (e.g. table(aian_full$region)). If supplied, each region's weights are
#        rescaled to sum to its target count before pooling (option a: target-
#        population scaling). If NULL, regions enter in proportion to their
#        linked n (option b: linked-sample scaling).
pool_regions = function(regional_data, regions, pop_n = NULL) {
  stopifnot(all(regions %in% names(regional_data)))
  parts = lapply(regions, function(r) {
    d = renorm(regional_data[[r]])
    if (!is.null(pop_n)) {
      stopifnot(r %in% names(pop_n))
      d$w_atc_norm = d$w_atc_norm * (pop_n[[r]] / sum(d$w_atc_norm))
    }
    d
  })
  renorm(dplyr::bind_rows(parts))
}

# Load the global weighted dataset and set macro_levels / meso_levels in the
# calling frame so pi_0() can find them without a warning.
load_global = function(path = "data/aian_weighted.rds") {
  assign("macro_levels", macro_order, envir = parent.frame())
  assign("meso_levels",  meso_order,          envir = parent.frame())
  readRDS(path)
}

# Load per-region weighted datasets.
load_regional = function(path = "data/aian_regional_weighted.rds") {
  readRDS(path)
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

# w_atc_norm must already be normalised in data (mean = 1). Call renorm() after subsetting.
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
  if (all(v <= 0)) {
    v = -v
  } else if (any(v < 0)) {
    stop("pi_star: eigenvector has mixed signs — P may not be ergodic")
  }
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
  pi_init = pi_init[rownames(P_mat)]
  stopifnot(!any(is.na(pi_init)))
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
  pi0 = pi0[rownames(P)]
  stopifnot(!any(is.na(pi0)))
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

# --- Bootstrap for global transition matrices — macro and meso in a single loop ---
#
# Runs compute_weights() once per draw, computes p_matrix() for both level pairs,
# so 2000 draws costs 2000 speedglm fits instead of 4000.
# tryCatch guards against rank-deficient resamples; dimension mismatches are also
# dropped. Discarded draw count is reported; stop() if >5% are lost.
#
# Return value: list with $macro and $meso, each a list with
#   $P  — tibble(dad, son, est, lo, hi)
#   $d1 — named vector c(est, lo, hi)  [log scale, matches dobrushin()$d1]

boot_pmatrix_ci_pair = function(
    data,
    level_dad1, level_son1,
    level_dad2, level_son2,
    df_linked, df_full,
    R = 2000, .seed = NULL,
    mc.cores = 1L,
    pid_col = "pid",
    alpha = 0.05) {

  if (!is.null(.seed)) set.seed(.seed)

  dad_nm1 = rlang::as_string(rlang::ensym(level_dad1))
  son_nm1 = rlang::as_string(rlang::ensym(level_son1))
  dad_nm2 = rlang::as_string(rlang::ensym(level_dad2))
  son_nm2 = rlang::as_string(rlang::ensym(level_son2))

  pid_idx = split(seq_len(nrow(df_linked)), df_linked[[pid_col]])
  fams    = names(pid_idx)

  w_full  = compute_weights(df_linked, df_full)
  d_full  = w_full$data
  P_hat1  = p_matrix(d_full, !!rlang::sym(dad_nm1), !!rlang::sym(son_nm1), matrix = TRUE)
  P_hat2  = p_matrix(d_full, !!rlang::sym(dad_nm2), !!rlang::sym(son_nm2), matrix = TRUE)

  boot_once = function() {
    fam_b = sample(fams, length(fams), replace = TRUE)
    idx   = unlist(pid_idx[fam_b], use.names = FALSE)
    w_b   = tryCatch(
      compute_weights(df_linked[idx, ], df_full),
      error = function(e) NULL)
    if (is.null(w_b)) return(NULL)
    d_b  = w_b$data
    P_b1 = tryCatch(
      p_matrix(d_b, !!rlang::sym(dad_nm1), !!rlang::sym(son_nm1), matrix = TRUE),
      error = function(e) NULL)
    P_b2 = tryCatch(
      p_matrix(d_b, !!rlang::sym(dad_nm2), !!rlang::sym(son_nm2), matrix = TRUE),
      error = function(e) NULL)
    if (is.null(P_b1) || is.null(P_b2)) return(NULL)
    if (!identical(dim(P_b1), dim(P_hat1)) || !identical(dim(P_b2), dim(P_hat2))) return(NULL)
    list(P1   = P_b1,
         P2   = P_b2,
         d1_1 = dobrushin(P_b1)$d1,
         d1_2 = dobrushin(P_b2)$d1)
  }

  boots = if (mc.cores > 1L)
    parallel::mclapply(seq_len(R), function(i) boot_once(), mc.cores = mc.cores)
  else
    lapply(seq_len(R), function(i) boot_once())

  n_null = sum(sapply(boots, is.null))
  if (n_null > 0)
    message(sprintf("bootstrap: %d / %d draws discarded (dimension mismatch or model failure)",
                    n_null, R))
  if (n_null / R > 0.05)
    stop(sprintf("bootstrap: %.1f%% of draws discarded — investigate before proceeding",
                 100 * n_null / R))
  boots = Filter(Negate(is.null), boots)

  lo_p = alpha / 2; hi_p = 1 - alpha / 2

  summarise_pair = function(P_hat, P_key, d1_key, dad_nm, son_nm) {
    rnames = rownames(P_hat); cnames = colnames(P_hat)
    nR = nrow(P_hat); nC = ncol(P_hat)
    arr_P = simplify2array(lapply(boots, `[[`, P_key))
    P_lo = apply(arr_P, c(1, 2), quantile, probs = lo_p, na.rm = TRUE)
    P_hi = apply(arr_P, c(1, 2), quantile, probs = hi_p, na.rm = TRUE)
    d1_draws = sapply(boots, `[[`, d1_key)
    tbl = tibble::tibble(
      dad = rep(rnames, times = nC),
      son = rep(cnames, each  = nR),
      est = as.vector(P_hat),
      lo  = as.vector(P_lo),
      hi  = as.vector(P_hi)
    )
    names(tbl)[1:2] = c(dad_nm, son_nm)
    list(
      P  = tbl,
      d1 = c(est = dobrushin(P_hat)$d1,
             lo  = quantile(d1_draws, lo_p, na.rm = TRUE),
             hi  = quantile(d1_draws, hi_p, na.rm = TRUE))
    )
  }

  list(
    macro = summarise_pair(P_hat1, "P1", "d1_1", dad_nm1, son_nm1),
    meso  = summarise_pair(P_hat2, "P2", "d1_2", dad_nm2, son_nm2)
  )
}

# --- Bootstrap for regional counterfactuals — LOO benchmark rebuilt per draw ---
#
# Resamples region k AND each benchmark region (clustered on pid throughout).
# Cluster indices precomputed once per region to avoid O(N) inner lookups.
# n_na: draws dropped because a benchmark state was absent from k's resample.
#
# Cost warning: 8 regions × R draws, each rebuilding a ~16,000-row benchmark.
# Run R = 200 first to time; R = 1000 is acceptable if 2000 is impractical.

boot_regional_cf = function(
    regional_data,
    level_dad, level_son,
    compare  = compare_regions,
    pop_n    = NULL,
    R        = 2000,
    .seed    = NULL,
    mc.cores = 1L,
    pid_col  = "pid",
    alpha    = 0.05) {

  dad_sym = rlang::ensym(level_dad)
  son_sym = rlang::ensym(level_son)
  if (!is.null(.seed)) set.seed(.seed)

  # Precompute cluster indices once per region — O(N) amortised across all draws
  idx_by_region = lapply(regional_data[compare], function(d)
    split(seq_len(nrow(d)), d[[pid_col]]))

  resample_region = function(r) {
    pid_idx = idx_by_region[[r]]
    fams    = names(pid_idx)
    fam_b   = sample(fams, length(fams), replace = TRUE)
    idx     = unlist(pid_idx[fam_b], use.names = FALSE)
    renorm(regional_data[[r]][idx, ])
  }

  purrr::map_dfr(compare, function(k) {
    bench = setdiff(compare, k)
    d_k   = renorm(regional_data[[k]])
    d_b   = pool_regions(regional_data, bench, pop_n)

    P_b0   = p_matrix(d_b, !!dad_sym, !!son_sym)
    pi0_b0 = pi_0(d_b, !!dad_sym)
    P_k0   = p_matrix(d_k, !!dad_sym, !!son_sym)
    pi0_k0 = pi_0(d_k, !!dad_sym)
    stopifnot(all(rownames(P_b0) %in% rownames(P_k0)))
    P_k0 = P_k0[rownames(P_b0), colnames(P_b0), drop = FALSE]

    regime_hat = sdm1(pi0_k0, P_k0, pi0_k0, P_b0)
    comp_hat   = sdm1(pi0_k0, P_b0, pi0_b0, P_b0)
    total_hat  = sdm1(pi0_k0, P_k0, pi0_b0, P_b0)

    boot_once = function() {
      dk = resample_region(k)
      bparts = lapply(bench, function(r) {
        d = resample_region(r)
        if (!is.null(pop_n)) d$w_atc_norm = d$w_atc_norm * (pop_n[[r]] / sum(d$w_atc_norm))
        d
      })
      db = renorm(dplyr::bind_rows(bparts))

      Pb = p_matrix(db, !!dad_sym, !!son_sym)
      Pk = p_matrix(dk, !!dad_sym, !!son_sym)
      if (!all(rownames(Pb) %in% rownames(Pk)))
        return(c(regime = NA_real_, comp = NA_real_, total = NA_real_, d1_k = NA_real_))
      Pk   = Pk[rownames(Pb), colnames(Pb), drop = FALSE]
      pk   = pi_0(dk, !!dad_sym)
      pb   = pi_0(db, !!dad_sym)
      c(regime = sdm1(pk, Pk, pk, Pb),
        comp   = sdm1(pk, Pb, pb, Pb),
        total  = sdm1(pk, Pk, pb, Pb),
        d1_k   = dobrushin(Pk)$d1)
    }

    boots = if (mc.cores > 1L)
      parallel::mclapply(seq_len(R), function(i) boot_once(), mc.cores = mc.cores)
    else
      lapply(seq_len(R), function(i) boot_once())

    arr  = do.call(rbind, boots)
    lo_p = alpha / 2; hi_p = 1 - alpha / 2
    q    = function(col, p) unname(quantile(arr[, col], p, na.rm = TRUE))

    tibble::tibble(
      region    = k,
      n         = nrow(d_k),
      n_bench   = nrow(d_b),
      n_na      = sum(is.na(arr[, "regime"])),
      regime    = regime_hat,
      regime_lo = q("regime", lo_p), regime_hi = q("regime", hi_p),
      comp      = comp_hat,
      comp_lo   = q("comp",   lo_p), comp_hi   = q("comp",   hi_p),
      total     = total_hat,
      total_lo  = q("total",  lo_p), total_hi  = q("total",  hi_p),
      d1_lo     = q("d1_k",   lo_p), d1_hi     = q("d1_k",   hi_p)
    )
  })
}

# --- Occupation labels and recoding ---

occ_labels = c(
  farming    = "Farming",
  farmer     = "Farming",
  farmworker = "Farmworker",
  nonemp     = "Non-employed",
  nonmanual  = "Non-manual",
  manual     = "Manual",
  crafts     = "Crafts",
  unskilled  = "Unskilled"
)

recode_occ_df = function(df, ...) {
  vars = rlang::ensyms(...)
  for (v in vars) {
    df = df |> mutate(!!v := dplyr::recode(as.character(!!v), !!!occ_labels))
  }
  df
}

recode_occ_vec = function(vec) {
  setNames(as.numeric(vec), dplyr::recode(names(vec), !!!occ_labels))
}

# Display-name level vectors (for axis labels in plots).
# Use these for factor levels on plot axes; use macro_order / meso_order for
# ordering computations and factor() calls on data columns.
macro_display_order = occ_labels[c("nonemp", "nonmanual", "manual", "farming")]
canonical_meso      = c("nonemp", "nonmanual", "crafts", "unskilled", "farmworker", "farmer")
meso_display_order  = occ_labels[canonical_meso]

# --- Dobrushin's contraction coefficient ---
#
# d1 = 1 - min_{i,j} sum_k min(P[i,k], P[j,k])
#    = max TV distance between any two rows of P.
# Equivalent to (1/2) * max_{i,j} sum_k |P[i,k] - P[j,k]|.
# Ranges in [0,1]; lower = faster mixing.

dobrushin = function(P) {
  idx      = combn(nrow(P), 2)
  overlaps = apply(idx, 2, function(ij) sum(pmin(P[ij[1], ], P[ij[2], ])))
  worst    = which.min(overlaps)
  ij       = idx[, worst]
  rn       = rownames(P)
  list(
    d1      = log(1 - overlaps[worst]),
    row1    = if (!is.null(rn)) rn[ij[1]] else ij[1],
    row2    = if (!is.null(rn)) rn[ij[2]] else ij[2],
    overlap = pmin(P[ij[1], ], P[ij[2], ])
  )
}

# --- Pairwise TV distances between rows of P ---
#
# row_dists(P)[i,j] = TV(P[i,], P[j,]).
# Sanity: max(row_dists(P)) == exp(dobrushin(P)$d1) to floating precision.

row_dists = function(P) {
  P = as.matrix(P)
  k = nrow(P)
  out = matrix(NA_real_, k, k, dimnames = list(rownames(P), rownames(P)))
  for (i in seq_len(k)) for (j in seq_len(k)) out[i, j] = tv_norm(P[i, ], P[j, ])
  out
}

# pi0-weighted average pairwise row distance; thin origins receive small weight
# and do not dominate the way they do in the maximand delta(P).
# Not in either anchor paper — define in prose if used.
avg_row_dist = function(P, pi0) {
  pi0 = pi0[rownames(as.matrix(P))]
  stopifnot(!any(is.na(pi0)))
  sum(outer(pi0, pi0) * row_dists(P))
}

# --- Single-period SDM (Blume et al. 2025, eq 7 at t = 1) ---
#
# TV distance between the son distributions generated by two chains:
# (pi0_a, P_a) and (pi0_b, P_b).  States must be identically ordered.

sdm1 = function(pi0_a, P_a, pi0_b, P_b) {
  P_a   = as.matrix(P_a);   P_b   = as.matrix(P_b)
  pi0_a = pi0_a[rownames(P_a)]; pi0_b = pi0_b[rownames(P_b)]
  stopifnot(!any(is.na(pi0_a)), !any(is.na(pi0_b)))
  stopifnot(identical(rownames(P_a), rownames(P_b)))
  tv_norm(as.numeric(pi0_a %*% P_a), as.numeric(pi0_b %*% P_b))
}

# --- Regional counterfactuals — leave-one-out benchmark ---
#
# Benchmark for region k = pool_regions() over all other compare regions.
#   regime_k: given the fathers region k had, how differently did its own regime
#             place sons vs the benchmark regime?
#   comp_k:   how much traces to having different fathers, holding regime at benchmark?
#   total_k:  SDM between region k's actual chain and the benchmark chain.
#
# NOT a decomposition: TV is a norm; regime + comp != total in general.
# Benchmark states must be a subset of region k's states (stopifnot catches gaps).

regional_counterfactuals = function(regional_data,
                                    level_dad, level_son,
                                    compare = compare_regions,
                                    pop_n   = NULL) {
  dad_sym = rlang::ensym(level_dad)
  son_sym = rlang::ensym(level_son)

  purrr::map_dfr(compare, function(k) {
    d_k = renorm(regional_data[[k]])
    d_b = pool_regions(regional_data, setdiff(compare, k), pop_n)

    P_b   = p_matrix(d_b, !!dad_sym, !!son_sym)
    pi0_b = pi_0(d_b, !!dad_sym)
    P_k   = p_matrix(d_k, !!dad_sym, !!son_sym)
    pi0_k = pi_0(d_k, !!dad_sym)

    stopifnot(all(rownames(P_b) %in% rownames(P_k)))
    P_k = P_k[rownames(P_b), colnames(P_b), drop = FALSE]
    stopifnot(identical(dim(P_k), dim(P_b)))

    tibble::tibble(
      region  = k,
      n       = nrow(d_k),
      n_bench = nrow(d_b),
      regime  = sdm1(pi0_k, P_k,   pi0_k,   P_b),
      comp    = sdm1(pi0_k, P_b,   pi0_b,   P_b),
      total   = sdm1(pi0_k, P_k,   pi0_b,   P_b)
    )
  })
}


# --- State FIPS to region lookup (1940 boundaries) ---

state_fips_1940 = tibble::tibble(
  state_name = tolower(c(
    "alabama","arizona","arkansas","california","colorado","connecticut","delaware",
    "florida","georgia","idaho","illinois","indiana","iowa","kansas","kentucky",
    "louisiana","maine","maryland","massachusetts","michigan","minnesota",
    "mississippi","missouri","montana","nebraska","nevada","new hampshire",
    "new jersey","new mexico","new york","north carolina","north dakota","ohio",
    "oklahoma","oregon","pennsylvania","rhode island","south carolina","south dakota",
    "tennessee","texas","utah","vermont","virginia","washington","west virginia",
    "wisconsin","wyoming")),
  statefip = c(
     1, 4, 5, 6, 8, 9,10,
    12,13,16,17,18,19,20,21,
    22,23,24,25,26,27,28,29,
    30,31,32,33,34,35,36,37,
    38,39,40,41,42,44,45,46,
    47,48,49,50,51,53,54,55,56)
) |>
  dplyr::mutate(region = assign_region(statefip))

# --- Plot helpers ---

# Heatmap for a transition matrix. boot_df is the tibble from boot_pmatrix_ci_pair()$macro
# or $meso ($P element). If lo/hi columns are present, labels show "est\n[lo,hi]".
plot_pmat = function(boot_df, dad_var, son_var,
                     levels = NULL, text_size = 5.5, title_expr = "P") {
  dad_sym = rlang::ensym(dad_var)
  son_sym = rlang::ensym(son_var)

  plot_df = boot_df
  if (!is.null(levels)) {
    plot_df = plot_df |>
      dplyr::mutate(!!dad_sym := factor(!!dad_sym, levels = levels),
                    !!son_sym := factor(!!son_sym, levels = rev(levels)))
  }
  has_ci = "lo" %in% names(plot_df) && !all(is.na(plot_df$lo))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = !!son_sym, y = !!dad_sym, fill = est)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = if (has_ci)
                     sprintf("%.2f\n[%.2f,%.2f]", est, lo, hi)
                   else
                     sprintf("%.2f", est)),
      vjust = 0.3, size = text_size) +
    ggplot2::scale_fill_gradient(low = "lightyellow", high = "firebrick",
                                 limits = c(0, 1)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(
      barwidth  = ggplot2::unit(7, "cm"),
      barheight = ggplot2::unit(0.5, "cm"))) +
    ggplot2::labs(x = "Son's occupation", y = NULL, fill = "Prob.",
                  title = title_expr) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1, size = 13),
      axis.text.y     = ggplot2::element_text(angle = 45, hjust = 1, size = 13),
      axis.ticks      = ggplot2::element_blank(),
      axis.title.x    = ggplot2::element_text(size = 14),
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 12),
      legend.title    = ggplot2::element_text(size = 13),
      plot.title      = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
      panel.grid      = ggplot2::element_blank()
    )
}

# Single-column tile for a named probability vector (pi_0 or pi*).
plot_pi = function(vec, title_expr, levels = NULL) {
  df = tibble::tibble(occ = names(vec), value = as.numeric(vec))
  if (!is.null(levels)) {
    df = df |> dplyr::mutate(occ = factor(occ, levels = levels))
  } else {
    df = df |> dplyr::mutate(occ = factor(occ, levels = rev(unique(occ))))
  }

  ggplot2::ggplot(df, ggplot2::aes(x = 1, y = occ, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), size = 5) +
    ggplot2::scale_fill_gradient(low = "lightyellow", high = "firebrick",
                                 limits = c(0, 1)) +
    ggplot2::labs(y = "Father's occupation", title = title_expr) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      axis.title.x    = ggplot2::element_blank(),
      axis.text.x     = ggplot2::element_blank(),
      axis.ticks.x    = ggplot2::element_blank(),
      axis.title.y    = ggplot2::element_text(size = 14),
      axis.text.y     = ggplot2::element_text(angle = 45, hjust = 1, size = 13),
      axis.ticks.y    = ggplot2::element_blank(),
      legend.position = "none",
      plot.title      = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
      panel.grid      = ggplot2::element_blank()
    )
}

pick_at = function(df, stub, year_col = "picked_year") {
  cols = grep(paste0("^", stub, "_pop_\\d{4}$"), names(df), value = TRUE)
  if (!length(cols)) {
    warning("no ", stub, "_pop_* columns found"); return(rep(NA_integer_, nrow(df)))
  }
  yrs = sort(as.integer(sub(".*_pop_", "", cols)))
  m   = as.matrix(df[paste0(stub, "_pop_", yrs)])
  idx = match(df[[year_col]], yrs)
  out = rep(NA_integer_, nrow(df))
  ok  = !is.na(idx)
  out[ok] = as.integer(m[cbind(which(ok), idx[ok])])
  out
}

xtab_pick = function(df, stub, occ_var = "meso_pop", by_year = FALSE) {
  d = df |> mutate(.val = pick_at(df, stub))

  cat("\n=====", stub, "at picked_year — NA coverage =====\n")
  print(d |> group_by(picked_year) |>
          summarise(n = n(), na = sum(is.na(.val)),
                    pct_na = round(100 * mean(is.na(.val)), 1), .groups = "drop"))

  grp = if (by_year) c("picked_year", occ_var) else occ_var

  d |>
    filter(!is.na(.val)) |>
    count(across(all_of(c(grp, ".val")))) |>
    group_by(across(all_of(grp))) |>
    mutate(pct = round(100 * n / sum(n), 1)) |>
    ungroup() |>
    rename(!!stub := .val)
}

xtab_son = function(df, var, occ_var = "meso_son") {
  cat("\n=====", var, "by", occ_var, "=====\n")
  cat("NA:", sum(is.na(df[[var]])), "of", nrow(df), "\n")
  df |>
    filter(!is.na(.data[[var]])) |>
    count(.data[[occ_var]], .data[[var]]) |>
    group_by(.data[[occ_var]]) |>
    mutate(pct = round(100 * n / sum(n), 1)) |>
    ungroup()
}