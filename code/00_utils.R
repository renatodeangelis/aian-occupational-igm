# Shared utilities: occupation classification, region mapping, weight estimation,
# transition matrix functions, mobility measures, and plot helpers.
# Sourced by all analysis scripts; never run directly.

# --- Occupation classification ---

occ_unclassified = c(975, 976, 977, 978, 979, 995, 997, 999)

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
    occ %in% 500:594 | occ %in% crafts_codes ~ "crafts",
    occ %in% 595:970 & !(occ %in% crafts_codes) & !(occ %in% farmwork_codes) ~ "unskilled",
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

# Canonical display orderings for plots (bottom → top on y-axis).
macro_order         = c("nonemp", "nonmanual", "manual", "farming")
meso_order          = c("nonemp", "nonmanual", "crafts", "unskilled", "farmworker", "farmer")

# Alphabetical order — matches p_matrix() output (which uses sort(union(...))).
# Use when indexing P by name to avoid positional errors.
macro_compute_order = c("farming", "manual", "nonmanual", "nonemp")

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
    summarise(
      birthyr_pop = first(birthyr_pop),
      birthyr_son = first(birthyr_son),
      res = {
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
      },
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

  out |>
    slice_head(n = 1) |>
    ungroup() |>
    transmute(pid, meso = meso_used, year, birthyr_pop, occ_pop)

  stopifnot(all(is.na(out$occ_pop) | classify_meso(out$occ_pop) == out$meso))
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

compute_weights = function(df_linked, df_full,
                           ps_formula = linked ~ cohort * region + education * region + as.factor(urban_1940)) {
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

# Load the global weighted dataset and set macro_levels / meso_levels in the
# calling frame so pi_0() can find them without a warning.
load_global = function(path = "data/aian_weighted.rds") {
  assign("macro_levels", macro_compute_order, envir = parent.frame())
  assign("meso_levels",  meso_order,          envir = parent.frame())
  readRDS(path) |> dplyr::mutate(w_atc_norm = w_trim_norm)
}

# Load per-region weighted datasets. Regional data already has w_atc_norm set
# to trimmed weights by 02_weighting.R — do NOT apply mutate(w_atc_norm = w_trim_norm).
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

boot_pmatrix_ci = function(
    data, level_dad, level_son,
    df_linked, df_full,
    R = 500, .seed = NULL,
    mc.cores = 1L) {

  if (!is.null(.seed)) set.seed(.seed)
  dad_sym = rlang::ensym(level_dad)
  son_sym = rlang::ensym(level_son)
  N = nrow(df_linked)

  w_full = compute_weights(df_linked, df_full)
  d_full = trim_weights_top1(w_full$data)
  P_hat  = p_matrix(d_full, !!dad_sym, !!son_sym, matrix = TRUE)
  rnames = rownames(P_hat); cnames = colnames(P_hat)
  nR = nrow(P_hat); nC = ncol(P_hat)

  boot_once = function() {
    idx = sample.int(N, N, replace = TRUE)
    w_b = compute_weights(df_linked[idx, ], df_full)
    d_b = trim_weights_top1(w_b$data)
    p_matrix(d_b, !!dad_sym, !!son_sym, matrix = TRUE)
  }

  boots = if (mc.cores > 1L) {
    parallel::mclapply(seq_len(R), function(i) boot_once(), mc.cores = mc.cores)
  } else {
    lapply(seq_len(R), function(i) boot_once())
  }
  arr   = simplify2array(boots)

  se_mat = apply(arr, c(1, 2), sd, na.rm = TRUE)

  tibble::tibble(
    !!dad_sym := rep(rnames, times = nC),
    !!son_sym := rep(cnames, each  = nR),
    est = as.vector(P_hat),
    se  = as.vector(se_mat)
  )
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
    d1      = 1 - overlaps[worst],
    row1    = if (!is.null(rn)) rn[ij[1]] else ij[1],
    row2    = if (!is.null(rn)) rn[ij[2]] else ij[2],
    overlap = pmin(P[ij[1], ], P[ij[2], ])
  )
}

# --- Long-format matrix helper ---
# Converts a point-estimate matrix to the same tibble format as boot_pmatrix_ci()
# output (se = NA), so it can be passed directly to plot_pmat().

pmat_long = function(P, dad_nm = "dad", son_nm = "son") {
  rn = rownames(P); cn = colnames(P)
  expand.grid(setNames(list(rn, cn), c(dad_nm, son_nm)),
              stringsAsFactors = FALSE) |>
    dplyr::mutate(est = as.vector(t(P)), se = NA_real_)
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

# Heatmap for a transition matrix. boot_df is the tibble from boot_pmatrix_ci()
# or pmat_long(). If se is non-NA, labels show "est\n(se)"; otherwise just "est".
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
  has_se = !all(is.na(plot_df$se))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = !!son_sym, y = !!dad_sym, fill = est)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = if (has_se)
                     sprintf("%.2f\n(%.3f)", est, se)
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

for (v in c("farm_1940", "classwkr_1940", "labforce_1940", "empstatd_1940",
            "relate_1940", "urban_1940", "gq_1940"))
  print(xtab_son(aian_merged, v), n = 60)
