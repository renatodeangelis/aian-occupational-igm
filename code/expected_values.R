# expected_values.R
# Tolerance and expected values for assertion checks across analysis scripts.
# Numeric values come from validated runs; fill NA entries after first clean run.
# Scripts source this file and check abs(got - expected) < TOL.

TOL = 0.002

EXPECTED = list(

  # --- Global macro pi_0 (w_trim_norm) ---
  # Source: 08_regional_table.R national check; compute_order = alphabetical
  pi0_global_macro = c(
    farming   = .659,
    manual    = .233,
    nonmanual = .046,
    nonemp    = .062
  ),

  # --- Regional macro statistics ---
  # Source: 03_estimate.R (region-specific PS model via load_regional())
  # Fields: farm_ret, pi_farming, pi_manual, pi_nonman, lambda2, relief, n
  # Fill NA entries after first clean run of 03_estimate.R.
  regional = list(
    sw     = list(farm_ret=NA_real_, pi_farming=NA_real_, pi_manual=NA_real_, pi_nonman=NA_real_,
                  lambda2=NA_real_, relief=NA_real_, n=NA_integer_),
    south  = list(farm_ret=NA_real_, pi_farming=NA_real_, pi_manual=NA_real_, pi_nonman=NA_real_,
                  lambda2=NA_real_, relief=NA_real_, n=NA_integer_),
    cali   = list(farm_ret=NA_real_, pi_farming=NA_real_, pi_manual=NA_real_, pi_nonman=NA_real_,
                  lambda2=NA_real_, relief=NA_real_, n=NA_integer_),
    ok     = list(farm_ret=NA_real_, pi_farming=NA_real_, pi_manual=NA_real_, pi_nonman=NA_real_,
                  lambda2=NA_real_, relief=NA_real_, n=NA_integer_),
    plains = list(farm_ret=NA_real_, pi_farming=NA_real_, pi_manual=NA_real_, pi_nonman=NA_real_,
                  lambda2=NA_real_, relief=NA_real_, n=NA_integer_),
    nw     = list(farm_ret=NA_real_, pi_farming=NA_real_, pi_manual=NA_real_, pi_nonman=NA_real_,
                  lambda2=NA_real_, relief=NA_real_, n=NA_integer_),
    north  = list(farm_ret=NA_real_, pi_farming=NA_real_, pi_manual=NA_real_, pi_nonman=NA_real_,
                  lambda2=NA_real_, relief=NA_real_, n=NA_integer_)
  ),

  # --- Slide 7 panel 1: farming father → macro_son exit shares ---
  # Source: 07_slide7_panels.R assertions
  panel1 = c(
    farming   = 0.507,
    manual    = 0.349,
    nonemp    = 0.109,
    nonmanual = 0.035
  ),

  # --- Slide 7 panel 2: farming→manual sons → empstatd_1940 ---
  # Source: 07_slide7_panels.R assertions
  panel2 = c(
    "10"    = 0.425,
    "11"    = 0.387,
    "21"    = 0.135,
    other   = 0.053
  ),

  # --- Global macro mobility scalars at t=0 (fill after re-run) ---
  macro_global = list(
    om_0  = NA_real_,
    sm_0  = NA_real_,
    em_0  = NA_real_,
    d1    = NA_real_,
    n     = NA_integer_
  ),

  # --- Global meso mobility scalars at t=0 (fill after re-run) ---
  meso_global = list(
    om_0  = NA_real_,
    sm_0  = NA_real_,
    em_0  = NA_real_,
    d1    = NA_real_,
    n     = NA_integer_
  )

)
