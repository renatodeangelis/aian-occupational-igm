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
  # Source: 08_regional_table.R assertions (global data filtered by region)
  # Fields: farm_ret, pi_farming, pi_manual, pi_nonman, lambda2, relief, n
  regional = list(
    sw     = list(farm_ret=.623, pi_farming=.469, pi_manual=.386, pi_nonman=.032,
                  lambda2=.313, relief=.110, n=2090L),
    south  = list(farm_ret=.739, pi_farming=.376, pi_manual=.482, pi_nonman=.072,
                  lambda2=.682, relief=.062, n=1139L),
    cali   = list(farm_ret=.442, pi_farming=.350, pi_manual=.513, pi_nonman=.016,
                  lambda2=.186, relief=.165, n= 645L),
    ok     = list(farm_ret=.459, pi_farming=.276, pi_manual=.405, pi_nonman=.113,
                  lambda2=.377, relief=.136, n=2531L),
    plains = list(farm_ret=.370, pi_farming=.268, pi_manual=.538, pi_nonman=.067,
                  lambda2=.170, relief=.255, n=2656L),
    nw     = list(farm_ret=.436, pi_farming=.243, pi_manual=.595, pi_nonman=.031,
                  lambda2=.295, relief=.174, n=1155L),
    north  = list(farm_ret=.234, pi_farming=.106, pi_manual=.698, pi_nonman=.061,
                  lambda2=.180, relief=.233, n=2030L)
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

  # --- Global macro mobility scalars at t=1 (fill after re-run) ---
  macro_global = list(
    om_1  = NA_real_,
    sm_1  = NA_real_,
    em_1  = NA_real_,
    d1    = NA_real_,
    n     = NA_integer_
  ),

  # --- Global meso mobility scalars at t=1 (fill after re-run) ---
  meso_global = list(
    om_1  = NA_real_,
    sm_1  = NA_real_,
    em_1  = NA_real_,
    d1    = NA_real_,
    n     = NA_integer_
  )

)
