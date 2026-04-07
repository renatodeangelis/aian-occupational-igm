# Shared classification functions and region mapping
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
    occ > 970 ~ "nilf"
  )
}

classify_macro = function(meso) {
  case_when(
    meso %in% c("farmer", "farmworker") ~ "farming",
    meso == "nonmanual" ~ "nonmanual",
    meso %in% c("crafts", "unskilled") ~ "manual",
    meso == "nilf" ~ "nilf")
}

macro_order = c("farming", "manual", "nonmanual", "nilf")
meso_order = c("farmworker", "farmer", "unskilled", "crafts", "nonmanual", "nilf")
meso_order_alt = c("farmer", "unskilled", "crafts", "nonmanual", "nilf")

# --- Region mapping ---
# Maps statefip codes to 12 regions. Works with both statefip_1940 and statefip columns
# — caller should pass the appropriate variable name.

assign_region = function(statefip) {
  case_when(
    statefip %in% c(8, 16, 32, 49, 56)  ~ "basin",
    statefip == 6                         ~ "cali",
    statefip %in% c(27, 55)               ~ "lakes",
    statefip %in% c(17, 18, 26, 39)       ~ "midwest",
    statefip == 37                        ~ "nc",
    statefip %in% c(9, 10, 23, 24, 25, 33, 34, 36, 42, 44, 50, 11) ~ "ne",
    statefip %in% c(30, 38, 46)           ~ "plains",
    statefip %in% c(41, 53)               ~ "nw",
    statefip == 40                        ~ "ok",
    statefip %in% c(19, 20, 29, 31)       ~ "prairie",
    statefip %in% c(1, 5, 12, 13, 21, 22, 28, 45, 47, 48, 51, 54) ~ "south",
    statefip %in% c(4, 35)                ~ "sw",
    TRUE ~ NA_character_)
}
