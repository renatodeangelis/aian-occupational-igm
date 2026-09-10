# 00_extract_diagnostic.R
#
# PURPOSE: decide whether to rebuild on MLP v2.0. This is NOT an analysis
# script and produces no analysis objects. It answers three questions and
# stops:
#
#   Q0  Is the extract valid? (did per-year case selection destroy it)
#   Q1  How many father-son pairs does MLP v2 alone recover, vs 12,246?
#   Q2  What share of sons coded AIAN in 1940 were coded otherwise earlier?
#   Q3  Does POPLOC recover fathers cleanly enough to retire the hand-rolled
#       reconstitution in 01_cleaning-script.R?
#
# Run it, read the VERDICT block at the bottom, then decide. Do not extend
# this file into a pipeline.

library(dplyr)
library(readr)
library(tidyr)
library(janitor)

# IPUMS linked extracts are LONG: all person records for census year X,
# then all person records for year Y, etc. The existing pipeline assumes
# WIDE (_1910, _1920 suffixes). Nothing in 01_cleaning-script.R transfers
# directly; the reshape IS the reconstitution.

EXTRACT_PATH = "https://www.dropbox.com/scl/fi/q4725zk5qw3ltruiucwgc/usa_00026.csv?rlkey=vzauffbwyj61upzhb4fezbq5c&st=y9467v1h&dl=1"
AIAN         = 3                              # IPUMS RACE code
OLD_N        = 12246                          # current sample, CLP u MLP v1.0

raw = read_csv(EXTRACT_PATH, col_types = cols(.default = col_character())) |>
  clean_names() |>
  type_convert(col_types = cols(histid = col_character(),
                                hik = col_character(),
                              .default = col_guess()))

# If this file is >~5 GB, swap to:
#   arrow::open_dataset(EXTRACT_PATH, format = "csv") |> filter(...) |> collect()

################################################################################
############################ 0. INVENTORY ######################################
################################################################################

cat("\n================ INVENTORY ================\n")
cat("Total records:", nrow(raw), "\n")
cat("Columns:", paste(sort(names(raw)), collapse = ", "), "\n\n")

print(raw |> count(year, name = "records"))

cat("\nUnique persons per year (HISTID):\n")
print(raw |> group_by(year) |> summarise(n_histid = n_distinct(histid)))

# Household members should have no HIK. If every record has one, the
# household-member option did not take and no fathers are recoverable.
cat("\nHIK coverage by year (linked persons vs household members):\n")
print(raw |> group_by(year) |>
        summarise(has_hik = sum(!is.na(hik) & hik != ""),
                  no_hik  = sum(is.na(hik) | hik == ""),
                  pct_hik = round(100 * mean(!is.na(hik) & hik != ""), 1)))

stopifnot("POPLOC missing - re-request extract" = "poploc" %in% names(raw))

################################################################################
##################### Q0. IS THE EXTRACT VALID? ################################
################################################################################
# If case selection applied per census year, every record in every year is
# AIAN, reclassifiers are gone, and Q2 is unanswerable from this file.

cat("\n================ Q0: EXTRACT VALIDITY ================\n")
race_by_year = raw |> count(year, race) |> group_by(year) |>
  mutate(pct = round(100 * n / sum(n), 1)) |> ungroup()
print(race_by_year)

nonaian_pre1940 = race_by_year |> filter(year < 1940, race != AIAN) |> pull(n) |> sum()

if (nonaian_pre1940 == 0) {
  cat("\n*** WARNING: zero non-AIAN records before 1940.\n")
  cat("*** Either the race filter applied per-year (extract is compromised for\n")
  cat("*** Q2 and understates Q1), or household members were excluded.\n")
  cat("*** Cross-check against the HIK coverage table above before proceeding.\n")
} else {
  cat("\nNon-AIAN records before 1940:", nonaian_pre1940,
      "- filter did not apply per-year. Q2 is answerable.\n")
}

################################################################################
###################### Q1. FATHER-SON RECOVERY #################################
################################################################################

sons = raw |>
  filter(year == 1940, sex == 1, between(age, 20, 44), race == AIAN,
         !is.na(hik), hik != "") |>
  select(son_hik = hik, son_histid = histid, son_age = age,
         son_occ = occ1950, son_statefip = statefip)

cat("\n================ Q1: FATHER-SON RECOVERY ================\n")
cat("Sons: AIAN males 20-44 in 1940 with a HIK:", nrow(sons), "\n")

# Every record belonging to a son, in every year he appears (incl. 1940).
son_records = raw |>
  filter(hik %in% sons$son_hik) |>
  select(son_hik = hik, year, serial, pernum, poploc, relate, age, race)

cat("Son-year records recovered:\n")
print(son_records |> count(year, name = "son_records"))

# Father = person at PERNUM == POPLOC in the same (year, serial).
father_records = son_records |>
  filter(poploc > 0) |>
  inner_join(raw |> select(year, serial, pernum, dad_histid = histid,
                           dad_hik = hik, dad_age = age, dad_race = race,
                           dad_occ = occ1950, dad_relate = relate),
             by = c("year", "serial", "poploc" = "pernum"))

cat("\nFather observations by year (via POPLOC):\n")
print(father_records |> count(year, name = "father_obs"))

pairs_poploc = n_distinct(father_records$son_hik)

cat("\nSons with >= 1 father observation:", pairs_poploc,
    sprintf("(%.1f%% of sons)\n", 100 * pairs_poploc / nrow(sons)))
cat("Current sample (CLP u MLP v1.0, post-all-filters):", OLD_N, "\n")

# NOTE ON THE COMPARISON BAR. pairs_poploc is pre-filter; 12,246 is
# post-filter (multiple-father drop, MAD <= 4, school == 1). To compare
# like with like, re-run 01_cleaning-script.R and use its "After
# multiple-father drop" cat() line, not the final n.
#
# Also: MLP v2 alone is one link source competing against two. The union
# is >= either alone by construction, so v2 starts at a disadvantage.
# Coming within ~10% of the old count is evidence v2 is the better source,
# not evidence it is worse.

# Do the same fathers reappear across years? HIK makes this checkable,
# which the old wide file could not do (HISTID is year-specific).
multi_dad = father_records |>
  group_by(son_hik) |>
  summarise(n_dad_hik = n_distinct(dad_hik[!is.na(dad_hik) & dad_hik != ""]),
            n_years   = n_distinct(year), .groups = "drop")

cat("\nSons by distinct father HIKs across years:\n")
print(multi_dad |> count(n_dad_hik, name = "sons"))
cat("(n_dad_hik > 1 = stepfather, remarriage, or a bad link. n_dad_hik == 0 =\n")
cat(" father observed but never himself linked - keep him, single-year occ.)\n")

################################################################################
################### Q2. RACIAL RECLASSIFICATION ################################
################################################################################

cat("\n================ Q2: RECLASSIFICATION ================\n")

if (nonaian_pre1940 == 0) {
  cat("Not answerable from this extract - see Q0.\n")
} else {
  son_race = son_records |>
    filter(year < 1940) |>
    group_by(son_hik) |>
    summarise(ever_nonaian = any(race != AIAN),
              n_years = n(), .groups = "drop")

  cat("Sons with >= 1 pre-1940 record:", nrow(son_race), "\n")
  cat("Ever coded non-AIAN before 1940:", sum(son_race$ever_nonaian),
      sprintf("(%.1f%%)\n", 100 * mean(son_race$ever_nonaian)))

  cat("\nFathers ever coded non-AIAN:\n")
  print(father_records |> group_by(son_hik) |>
          summarise(ever = any(dad_race != AIAN), .groups = "drop") |>
          summarise(n = sum(ever), pct = round(100 * mean(ever), 1)))

  cat("\nThis is the number that justifies a rebuild. If it is large, the\n")
  cat("current CLP race-matched sample excludes reclassifiers by construction\n")
  cat("and the selection caveat in the methods section is understated.\n")
}

################################################################################
######################## Q3. POPLOC VALIDATION #################################
################################################################################

cat("\n================ Q3: POPLOC QUALITY ================\n")

cat("POPLOC == 0 among son-year records (no co-resident father):\n")
print(son_records |> group_by(year) |>
        summarise(poploc_zero = sum(poploc == 0),
                  pct = round(100 * mean(poploc == 0), 1)))

implausible = father_records |>
  mutate(gap = dad_age - age) |>
  summarise(n           = n(),
            gap_lt_15   = sum(gap < 15, na.rm = TRUE),
            gap_gt_60   = sum(gap > 60, na.rm = TRUE),
            gap_median  = median(gap, na.rm = TRUE))
print(implausible)

cat("\nIf gap_lt_15 + gap_gt_60 is under ~1%, POPLOC is trustworthy and the\n")

cat("multiple-father paste-and-drop logic in 01_cleaning-script.R can go.\n")
cat("Add a SEX check on the father record if SEX was included in the extract.\n")

################################################################################
############################## VERDICT #########################################
################################################################################

cat("\n================ VERDICT ================\n")
cat("Q0 extract valid (non-AIAN pre-1940 present): ", nonaian_pre1940 > 0, "\n")
cat("Q1 pairs recovered (pre-filter):              ", pairs_poploc, "\n")
cat("   vs current post-filter sample:             ", OLD_N, "\n")
cat("Q3 median father-son age gap:                 ", implausible$gap_median, "\n")
cat("\nRebuild only if Q0 is TRUE and either Q1 clears the adjusted bar or Q2\n")
cat("is large. Otherwise: document the version gap in one footnote and spend\n")
cat("the remaining time on the section 10 backlog.\n")

