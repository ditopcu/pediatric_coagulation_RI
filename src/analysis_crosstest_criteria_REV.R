# ==============================================================================
# CROSS-TEST EXCLUSION, EXTERNAL CRITERION (reviewer R2-4)
# ==============================================================================
# Reviewer question: were results excluded when another coagulation test on the
# same sample was pathological -- for example, were PT and aPTT dropped when the
# co-measured fibrinogen was low? They were not; each analyte was modeled
# independently. This script quantifies what would change if they had been.
#
# It supersedes the exploratory src/analysis_crosstest_sensitivity_REV.R in one
# respect: "pathological" is defined by an EXTERNAL reference range rather than
# by this study's own refineR limits. Using the study's own limits to decide
# which samples enter the estimation of those same limits is circular.
#
# Criterion 1 -- MANUFACTURER range. This is the range the laboratory reported
# against during the study period, so a result outside it is what the clinician
# would have seen flagged at the time. Single range per analyte; the
# manufacturer does not publish an age-partitioned aPTT range, so it is applied
# irrespective of age.
#
# Criterion 2 -- DIRECT (a posteriori) limits, CLSI EP28-A3c, from the separately
# recruited cohorts. Reported beside the manufacturer criterion so the result can
# be shown not to depend on which external range defines "pathological". These
# limits come from cohorts that are not part of the data being filtered, so they
# are free of the circularity as well.
#
# Exclusion rule (both criteria, all combinations): a sample is dropped from
# analyte T when ANY co-measured OTHER analyte on that sample falls outside the
# criterion range, in either direction. The target analyte's own value is never
# used to exclude it.
#
# Computed on the final analysis data sets (first result per patient), the same
# population as the manuscript, so the values sit directly beside the published
# ones. refineR point estimates only (NBootstrap = 1), no confidence intervals.
# Main RIs and tables are NOT changed; this is supplementary material.
#
# Outputs (new files only):
#   data/processed/REVIEW_A5_crosstest_criteria_RI_REV.csv          (both criteria)
#   data/processed/REVIEW_A5_crosstest_criteria_exclusions_REV.csv  (both criteria)
#   data/processed/REVIEW_A5_crosstest_criteria_flagging_REV.csv    (both criteria)
#   data/processed/tables/SUPPLEMENT_crosstest_REV.xlsx             (both criteria)
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(janitor); library(readr); library(tidyr)
  library(tibble); library(openxlsx); library(refineR)
})

# --- Reference limits --------------------------------------------------------
# Reference limits come from src/common.R, which is the single source of truth
# for `ri_refiner`, `ri_manufacturer` and `ri_direct`. Sourcing it also loads
# the analysis data and the GAMLSS continuous curves from cont_out/, so that
# directory has to hold a complete run before this script will start.
source("src/common.R")

# Expand an age-partitioned limit table into one row per integer age. The
# partition labels follow common.R: "1-12" is age_year < 12, so age_int 1-11.
expand_by_age <- function(ri, ages = 1:18) {
  ri |>
    rowwise() |>
    mutate(age_int = list(switch(age_group,
                                 "1-18"  = ages,
                                 "1-12"  = ages[ages <  12],
                                 "12-18" = ages[ages >= 12]))) |>
    ungroup() |>
    unnest(age_int) |>
    select(test, age_int, lower, upper)
}

# Manufacturer: one range per analyte, no age partition.
crit_manufacturer <- ri_manufacturer |>
  crossing(age_int = 1:18) |>
  select(test, age_int, lower, upper)

# Direct: aPTT partitioned at 12 years, as elsewhere in the manuscript.
crit_direct <- expand_by_age(ri_direct)

criteria <- list(Manufacturer = crit_manufacturer, Direct = crit_direct)

# --- Data --------------------------------------------------------------------
main_files <- c(
  PT         = "data/coa_results/tce 2026 1-18 pt.xlsx",
  aPTT       = "data/coa_results/tce 2026 1-18 aptt.xlsx",
  Fibrinogen = "data/coa_results/tce 2026 1-18 fib.xlsx"
)

load_main <- function(path, tname) {
  suppressMessages(read_excel(path)) |>
    clean_names() |>
    mutate(test = tname, age_int = as.integer(floor(age_year))) |>
    select(sample_id, test, age_int, result_num)
}

dat <- bind_rows(
  load_main(main_files["PT"],         "PT"),
  load_main(main_files["aPTT"],       "aPTT"),
  load_main(main_files["Fibrinogen"], "Fibrinogen")
)

cat("\n=== INPUT ===\n")
print(dat |> count(test))

# --- refineR helper ----------------------------------------------------------
refiner_point <- function(values, test_name, variant, age_grp) {
  n <- length(values)
  if (n < 100) {
    return(tibble(test = test_name, variant = variant, age_group = age_grp,
                  n = n, lower = NA_real_, median = NA_real_, upper = NA_real_,
                  note = "n<100 skipped"))
  }
  fit <- try(refineR::findRI(values, NBootstrap = 1), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(tibble(test = test_name, variant = variant, age_group = age_grp,
                  n = n, lower = NA_real_, median = NA_real_, upper = NA_real_,
                  note = "refineR failed"))
  }
  ri <- as.data.frame(getRI(fit, RIperc = c(0.025, 0.50, 0.975)))
  tibble(test = test_name, variant = variant, age_group = age_grp, n = n,
         lower  = round(ri$PointEst[ri$Percentile == 0.025], 2),
         median = round(ri$PointEst[ri$Percentile == 0.500], 2),
         upper  = round(ri$PointEst[ri$Percentile == 0.975], 2),
         note   = "")
}

# RI for one analyte from a given row set
ri_for <- function(d, target, variant) {
  if (target == "aPTT") {
    bind_rows(
      refiner_point(d$result_num[d$age_int <  12], "aPTT", variant, "1-12 y"),
      refiner_point(d$result_num[d$age_int >= 12], "aPTT", variant, "12-18 y")
    )
  } else {
    refiner_point(d$result_num, target, variant, "1-18 y")
  }
}

# Flagging at the published refineR RI, for a given row set
flag_at_published <- function(d, target, variant) {
  lut <- bind_rows(
    ri_refiner |> filter(test == "PT")                         |> crossing(age_int = 1:18),
    ri_refiner |> filter(test == "Fibrinogen")                 |> crossing(age_int = 1:18),
    ri_refiner |> filter(test == "aPTT", age_group == "1-12")  |> crossing(age_int = 1:11),
    ri_refiner |> filter(test == "aPTT", age_group == "12-18") |> crossing(age_int = 12:18)
  ) |> select(test, age_int, pub_lower = lower, pub_upper = upper)

  d |>
    left_join(lut, by = c("test", "age_int")) |>
    mutate(grp = case_when(
      test == "PT"                   ~ "PT 1-18 y",
      test == "aPTT" & age_int <  12 ~ "aPTT 1-12 y",
      test == "aPTT" & age_int >= 12 ~ "aPTT 12-18 y",
      test == "Fibrinogen"           ~ "Fibrinogen 1-18 y"
    ),
    outside = result_num < pub_lower | result_num > pub_upper) |>
    group_by(grp) |>
    summarise(variant = variant, n = n(),
              flag_pct = round(mean(outside) * 100, 2), .groups = "drop")
}

# --- Baseline ----------------------------------------------------------------
set.seed(42)
message("[CROSSTEST] Baseline (no cross-test exclusion)...")
ri_all   <- list()
fp_all   <- list()
excl_all <- list()

for (t in c("PT", "aPTT", "Fibrinogen")) {
  d_t <- dat |> filter(test == t)
  ri_all[[length(ri_all) + 1]] <- ri_for(d_t, t, "V0_baseline")
  fp_all[[length(fp_all) + 1]] <- flag_at_published(d_t, t, "V0_baseline")
}

# --- One pass per criterion --------------------------------------------------
for (cname in names(criteria)) {
  crit <- criteria[[cname]]
  message("[CROSSTEST] Criterion: ", cname, " ...")

  # Abnormality of every measurement under this criterion
  ab <- dat |>
    left_join(crit, by = c("test", "age_int")) |>
    mutate(abn = result_num < lower | result_num > upper) |>
    select(sample_id, test, abn)

  for (t in c("PT", "aPTT", "Fibrinogen")) {
    others <- setdiff(c("PT", "aPTT", "Fibrinogen"), t)

    # Per-other-analyte flags on the target's samples
    d_t <- dat |> filter(test == t)
    flags <- d_t |> select(sample_id)
    for (o in others) {
      f_o <- ab |> filter(test == o) |> select(sample_id, !!o := abn)
      flags <- flags |> left_join(f_o, by = "sample_id")
    }
    flags <- flags |>
      mutate(across(all_of(others), ~ coalesce(.x, FALSE))) |>
      mutate(drop_any = if_any(all_of(others), identity))

    n_by_other <- vapply(others, function(o) sum(flags[[o]]), integer(1))
    excl_all[[paste(cname, t)]] <- tibble(
      criterion = cname, test = t,
      n_before = nrow(d_t),
      n_excluded = sum(flags$drop_any),
      pct_excluded = round(sum(flags$drop_any) / nrow(d_t) * 100, 2),
      n_after = nrow(d_t) - sum(flags$drop_any),
      due_to_1 = paste0(others[1], ": ", n_by_other[1]),
      due_to_2 = paste0(others[2], ": ", n_by_other[2])
    )

    d_keep <- d_t |> filter(!flags$drop_any)
    vlabel <- paste0("V_", cname)
    ri_all[[length(ri_all) + 1]] <- ri_for(d_keep, t, vlabel)
    fp_all[[length(fp_all) + 1]] <- flag_at_published(d_keep, t, vlabel)
  }
}

ri_out   <- bind_rows(ri_all)   |> arrange(test, age_group, variant)
fp_out   <- bind_rows(fp_all)   |> arrange(grp, variant)
excl_out <- bind_rows(excl_all) |> arrange(criterion, test)

cat("\n=== EXCLUSION COUNTS ===\n");            print(excl_out, n = Inf)
cat("\n=== RI BY VARIANT ===\n");               print(ri_out,   n = Inf)
cat("\n=== FLAGGING AT PUBLISHED RI ===\n");    print(fp_out,   n = Inf)

# --- Save --------------------------------------------------------------------
write_csv(ri_out,   "data/processed/REVIEW_A5_crosstest_criteria_RI_REV.csv")
write_csv(excl_out, "data/processed/REVIEW_A5_crosstest_criteria_exclusions_REV.csv")
write_csv(fp_out,   "data/processed/REVIEW_A5_crosstest_criteria_flagging_REV.csv")

# Supplement carries both criteria side by side.
out_xlsx <- "data/processed/tables/SUPPLEMENT_crosstest_REV.xlsx"
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                  halign = "center", border = "bottom")
add_sheet <- function(nm, df) {
  addWorksheet(wb, nm); writeData(wb, nm, df, headerStyle = hs)
  freezePane(wb, nm, firstRow = TRUE)
  setColWidths(wb, nm, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet("Exclusions", excl_out)
add_sheet("RI",         ri_out)
add_sheet("Flagging",   fp_out)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("\n[OK] Saved: %s\n", out_xlsx))
