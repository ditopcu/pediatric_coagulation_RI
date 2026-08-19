# ==============================================================================
# CROSS-TEST EXCLUSION SENSITIVITY (reviewer R2-4) -- exploratory only
# ==============================================================================
# Reviewer question: were PT/aPTT results excluded when another coagulation
# test on the SAME sample was pathological (e.g. low fibrinogen)?
# They were not (each analyte was modeled independently). This script
# quantifies what would change if such exclusions had been applied.
#
# Feasibility (verified): the three final data sets share the accession number
# (sample_id): PT & aPTT co-measured on 14,934 samples; PT & Fib on 2,467;
# aPTT & Fib on 2,246; all three on 2,049.
#
# Primary variant (A): a sample is "co-abnormal" for analyte T when ANY
# co-measured OTHER analyte falls outside its own study refineR RI
# (age-appropriate partition for aPTT). Those samples are excluded from T's
# data set and T's RI is re-derived (refineR point estimate, NBootstrap = 1,
# matching analysis4_sensitivity.R) and flagging at the PUBLISHED RI is
# re-computed on the restricted set.
#
# Secondary variant (B): clinical-threshold exclusion -- co-measured
# fibrinogen < 1.5 g/L excludes the sample from PT and aPTT.
#
# NOTE: computed on the final analysis data sets (first result per patient),
# i.e. the same population as the manuscript, so deltas are directly
# interpretable against the published values. Main RIs/tables are NOT changed.
#
# Outputs (new files only):
#   data/processed/REVIEW_A5_crosstest_RI_REV.csv
#   data/processed/REVIEW_A5_crosstest_FP_REV.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(refineR)
})

if (!exists("PAL")) source("src/common.R")

dat <- all_data |>
  mutate(age_int = as.integer(as.character(age_int)))

# --- Own-RI limits per test/age (published refineR RIs) ----------------------
own_ri_lookup <- bind_rows(
  ri_refiner |> filter(test == "PT") |> crossing(age_int = 1:18),
  ri_refiner |> filter(test == "Fibrinogen") |> crossing(age_int = 1:18),
  ri_refiner |> filter(test == "aPTT", age_group == "1-12") |> crossing(age_int = 1:11),
  ri_refiner |> filter(test == "aPTT", age_group == "12-18") |> crossing(age_int = 12:18)
) |>
  select(test, age_int, own_lower = lower, own_upper = upper)

dat <- dat |>
  left_join(own_ri_lookup, by = c("test", "age_int")) |>
  mutate(own_abnormal = result_num < own_lower | result_num > own_upper)

# --- Per-sample abnormality of each analyte ----------------------------------
sample_flags <- dat |>
  select(sample_id, test, own_abnormal, result_num) |>
  pivot_wider(names_from = test,
              values_from = c(own_abnormal, result_num),
              values_fn = first)

# Variant A: co-abnormal = any co-measured OTHER analyte outside its own RI
co_abnormal_for <- function(target) {
  others <- setdiff(c("PT", "aPTT", "Fibrinogen"), target)
  flags <- sample_flags |>
    mutate(co_ab = rowSums(across(all_of(paste0("own_abnormal_", others)),
                                  ~coalesce(.x, FALSE))) > 0) |>
    select(sample_id, co_ab)
  flags
}

# Variant B: co-measured fibrinogen < 1.5 g/L (clinical threshold)
fib_low_flags <- sample_flags |>
  mutate(fib_low = coalesce(result_num_Fibrinogen < 1.5, FALSE)) |>
  select(sample_id, fib_low)

# --- refineR helper (point estimate) -----------------------------------------
refiner_point <- function(values, test_name, variant, age_grp) {
  n <- length(values)
  fit <- try(refineR::findRI(values, NBootstrap = 1), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(tibble(test = test_name, variant = variant, age_group = age_grp,
                  n = n, lower = NA_real_, median = NA_real_, upper = NA_real_))
  }
  ri <- as.data.frame(getRI(fit, RIperc = c(0.025, 0.50, 0.975)))
  tibble(test = test_name, variant = variant, age_group = age_grp, n = n,
         lower  = round(ri$PointEst[ri$Percentile == 0.025], 2),
         median = round(ri$PointEst[ri$Percentile == 0.500], 2),
         upper  = round(ri$PointEst[ri$Percentile == 0.975], 2))
}

run_variant <- function(dat_v, variant_label) {
  out <- list()
  d <- dat_v |> filter(test == "PT")
  out[[1]] <- refiner_point(d$result_num, "PT", variant_label, "1-18 y")
  d <- dat_v |> filter(test == "aPTT", age_int < 12)
  out[[2]] <- refiner_point(d$result_num, "aPTT", variant_label, "1-12 y")
  d <- dat_v |> filter(test == "aPTT", age_int >= 12)
  out[[3]] <- refiner_point(d$result_num, "aPTT", variant_label, "12-18 y")
  d <- dat_v |> filter(test == "Fibrinogen")
  out[[4]] <- refiner_point(d$result_num, "Fibrinogen", variant_label, "1-18 y")
  bind_rows(out)
}

flagging_at_published <- function(dat_v, variant_label) {
  dat_v |>
    mutate(grp = case_when(
      test == "PT"                    ~ "PT 1-18 y",
      test == "aPTT" & age_int < 12   ~ "aPTT 1-12 y",
      test == "aPTT" & age_int >= 12  ~ "aPTT 12-18 y",
      test == "Fibrinogen"            ~ "Fibrinogen 1-18 y"
    )) |>
    group_by(grp) |>
    summarise(variant = variant_label,
              n = n(),
              flag_pct = round(mean(own_abnormal) * 100, 2),
              .groups = "drop")
}

set.seed(42)
message("[CROSSTEST] Baseline (no cross-test exclusion)...")
ri_all <- list(run_variant(dat, "V0_baseline"))
fp_all <- list(flagging_at_published(dat, "V0_baseline"))
excl_summary <- list()

message("[CROSSTEST] Variant A: exclude samples with co-abnormal other analyte...")
for (t in c("PT", "aPTT", "Fibrinogen")) {
  flags <- co_abnormal_for(t)
  d_t <- dat |> filter(test == t) |> left_join(flags, by = "sample_id")
  n_excl <- sum(d_t$co_ab, na.rm = TRUE)
  excl_summary[[paste0("A_", t)]] <- tibble(
    variant = "VA_co_abnormal", test = t,
    n_before = nrow(d_t), n_excluded = n_excl,
    pct_excluded = round(n_excl / nrow(d_t) * 100, 2))
  d_keep <- d_t |> filter(!co_ab)
  ri_all[[length(ri_all) + 1]] <- run_variant(d_keep, "VA_co_abnormal") |>
    filter(test == t)
  fp_all[[length(fp_all) + 1]] <- flagging_at_published(d_keep, "VA_co_abnormal") |>
    filter(grepl(t, grp))
}

message("[CROSSTEST] Variant B: exclude PT/aPTT samples with co-measured fibrinogen < 1.5 g/L...")
for (t in c("PT", "aPTT")) {
  d_t <- dat |> filter(test == t) |> left_join(fib_low_flags, by = "sample_id")
  n_excl <- sum(d_t$fib_low, na.rm = TRUE)
  excl_summary[[paste0("B_", t)]] <- tibble(
    variant = "VB_fib_below_1.5", test = t,
    n_before = nrow(d_t), n_excluded = n_excl,
    pct_excluded = round(n_excl / nrow(d_t) * 100, 2))
  d_keep <- d_t |> filter(!fib_low)
  ri_all[[length(ri_all) + 1]] <- run_variant(d_keep, "VB_fib_below_1.5") |>
    filter(test == t)
  fp_all[[length(fp_all) + 1]] <- flagging_at_published(d_keep, "VB_fib_below_1.5") |>
    filter(grepl(t, grp))
}

ri_out <- bind_rows(ri_all) |> arrange(test, age_group, variant)
fp_out <- bind_rows(fp_all) |> arrange(grp, variant)
excl_out <- bind_rows(excl_summary)

cat("\n=== EXCLUSION COUNTS ===\n");  print(excl_out, n = Inf)
cat("\n=== RI SENSITIVITY ===\n");    print(ri_out, n = Inf)
cat("\n=== FLAGGING AT PUBLISHED RI ===\n"); print(fp_out, n = Inf)

write_csv(bind_rows(ri_out |> mutate(section = "RI"),
                    excl_out |> mutate(section = "exclusions")),
          "data/processed/REVIEW_A5_crosstest_RI_REV.csv")
write_csv(fp_out, "data/processed/REVIEW_A5_crosstest_FP_REV.csv")
message("[OK] Cross-test sensitivity outputs saved.")
