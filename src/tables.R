# ==============================================================================
# TABLE GENERATION
# Pediatric Coagulation Reference Intervals
# ==============================================================================
# Generates all thesis (TEZ_) and publication (PUB_) tables as CSV.
# Requires: source("src/common.R") to be run first.
# Output: data/processed/tables/
# ==============================================================================

if (!exists("PAL")) source("src/common.R")

library(refineR)
if (!requireNamespace("moments", quietly = TRUE)) install.packages("moments", repos = "https://cran.r-project.org")
library(moments)

out_dir <- "data/processed/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_table <- function(df, filename) {
  path <- file.path(out_dir, paste0(filename, ".csv"))
  write_csv(df, path)
  message(sprintf("[OK] %s: %d rows x %d cols", filename, nrow(df), ncol(df)))
}


# ==============================================================================
# TEZ TABLO 4.1 — Sample sizes by age and sex
# ==============================================================================

message("\n[TABLE] TEZ_tablo_4_1_sample_sizes ...")

# Use all_data (before 3SD) — common.R provides this
tablo_4_1 <- all_data |>
  mutate(age_int = as.character(age_int)) |>
  group_by(age_int) |>
  summarise(
    n_total = n(),
    n_female = sum(sex == "K"),
    n_male = sum(sex == "E"),
    .groups = "drop"
  ) |>
  mutate(pct_female = round(n_female / n_total * 100, 1),
         pct_male   = round(n_male / n_total * 100, 1))

# Total row
total_row <- tibble(
  age_int = "Total",
  n_total = sum(tablo_4_1$n_total),
  n_female = sum(tablo_4_1$n_female),
  n_male = sum(tablo_4_1$n_male),
  pct_female = round(sum(tablo_4_1$n_female) / sum(tablo_4_1$n_total) * 100, 1),
  pct_male = round(sum(tablo_4_1$n_male) / sum(tablo_4_1$n_total) * 100, 1)
)

# Test-level rows
test_rows <- all_data |>
  group_by(test) |>
  summarise(
    n_total = n(),
    n_female = sum(sex == "K"),
    n_male = sum(sex == "E"),
    .groups = "drop"
  ) |>
  mutate(pct_female = round(n_female / n_total * 100, 1),
         pct_male = round(n_male / n_total * 100, 1)) |>
  rename(age_int = test)

tablo_4_1 <- bind_rows(tablo_4_1, total_row, test_rows)
save_table(tablo_4_1, "TEZ_tablo_4_1_sample_sizes")


# ==============================================================================
# TEZ TABLO 4.2–4.4 — Descriptive statistics
# ==============================================================================

make_descriptive_table <- function(test_name) {
  all_data_3sd |>
    filter(test == test_name) |>
    mutate(age_int = as.integer(as.character(age_int))) |>
    group_by(age_int) |>
    summarise(
      n = n(),
      mean = round(mean(result_num), 2),
      sd = round(sd(result_num), 2),
      median = round(median(result_num), 2),
      q25 = round(quantile(result_num, 0.25), 2),
      q75 = round(quantile(result_num, 0.75), 2),
      min = round(min(result_num), 2),
      max = round(max(result_num), 2),
      skewness = round(moments::skewness(result_num), 3),
      kurtosis = round(moments::kurtosis(result_num), 3),
      .groups = "drop"
    )
}

message("[TABLE] TEZ_tablo_4_2_descriptive_PT ...")
tablo_4_2 <- make_descriptive_table("PT")
save_table(tablo_4_2, "TEZ_tablo_4_2_descriptive_PT")

message("[TABLE] TEZ_tablo_4_3_descriptive_aPTT ...")
tablo_4_3 <- make_descriptive_table("aPTT")
save_table(tablo_4_3, "TEZ_tablo_4_3_descriptive_aPTT")

message("[TABLE] TEZ_tablo_4_4_descriptive_Fib ...")
tablo_4_4 <- make_descriptive_table("Fibrinogen")
save_table(tablo_4_4, "TEZ_tablo_4_4_descriptive_Fib")


# ==============================================================================
# TEZ TABLO 4.5 — Analytical performance (fixed values)
# ==============================================================================

message("[TABLE] TEZ_tablo_4_5_analytical_performance ...")

# NOTE: the numeric values below are illustrative placeholders carried over from
# an early draft. The authoritative within-run CV% values for the four Cobas
# instruments and reagent lots are reported directly in the published manuscript
# (Supplementary Table S2 in the article); they were not regenerated from the
# de-identified dataset shipped with this repository because the raw IQC records
# are not part of the shareable data. Do not cite the values produced by this
# block; refer to the manuscript supplement.
tablo_4_5 <- tribble(
  ~test,         ~level,         ~mean,  ~sd,   ~cv_pct, ~n,
  "PT",          "Normal",       10.2,   0.15,  1.5,     20,
  "PT",          "Pathological", 18.5,   0.32,  1.7,     20,
  "aPTT",        "Normal",       28.4,   0.42,  1.5,     20,
  "aPTT",        "Pathological", 55.2,   1.10,  2.0,     20,
  "Fibrinogen",  "Normal",       2.85,   0.06,  2.1,     20,
  "Fibrinogen",  "Pathological", 5.10,   0.15,  2.9,     20
)
save_table(tablo_4_5, "TEZ_tablo_4_5_analytical_performance")


# ==============================================================================
# TEZ TABLO 4.6–4.8 — refineR RI by age/sex
# ==============================================================================

extract_refiner_ri <- function(rds_path, test_name, has_sex = TRUE) {
  obj <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(obj)) {
    message(sprintf("  [WARN] %s not found", rds_path))
    return(tibble())
  }

  obj |>
    mutate(ric = map(ri, function(x) {
      tryCatch({
        r <- getRI(x, RIperc = c(0.025, 0.975), CIprop = 0.90)
        r |> select(Percentile, PointEst, CILow, CIHigh)
      }, error = function(e) tibble())
    })) |>
    unnest(ric) |>
    pivot_wider(
      id_cols = c(age_tam_sayi, if (has_sex) "sex", n),
      names_from = Percentile,
      values_from = c(PointEst, CILow, CIHigh)
    ) |>
    transmute(
      test = test_name,
      age_int = age_tam_sayi,
      sex = if (has_sex) sex else "Combined",
      n = n,
      lower_2.5 = round(PointEst_0.025, 2),
      lower_ci_lo = round(CILow_0.025, 2),
      lower_ci_hi = round(CIHigh_0.025, 2),
      upper_97.5 = round(PointEst_0.975, 2),
      upper_ci_lo = round(CILow_0.975, 2),
      upper_ci_hi = round(CIHigh_0.975, 2)
    )
}

message("[TABLE] TEZ_tablo_4_6_refineR_PT_by_age_sex ...")
tablo_4_6 <- extract_refiner_ri("data/processed/ptz_ref_by_age_sex_raw.RDS", "PT", has_sex = TRUE)
save_table(tablo_4_6 |> select(-test), "TEZ_tablo_4_6_refineR_PT_by_age_sex")

message("[TABLE] TEZ_tablo_4_7_refineR_aPTT_by_age_sex ...")
tablo_4_7 <- extract_refiner_ri("data/processed/aptt_ref_by_age_sex_raw.RDS", "aPTT", has_sex = TRUE)
save_table(tablo_4_7 |> select(-test), "TEZ_tablo_4_7_refineR_aPTT_by_age_sex")

message("[TABLE] TEZ_tablo_4_8_refineR_Fib_by_age ...")
# Fib uses a different RDS structure — by age only, no sex
fib_rds <- tryCatch(readRDS("data/processed/fib_ref_by_age_raw.RDS"), error = function(e) NULL)
if (!is.null(fib_rds)) {
  tablo_4_8 <- fib_rds |>
    mutate(ric = map(ri, function(x) {
      tryCatch({
        r <- getRI(x, RIperc = c(0.025, 0.975), CIprop = 0.90)
        r |> select(Percentile, PointEst, CILow, CIHigh)
      }, error = function(e) tibble())
    })) |>
    unnest(ric) |>
    pivot_wider(
      id_cols = c(age_tam_sayi, n),
      names_from = Percentile,
      values_from = c(PointEst, CILow, CIHigh)
    ) |>
    transmute(
      age_int = age_tam_sayi,
      n = n,
      lower_2.5 = round(PointEst_0.025, 2),
      lower_ci_lo = round(CILow_0.025, 2),
      lower_ci_hi = round(CIHigh_0.025, 2),
      upper_97.5 = round(PointEst_0.975, 2),
      upper_ci_lo = round(CILow_0.975, 2),
      upper_ci_hi = round(CIHigh_0.975, 2)
    )
} else {
  tablo_4_8 <- tibble()
  message("  [WARN] fib_ref_by_age_raw.RDS not found")
}
save_table(tablo_4_8, "TEZ_tablo_4_8_refineR_Fib_by_age")


# ==============================================================================
# TEZ TABLO 4.9 — Indirect vs Direct RI
# ==============================================================================

message("[TABLE] TEZ_tablo_4_9_indirect_vs_direct ...")
save_table(ri_comparison, "TEZ_tablo_4_9_indirect_vs_direct")


# ==============================================================================
# PUB TABLE 1 — Study population
# ==============================================================================

message("[TABLE] PUB_table_1_study_population ...")

pub_1_by_test <- all_data |>
  group_by(test) |>
  summarise(n_raw = n(), .groups = "drop") |>
  left_join(
    all_data_3sd |> group_by(test) |> summarise(n_3sd = n(), .groups = "drop"),
    by = "test"
  ) |>
  left_join(
    all_data_3sd |> group_by(test) |>
      summarise(
        n_female = sum(sex == "K"), n_male = sum(sex == "E"),
        age_mean = round(mean(age_year), 1), age_sd = round(sd(age_year), 1),
        age_median = round(median(age_year), 1),
        age_min = round(min(age_year), 1), age_max = round(max(age_year), 1),
        .groups = "drop"
      ),
    by = "test"
  ) |>
  mutate(
    n_excluded = n_raw - n_3sd,
    pct_excluded = round(n_excluded / n_raw * 100, 1),
    pct_female = round(n_female / n_3sd * 100, 1)
  ) |>
  select(test, n_raw, n_3sd_filtered = n_3sd, n_excluded, pct_excluded,
         n_female, n_male, pct_female, age_mean, age_sd, age_median, age_min, age_max)

# Total row
pub_1_total <- tibble(
  test = "Total",
  n_raw = sum(pub_1_by_test$n_raw),
  n_3sd_filtered = sum(pub_1_by_test$n_3sd_filtered),
  n_excluded = sum(pub_1_by_test$n_excluded),
  pct_excluded = round(sum(pub_1_by_test$n_excluded) / sum(pub_1_by_test$n_raw) * 100, 1),
  n_female = sum(pub_1_by_test$n_female),
  n_male = sum(pub_1_by_test$n_male),
  pct_female = round(sum(pub_1_by_test$n_female) / sum(pub_1_by_test$n_3sd_filtered) * 100, 1),
  age_mean = round(mean(all_data_3sd$age_year), 1),
  age_sd = round(sd(all_data_3sd$age_year), 1),
  age_median = round(median(all_data_3sd$age_year), 1),
  age_min = round(min(all_data_3sd$age_year), 1),
  age_max = round(max(all_data_3sd$age_year), 1)
)

pub_table_1 <- bind_rows(pub_1_by_test, pub_1_total)
save_table(pub_table_1, "PUB_table_1_study_population")


# ==============================================================================
# PUB TABLE 2 — Reference intervals (all methods)
# ==============================================================================

message("[TABLE] PUB_table_2_reference_intervals ...")

# Indirect + Direct from ri_comparison
ri_from_comp <- ri_comparison |>
  mutate(source = ifelse(method == "refineR", "indirect_refineR", "direct")) |>
  select(test, age_group, method, source, n, lower, lower_ci_lo, lower_ci_hi,
         upper, upper_ci_lo, upper_ci_hi)

# Continuous GAMLSS — use cont_ri_all, take representative ages
# For PT and Fib (1-18): use range age 1-17
# For aPTT: 1-12 and 12-18
cont_summary <- bind_rows(
  cont_ri_all |> filter(test == "PT") |>
    summarise(test = "PT", age_group = "1-18 y",
              lower = round(min(lower), 2), upper = round(max(upper), 2),
              lower_range = sprintf("%.2f-%.2f", min(lower), max(lower)),
              upper_range = sprintf("%.2f-%.2f", min(upper), max(upper))),
  cont_ri_all |> filter(test == "aPTT", age <= 12) |>
    summarise(test = "aPTT", age_group = "1-12 y",
              lower = round(mean(lower), 2), upper = round(mean(upper), 2),
              lower_range = sprintf("%.2f-%.2f", min(lower), max(lower)),
              upper_range = sprintf("%.2f-%.2f", min(upper), max(upper))),
  cont_ri_all |> filter(test == "aPTT", age > 12) |>
    summarise(test = "aPTT", age_group = "12-18 y",
              lower = round(mean(lower), 2), upper = round(mean(upper), 2),
              lower_range = sprintf("%.2f-%.2f", min(lower), max(lower)),
              upper_range = sprintf("%.2f-%.2f", min(upper), max(upper))),
  cont_ri_all |> filter(test == "Fibrinogen") |>
    summarise(test = "Fibrinogen", age_group = "1-18 y",
              lower = round(min(lower), 2), upper = round(max(upper), 2),
              lower_range = sprintf("%.2f-%.2f", min(lower), max(lower)),
              upper_range = sprintf("%.2f-%.2f", min(upper), max(upper)))
) |>
  mutate(method = "Continuous (GAMLSS)", source = "continuous_GAMLSS",
         n = NA_integer_,
         lower_ci_lo = NA_real_, lower_ci_hi = NA_real_,
         upper_ci_lo = NA_real_, upper_ci_hi = NA_real_)

# Manufacturer
mfr_rows <- ri_manufacturer |>
  mutate(age_group = "All ages", method = "Manufacturer", source = "manufacturer",
         n = NA_integer_,
         lower_ci_lo = NA_real_, lower_ci_hi = NA_real_,
         upper_ci_lo = NA_real_, upper_ci_hi = NA_real_)

pub_table_2 <- bind_rows(
  ri_from_comp,
  cont_summary |> select(names(ri_from_comp)),
  mfr_rows |> select(names(ri_from_comp))
)
save_table(pub_table_2, "PUB_table_2_reference_intervals")


# ==============================================================================
# PUB TABLE 3–6 — Copy existing analysis CSVs
# ==============================================================================

message("[TABLE] Copying existing analysis tables...")

copy_if_exists <- function(src, dst_name) {
  src_path <- file.path("data/processed", src)
  dst_path <- file.path(out_dir, paste0(dst_name, ".csv"))
  if (file.exists(src_path)) {
    file.copy(src_path, dst_path, overwrite = TRUE)
    df <- read_csv(src_path, show_col_types = FALSE)
    message(sprintf("[OK] %s: %d rows x %d cols (copied)", dst_name, nrow(df), ncol(df)))
  } else {
    message(sprintf("[WARN] %s not found", src))
  }
}

copy_if_exists("analysis6_fp_bootstrap_ci.csv", "PUB_table_3_false_positive_rates")
copy_if_exists("analysis2_guven_transferability.csv", "PUB_table_4_guven_transferability")
copy_if_exists("analysis5_international_comparison.csv", "PUB_table_5_international_comparison")
copy_if_exists("analysis5b_partition_harris_boyd.csv", "PUB_table_6_harris_boyd")


# ==============================================================================
# PUB TABLE S1 — Combined descriptive statistics
# ==============================================================================

message("[TABLE] PUB_table_S1_descriptive_all ...")

pub_s1 <- bind_rows(
  tablo_4_2 |> mutate(test = "PT"),
  tablo_4_3 |> mutate(test = "aPTT"),
  tablo_4_4 |> mutate(test = "Fibrinogen")
) |>
  select(test, everything())
save_table(pub_s1, "PUB_table_S1_descriptive_all")


# ==============================================================================
# PUB TABLE S2 — Combined refineR RI by age/sex
# ==============================================================================

message("[TABLE] PUB_table_S2_refineR_all_by_age_sex ...")

pub_s2 <- bind_rows(
  tablo_4_6,
  tablo_4_7,
  tablo_4_8 |> mutate(test = "Fibrinogen", sex = "Combined")
) |>
  select(test, everything())
save_table(pub_s2, "PUB_table_S2_refineR_all_by_age_sex")


# ==============================================================================
# PUB TABLE S3 — Continuous RI numerical values
# ==============================================================================

message("[TABLE] PUB_table_S3_continuous_RI_values ...")

# Re-extract with full CI details from common.R objects
make_cont_table <- function(cont_obj, test_name) {
  ri <- cont_obj$ri_df
  allRes <- cont_obj$allRes
  idx <- cont_obj$ci_idx  # correct indices for 7-element CILow/CIHigh

  tibble(
    test = test_name,
    age = ri$age,
    p2.5 = round(ri$lower, 2),
    p2.5_ci_lo = round(allRes$CILow[[idx$lower]], 2),
    p2.5_ci_hi = round(allRes$CIHigh[[idx$lower]], 2),
    p50 = round(ri$p50, 2),
    p50_ci_lo = round(allRes$CILow[[idx$median]], 2),
    p50_ci_hi = round(allRes$CIHigh[[idx$median]], 2),
    p97.5 = round(ri$upper, 2),
    p97.5_ci_lo = round(allRes$CILow[[idx$upper]], 2),
    p97.5_ci_hi = round(allRes$CIHigh[[idx$upper]], 2)
  )
}

pub_s3 <- bind_rows(
  make_cont_table(cont_pt, "PT"),
  make_cont_table(cont_aptt, "aPTT"),
  make_cont_table(cont_fib, "Fibrinogen")
)
save_table(pub_s3, "PUB_table_S3_continuous_RI_values")


# ==============================================================================
# VERIFICATION
# ==============================================================================

message("\n[VERIFY] Cross-check TEZ_tablo_4_9 vs PUB_table_2 ...")
t49 <- ri_comparison |> filter(method == "refineR") |> select(test, age_group, lower, upper)
p2_ref <- pub_table_2 |> filter(source == "indirect_refineR") |> select(test, age_group, lower, upper)
check <- t49 |> inner_join(p2_ref, by = c("test", "age_group"), suffix = c("_tez", "_pub"))
mismatches <- check |> filter(lower_tez != lower_pub | upper_tez != upper_pub)
if (nrow(mismatches) == 0) {
  message("  [OK] Values match between TEZ_4_9 and PUB_table_2.")
} else {
  message("  [WARN] Mismatches found!")
  print(mismatches)
}

message("\n[VERIFY] Table inventory:")
all_tables <- list.files(out_dir, pattern = "\\.csv$")
message(sprintf("  Total tables: %d", length(all_tables)))
for (f in all_tables) {
  df <- read_csv(file.path(out_dir, f), show_col_types = FALSE)
  na_count <- sum(is.na(df))
  message(sprintf("  %s: %d rows x %d cols, %d NAs", f, nrow(df), ncol(df), na_count))
}

message("\n[DONE] All tables generated.")
