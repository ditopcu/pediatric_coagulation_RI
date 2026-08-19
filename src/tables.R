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

make_descriptive_table <- function(test_name, data = all_data_3sd) {
  data |>
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
    all_data |> group_by(test) |>
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
    pct_female = round(n_female / n_raw * 100, 1)
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
  pct_female = round(sum(pub_1_by_test$n_female) / sum(pub_1_by_test$n_raw) * 100, 1),
  age_mean = round(mean(all_data$age_year), 1),
  age_sd = round(sd(all_data$age_year), 1),
  age_median = round(median(all_data$age_year), 1),
  age_min = round(min(all_data$age_year), 1),
  age_max = round(max(all_data$age_year), 1)
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

# Continuous GAMLSS — one row per reported partition, summarising the per-age
# curves in PUB_table_S3. Each row is the envelope of the ages it covers: the
# lowest lower limit and the highest upper limit, so the interval contains
# every age in the partition. The same statistic is used for all four rows.
# The aPTT partitions follow the refineR ones, 1 to <12 and 12 to <18, so
# age 12 belongs to the second row.
cont_envelope <- function(tname, group, ages) {
  cont_ri_all |>
    filter(test == tname, age %in% ages) |>
    summarise(test = tname, age_group = group,
              lower = round(min(lower), 2), upper = round(max(upper), 2),
              lower_range = sprintf("%.2f-%.2f", min(lower), max(lower)),
              upper_range = sprintf("%.2f-%.2f", min(upper), max(upper)))
}

cont_partitions <- list(
  list("PT",         "1-18 y",  1:18),
  list("aPTT",       "1-12 y",  1:11),
  list("aPTT",       "12-18 y", 12:18),
  list("Fibrinogen", "1-18 y",  1:18)
)

cont_summary <- bind_rows(
  lapply(cont_partitions, function(p) cont_envelope(p[[1]], p[[2]], p[[3]]))
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

copy_if_exists("analysis6_fp_bootstrap_ci_RAW.csv", "PUB_table_3_false_positive_rates")
copy_if_exists("analysis2_guven_transferability_RAW.csv", "PUB_table_4_guven_transferability")
copy_if_exists("analysis5_international_comparison.csv", "PUB_table_5_international_comparison")
copy_if_exists("analysis5b_partition_harris_boyd.csv", "PUB_table_6_harris_boyd")


# ==============================================================================
# PUB TABLE S1 — Combined descriptive statistics
# ==============================================================================

message("[TABLE] PUB_table_S1_descriptive_all ...")

# Supplementary Table S1 describes the full analysis set. The TEZ tables above
# and Figure 1 use the 3SD display subset, which is a plotting aid and never a
# reported statistic, so it must not feed the supplement.
pub_s1 <- bind_rows(
  make_descriptive_table("PT",         all_data) |> mutate(test = "PT"),
  make_descriptive_table("aPTT",       all_data) |> mutate(test = "aPTT"),
  make_descriptive_table("Fibrinogen", all_data) |> mutate(test = "Fibrinogen")
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
# Consistency guards. Every check runs, then the script stops if any failed, so
# a single run reports everything that has drifted rather than the first item.
# docs/PROVENANCE.md records which output file feeds which manuscript table.

.checks <- list()
check_that <- function(label, ok, detail = "") {
  .checks[[length(.checks) + 1]] <<- list(label = label, ok = isTRUE(ok), detail = detail)
  message(sprintf("  [%s] %s%s", if (isTRUE(ok)) "OK" else "FAIL", label,
                  if (!isTRUE(ok) && nzchar(detail)) paste0(" -- ", detail) else ""))
}

message("\n[VERIFY] Reference limits ...")

# TEZ_tablo_4_9 is ri_comparison written out unchanged.
check_that("TEZ_tablo_4_9 equals ri_comparison",
           isTRUE(all.equal(as.data.frame(ri_comparison),
                            as.data.frame(read_csv(file.path(out_dir,
                              "TEZ_tablo_4_9_indirect_vs_direct.csv"),
                              show_col_types = FALSE)))))

# PUB_table_2 carries ri_comparison's refineR and Direct rows unchanged.
p2_cmp <- pub_table_2 |>
  filter(source %in% c("indirect_refineR", "direct")) |>
  select(test, age_group, method, n, lower, lower_ci_lo, lower_ci_hi,
         upper, upper_ci_lo, upper_ci_hi) |>
  arrange(test, age_group, method)
ri_cmp <- ri_comparison |>
  select(test, age_group, method, n, lower, lower_ci_lo, lower_ci_hi,
         upper, upper_ci_lo, upper_ci_hi) |>
  arrange(test, age_group, method)
d2 <- anti_join(ri_cmp, p2_cmp, by = names(ri_cmp))
check_that("PUB_table_2 matches ri_comparison", nrow(d2) == 0,
           sprintf("%d row(s) differ", nrow(d2)))
if (nrow(d2) > 0) print(d2)

# Manufacturer rows in PUB_table_2 come from ri_manufacturer.
p2_mfr <- pub_table_2 |> filter(source == "manufacturer") |>
  select(test, lower, upper) |> arrange(test)
check_that("PUB_table_2 manufacturer rows match ri_manufacturer",
           isTRUE(all.equal(as.data.frame(p2_mfr),
                            as.data.frame(arrange(select(ri_manufacturer, test, lower, upper),
                                                  test)))))

# common.R is the only place the reference limits are defined. A second copy in
# another script drifts silently the first time one side is edited.
src_files <- setdiff(list.files("src", pattern = "\\.R$", full.names = TRUE),
                     "src/common.R")
dupes <- src_files[vapply(src_files, function(f)
  any(grepl("^\\s*ri_(refiner|manufacturer|direct)\\s*<-", readLines(f, warn = FALSE))),
  logical(1))]
check_that("reference limits defined only in common.R", length(dupes) == 0,
           paste(basename(dupes), collapse = ", "))

message("\n[VERIFY] Sample sizes ...")

expected_n <- c(PT = 17000L, aPTT = 17001L, Fibrinogen = 3166L)
actual_n <- vapply(names(expected_n),
                   function(t) sum(all_data$test == t), integer(1))
for (t in names(expected_n)) {
  check_that(sprintf("%s n = %d", t, expected_n[[t]]),
             actual_n[[t]] == expected_n[[t]],
             sprintf("found %d", actual_n[[t]]))
}
check_that("total across analytes = 37167", nrow(all_data) == 37167L,
           sprintf("found %d", nrow(all_data)))

message("\n[VERIFY] Continuous reference intervals ...")

check_that("PUB_table_S3 covers ages 1-18 for all three analytes",
           nrow(pub_s3) == 54 && setequal(unique(pub_s3$test),
                                          c("PT", "aPTT", "Fibrinogen")),
           sprintf("%d rows", nrow(pub_s3)))

# Every continuous row in the written file must be the envelope of the ages its
# own partition covers, on the same definition for all four. This fails if one
# row is switched to a different summary statistic or a partition boundary
# moves away from the refineR one.
p2_cont <- read_csv(file.path(out_dir, "PUB_table_2_reference_intervals.csv"),
                    show_col_types = FALSE) |>
  filter(source == "continuous_GAMLSS") |>
  select(test, age_group, lower, upper) |> arrange(test, age_group)
expected_cont <- bind_rows(
  lapply(cont_partitions, function(p) cont_envelope(p[[1]], p[[2]], p[[3]]))) |>
  select(test, age_group, lower, upper) |> arrange(test, age_group)
check_that("PUB_table_2 continuous rows are the envelope of their own age range",
           isTRUE(all.equal(as.data.frame(p2_cont), as.data.frame(expected_cont))))

message("\n[VERIFY] Flagging tables ...")

# PUB_table_3 and PUB_table_4 carry the manuscript's Table 3 and Table S4, so
# they copy the _RAW outputs: those describe the full analysis set, whereas the
# unsuffixed outputs describe the 3SD display subset.
for (p in list(c("PUB_table_3_false_positive_rates.csv", "analysis6_fp_bootstrap_ci_RAW.csv"),
               c("PUB_table_4_guven_transferability.csv", "analysis2_guven_transferability_RAW.csv"))) {
  a <- read_csv(file.path(out_dir, p[1]), show_col_types = FALSE)
  b <- read_csv(file.path("data/processed", p[2]), show_col_types = FALSE)
  check_that(sprintf("%s is a faithful copy of %s", p[1], p[2]),
             isTRUE(all.equal(as.data.frame(a), as.data.frame(b))))
}

message("\n[VERIFY] Table inventory:")
all_tables <- list.files(out_dir, pattern = "\\.csv$")
message(sprintf("  Total tables: %d", length(all_tables)))
for (f in all_tables) {
  df <- read_csv(file.path(out_dir, f), show_col_types = FALSE)
  na_count <- sum(is.na(df))
  message(sprintf("  %s: %d rows x %d cols, %d NAs", f, nrow(df), ncol(df), na_count))
}

failed <- Filter(function(x) !x$ok, .checks)
message(sprintf("\n[VERIFY] %d checks run, %d failed.", length(.checks), length(failed)))
if (length(failed) > 0) {
  stop(sprintf("Consistency checks failed:\n  - %s",
               paste(vapply(failed, function(x)
                 paste0(x$label, if (nzchar(x$detail)) paste0(" (", x$detail, ")") else ""),
                 character(1)), collapse = "\n  - ")),
       call. = FALSE)
}

message("\n[DONE] All tables generated.")
