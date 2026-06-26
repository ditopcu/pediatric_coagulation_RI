# ==============================================================================
# MANUSCRIPT ANALYSES
# Pediatric Coagulation Reference Intervals
# ==============================================================================
# Analysis 1: Misclassification / false-positive rates (Table 3)
# Analysis 2: Guven et al. RI transferability (CLSI EP28-A3c)
# Analysis 3: Age-specific false-positive curves (Figure)
# Analysis 5: International comparison table
# Analysis 5b: Age partitioning justification (Harris-Boyd)
# Analysis 6: Bootstrap CI for FP rates
# ==============================================================================
# Requires: source("src/common.R") to be run first (via run_all.R or manually).
# ==============================================================================

if (!exists("PAL")) source("src/common.R")

# Ensure age_int is integer (pub.R may have converted to ordered factor)
all_data_3sd <- all_data_3sd |> mutate(age_int = as.integer(as.character(age_int)))


# ==============================================================================
# ANALYSIS 1: MISCLASSIFICATION / FALSE-POSITIVE RATES (Table 3)
# ==============================================================================
# For each test and age: what % of 3SD-filtered data falls outside each RI set?
# "Outside" = below lower OR above upper limit

message("[ANALYSIS 1] Computing false-positive rates by RI source...")

# Build a lookup of all RI limits per test/age
# 1) Manufacturer: same limits for all ages
ri_mfr_lookup <- ri_manufacturer |>
  crossing(age_int = 1:18) |>
  rename(mfr_lower = lower, mfr_upper = upper)

# 2) refineR: aPTT has age-group split at 12
ri_ref_lookup <- bind_rows(
  ri_refiner |> filter(test == "PT") |> crossing(age_int = 1:18),
  ri_refiner |> filter(test == "Fibrinogen") |> crossing(age_int = 1:18),
  ri_refiner |> filter(test == "aPTT", age_group == "1-12") |> crossing(age_int = 1:11),
  ri_refiner |> filter(test == "aPTT", age_group == "12-18") |> crossing(age_int = 12:18)
) |>
  rename(ref_lower = lower, ref_upper = upper) |>
  select(test, age_int, ref_lower, ref_upper)

# 3) Continuous: already per-age
ri_cont_lookup <- cont_ri_all |>
  rename(age_int = age, cont_lower = lower, cont_upper = upper) |>
  select(test, age_int, cont_lower, cont_upper)

# 4) Guven: per-age for PT and aPTT only
ri_guven_lookup <- guven_ri |>
  rename(age_int = age, guven_lower = lower, guven_upper = upper) |>
  select(test, age_int, guven_lower, guven_upper)

# Join all limits to individual observations, compute outside flags
fp_detail <- all_data_3sd |>
  left_join(ri_mfr_lookup,   by = c("test", "age_int")) |>
  left_join(ri_ref_lookup,   by = c("test", "age_int")) |>
  left_join(ri_cont_lookup,  by = c("test", "age_int")) |>
  left_join(ri_guven_lookup, by = c("test", "age_int")) |>
  mutate(
    out_mfr  = (result_num < mfr_lower  | result_num > mfr_upper),
    out_ref  = (result_num < ref_lower  | result_num > ref_upper),
    out_cont = if_else(!is.na(cont_lower),
                       result_num < cont_lower | result_num > cont_upper,
                       NA),
    out_guven = if_else(!is.na(guven_lower),
                        result_num < guven_lower | result_num > guven_upper,
                        NA)
  )

# Aggregate per test/age
fp_by_age <- fp_detail |>
  group_by(test, age_int) |>
  summarise(
    n = n(),
    fp_manufacturer = mean(out_mfr, na.rm = TRUE) * 100,
    fp_refiner      = mean(out_ref, na.rm = TRUE) * 100,
    fp_continuous   = mean(out_cont, na.rm = TRUE) * 100,
    fp_guven        = mean(out_guven, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Summary table (overall per test)
fp_summary <- fp_by_age |>
  group_by(test) |>
  summarise(
    N = sum(n),
    FP_Manufacturer = weighted.mean(fp_manufacturer, n, na.rm = TRUE),
    FP_refineR      = weighted.mean(fp_refiner,      n, na.rm = TRUE),
    FP_Continuous   = weighted.mean(fp_continuous,    n, na.rm = TRUE),
    FP_Guven        = weighted.mean(fp_guven,         n, na.rm = TRUE),
    .groups = "drop"
  )

message("\n=== TABLE 3: Overall False-Positive Rates (%) ===")
print(fp_summary, n = Inf)

# Save CSV
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write_csv(fp_by_age, "data/processed/analysis1_fp_by_age.csv")
write_csv(fp_summary, "data/processed/analysis1_fp_summary.csv")
message("[OK] Analysis 1 saved to data/processed/\n")


# ==============================================================================
# ANALYSIS 2: GUVEN et al. TRANSFERABILITY (CLSI EP28-A3c)
# ==============================================================================
# CLSI EP28-A3c validation: draw 20 random samples from our data per age.
# If <= 2/20 (10%) fall outside Guven RI => transferable for that age.
# Repeat B=200 bootstrap iterations to get proportion of "pass".

message("[ANALYSIS 2] Guven et al. RI transferability (CLSI EP28-A3c, 20-sample)...")

set.seed(42)
B <- 200
n_validate <- 20
threshold <- 2  # max allowed outside

validate_one_age <- function(test_val, age_val, lo, hi) {
  d <- all_data_3sd |> filter(test == test_val, age_int == age_val)
  n_total <- nrow(d)
  if (n_total < n_validate) {
    return(tibble(n_total = n_total, pass_rate = NA_real_,
                  n_outside = NA_integer_, prop_outside = NA_real_,
                  binom_p = NA_real_, binom_reject = NA_character_))
  }
  passes <- replicate(B, {
    samp <- slice_sample(d, n = n_validate, replace = TRUE)
    n_outside <- sum(samp$result_num < lo | samp$result_num > hi)
    n_outside <= threshold
  })
  # Exact binomial test on full age-group data: H0: p <= 0.05 vs H1: p > 0.05
  # Rejects when observed out-of-RI proportion significantly exceeds nominal 5%.
  n_out_full <- sum(d$result_num < lo | d$result_num > hi)
  bt <- binom.test(n_out_full, n_total, p = 0.05, alternative = "greater")
  tibble(n_total = n_total, pass_rate = mean(passes) * 100,
         n_outside = n_out_full, prop_outside = n_out_full / n_total * 100,
         binom_p = bt$p.value,
         binom_reject = ifelse(bt$p.value < 0.05, "Yes", "No"))
}

transfer_results <- guven_ri |>
  mutate(res = pmap(list(test, age, lower, upper), validate_one_age)) |>
  unnest(res) |>
  mutate(transferable = ifelse(is.na(pass_rate), NA_character_,
                               ifelse(pass_rate >= 90, "Yes", "No")))

message("\n=== TABLE: Guven RI Transferability ===")
transfer_results |>
  select(test, age, n = n_total, guven_lower = lower, guven_upper = upper,
         pass_rate, transferable) |>
  print(n = Inf)

# Summary
transfer_summary <- transfer_results |>
  group_by(test) |>
  summarise(
    ages_tested = n(),
    ages_pass   = sum(transferable == "Yes", na.rm = TRUE),
    ages_fail   = sum(transferable == "No",  na.rm = TRUE),
    ages_na     = sum(is.na(transferable)),
    mean_pass_rate = mean(pass_rate, na.rm = TRUE),
    ages_binom_reject = sum(binom_reject == "Yes", na.rm = TRUE),
    .groups = "drop"
  )
message("\n=== Transferability Summary ===")
print(transfer_summary)

write_csv(transfer_results, "data/processed/analysis2_guven_transferability.csv")
message("[OK] Analysis 2 saved.\n")


# ==============================================================================
# ANALYSIS 3: AGE-SPECIFIC FALSE-POSITIVE CURVES (Figure)
# ==============================================================================

message("[ANALYSIS 3] Generating age-specific FP rate curves...")

fp_long <- fp_by_age |>
  pivot_longer(cols = starts_with("fp_"),
               names_to = "ri_source", values_to = "fp_rate",
               names_prefix = "fp_") |>
  mutate(
    ri_source = case_when(
      ri_source == "manufacturer" ~ "Manufacturer",
      ri_source == "refiner"      ~ "refineR",
      ri_source == "continuous"   ~ "Continuous (GAMLSS)",
      ri_source == "guven"        ~ "Guven et al.",
      TRUE ~ ri_source
    ),
    ri_source = factor(ri_source,
                       levels = c("Manufacturer", "Guven et al.",
                                  "refineR", "Continuous (GAMLSS)"))
  )

make_fp_panel <- function(test_name, tag_label) {
  fp_long |>
    filter(test == test_name, !is.na(fp_rate)) |>
    ggplot(aes(x = age_int, y = fp_rate, colour = ri_source, linetype = ri_source)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 5, linetype = "dotted", colour = PAL$base2) +
    annotate("text", x = 17, y = 5.8, label = "Expected 5%",
             colour = PAL$base2, size = 2.8, hjust = 1) +
    scale_colour_manual(values = c(
      "Manufacturer"       = PAL$base2,
      "Guven et al."       = PAL$accent3,
      "refineR"            = PAL$refiner,
      "Continuous (GAMLSS)" = PAL$contin
    )) +
    scale_linetype_manual(values = c(
      "Manufacturer"       = "dashed",
      "Guven et al."       = "dotdash",
      "refineR"            = "solid",
      "Continuous (GAMLSS)" = "solid"
    )) +
    scale_x_continuous(breaks = seq(1, 18, 1)) +
    labs(x = "Age, years",
         y = "False-positive rate, %",
         tag = tag_label) +
    theme_tufte_academic() +
    theme(
      plot.tag = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
}

p_fp_pt   <- make_fp_panel("PT", "A")
p_fp_aptt <- make_fp_panel("aPTT", "B")
p_fp_fib  <- make_fp_panel("Fibrinogen", "C") +
  theme(legend.position = "bottom")

fig_fp_curves <- p_fp_pt / p_fp_aptt / p_fp_fib +
  plot_annotation(
    theme = theme_tufte_academic()
  )

save_fig(fig_fp_curves, "figure_fp_curves_by_age", width = 180, height = 280)

write_csv(fp_long, "data/processed/analysis3_fp_curves_data.csv")
message("[OK] Analysis 3 saved.\n")


# ==============================================================================
# ANALYSIS 5: INTERNATIONAL COMPARISON TABLE
# ==============================================================================
# Columns: Study, Year, Country, Analyzer, Method, Population, N, Age_group,
#          Test, Lower, Upper, Percentiles, Unit

message("[ANALYSIS 5] Building international comparison table...")

intl_comparison <- tribble(
  ~study, ~year, ~country, ~analyzer, ~method, ~n, ~age_group, ~test, ~lower, ~upper, ~percentiles, ~unit,

  # --- Our study ---
  "This study (refineR)", 2026, "Turkey", "Cobas t511/t711", "Indirect (refineR)", 17000, "1-18 y", "PT", 8.56, 10.90, "2.5-97.5", "s",
  "This study (refineR)", 2026, "Turkey", "Cobas t511/t711", "Indirect (refineR)", 10091, "1-12 y", "aPTT", 24.30, 35.30, "2.5-97.5", "s",
  "This study (refineR)", 2026, "Turkey", "Cobas t511/t711", "Indirect (refineR)", 6910, "12-18 y", "aPTT", 23.30, 33.90, "2.5-97.5", "s",
  "This study (refineR)", 2026, "Turkey", "Cobas t511/t711", "Indirect (refineR)", 3166, "1-18 y", "Fibrinogen", 1.89, 3.79, "2.5-97.5", "g/L",
  "This study (Direct)", 2026, "Turkey", "Cobas t511/t711", "Direct (a posteriori)", 305, "1-18 y", "PT", 8.69, 10.80, "2.5-97.5", "s",
  "This study (Direct)", 2026, "Turkey", "Cobas t511/t711", "Direct (a posteriori)", 170, "1-12 y", "aPTT", 22.96, 34.74, "2.5-97.5", "s",
  "This study (Direct)", 2026, "Turkey", "Cobas t511/t711", "Direct (a posteriori)", 150, "12-18 y", "aPTT", 24.10, 32.69, "2.5-97.5", "s",
  "This study (Direct)", 2026, "Turkey", "Cobas t511/t711", "Direct (a posteriori)", 198, "1-18 y", "Fibrinogen", 2.10, 3.56, "2.5-97.5", "g/L",

  # --- Guven et al. 2026 (age-pooled approximation from Table 1/2, ages 1-17) ---
  "Guven et al.", 2026, "Turkey", "Cobas t511", "Indirect (refineR)", 1056, "1-18 y", "PT", 7.8, 10.7, "2.5-97.5", "s",
  "Guven et al.", 2026, "Turkey", "Cobas t511", "Indirect (refineR)", 1044, "1-18 y", "aPTT", 24.3, 36.2, "2.5-97.5", "s",

  # --- Luo et al. 2026 (Sysmex CN-6000, Chinese children, Table 1) ---
  "Luo et al.", 2026, "China", "Sysmex CN-6000", "Direct", 500, "6m-3 y", "PT", 10.1, 12.4, "2.5-97.5", "s",
  "Luo et al.", 2026, "China", "Sysmex CN-6000", "Direct", 264, "3-12 y", "PT", 10.3, 12.4, "2.5-97.5", "s",
  "Luo et al.", 2026, "China", "Sysmex CN-6000", "Direct", 270, "12-18 y", "PT", 10.4, 13.0, "2.5-97.5", "s",
  "Luo et al.", 2026, "China", "Sysmex CN-6000", "Direct", 500, "6m-3 y", "aPTT", 23.4, 32.6, "2.5-97.5", "s",
  "Luo et al.", 2026, "China", "Sysmex CN-6000", "Direct", 534, "3-18 y", "aPTT", 24.7, 32.9, "2.5-97.5", "s",
  "Luo et al.", 2026, "China", "Sysmex CN-6000", "Direct", 258, "1-3 y", "Fibrinogen", 1.52, 3.58, "2.5-97.5", "g/L",
  "Luo et al.", 2026, "China", "Sysmex CN-6000", "Direct", 534, "3-18 y", "Fibrinogen", 1.74, 3.65, "2.5-97.5", "g/L",

  # --- Weidhofer et al. 2018 (STA-Compact/BCS-XP, Austrian, indirect continuous) ---
  # PT reported in %, not seconds — include for completeness with note
  "Weidhofer et al.", 2018, "Austria", "STA-Compact / BCS-XP", "Indirect (continuous)", 35492, "1-18 y", "PT", NA, NA, "2.5-97.5", "% (not comparable)",
  "Weidhofer et al.", 2018, "Austria", "STA-Compact / BCS-XP", "Indirect (continuous)", 55100, "1-18 y", "aPTT", 26.0, 42.0, "2.5-97.5", "s",
  "Weidhofer et al.", 2018, "Austria", "STA-Compact / BCS-XP", "Indirect (continuous)", 49789, "1-18 y", "Fibrinogen", 150, 450, "2.5-97.5", "mg/dL",

  # --- Manufacturer (Roche Cobas t511/t711 package insert) ---
  "Manufacturer (Roche)", NA, NA, "Cobas t511/t711", "Package insert", NA, "All ages", "PT", 8.4, 10.6, NA, "s",
  "Manufacturer (Roche)", NA, NA, "Cobas t511/t711", "Package insert", NA, "All ages", "aPTT", 23.6, 30.6, NA, "s",
  "Manufacturer (Roche)", NA, NA, "Cobas t511/t711", "Package insert", NA, "All ages", "Fibrinogen", 1.93, 4.12, NA, "g/L"
)

write_csv(intl_comparison, "data/processed/analysis5_international_comparison.csv")
message("[OK] Analysis 5: International comparison table saved.\n")


# ==============================================================================
# ANALYSIS 5b: AGE PARTITIONING JUSTIFICATION (aPTT at age 12)
# ==============================================================================
# Harris-Boyd method: compare two adjacent groups.
# z = |mean1 - mean2| / sqrt(sd1^2/n1 + sd2^2/n2)
# If z > z* (critical value based on desired significance), partition is justified.
# Standard threshold: z* = 3 (conservative) or z* from normal distribution.
# Also report: CI overlap for 2.5th and 97.5th percentiles from refineR by age/sex.

message("[ANALYSIS 5b] aPTT age partitioning justification...")

# Load refineR by-age-sex results
aptt_ref_sex <- tryCatch(readRDS("data/processed/aptt_ref_by_age_sex_raw.RDS"), error = function(e) NULL)

# refineR partition data is optional (requires refineR package loaded)
# Harris-Boyd test below uses raw data directly

# Harris-Boyd z-test: compare ages <12 vs >=12
aptt_young <- all_data_3sd |> filter(test == "aPTT", age_int < 12)
aptt_old   <- all_data_3sd |> filter(test == "aPTT", age_int >= 12)

hb_stats <- tibble(
  group = c("<12 y", ">=12 y"),
  n     = c(nrow(aptt_young), nrow(aptt_old)),
  mean  = c(mean(aptt_young$result_num), mean(aptt_old$result_num)),
  sd    = c(sd(aptt_young$result_num), sd(aptt_old$result_num))
)

z_harris_boyd <- abs(hb_stats$mean[1] - hb_stats$mean[2]) /
  sqrt(hb_stats$sd[1]^2 / hb_stats$n[1] + hb_stats$sd[2]^2 / hb_stats$n[2])

# Critical z for partitioning (Harris & Boyd 1990: z* = 3 * sqrt(n_total) / (n1 + n2))
# Simplified: if z > 3, partitioning is statistically justified
partition_justified <- z_harris_boyd > 3

message(sprintf("\n=== Harris-Boyd Partitioning Test (aPTT at age 12) ==="))
message(sprintf("  <12 y: n=%d, mean=%.2f, sd=%.2f", hb_stats$n[1], hb_stats$mean[1], hb_stats$sd[1]))
message(sprintf("  >=12 y: n=%d, mean=%.2f, sd=%.2f", hb_stats$n[2], hb_stats$mean[2], hb_stats$sd[2]))
message(sprintf("  z = %.2f (threshold = 3)", z_harris_boyd))
message(sprintf("  Partition justified: %s", ifelse(partition_justified, "YES", "NO")))

# Also test PT (should NOT require partitioning)
pt_young <- all_data_3sd |> filter(test == "PT", age_int < 12)
pt_old   <- all_data_3sd |> filter(test == "PT", age_int >= 12)

z_pt <- abs(mean(pt_young$result_num) - mean(pt_old$result_num)) /
  sqrt(sd(pt_young$result_num)^2 / nrow(pt_young) + sd(pt_old$result_num)^2 / nrow(pt_old))

message(sprintf("\n=== Harris-Boyd Test (PT at age 12) ==="))
message(sprintf("  z = %.2f — Partition justified: %s", z_pt, ifelse(z_pt > 3, "YES", "NO")))

# Save
partition_results <- tibble(
  test = c("aPTT", "PT"),
  cutpoint = c(12, 12),
  n_young = c(hb_stats$n[1], nrow(pt_young)),
  n_old   = c(hb_stats$n[2], nrow(pt_old)),
  mean_young = c(hb_stats$mean[1], mean(pt_young$result_num)),
  mean_old   = c(hb_stats$mean[2], mean(pt_old$result_num)),
  sd_young = c(hb_stats$sd[1], sd(pt_young$result_num)),
  sd_old   = c(hb_stats$sd[2], sd(pt_old$result_num)),
  z_score  = c(z_harris_boyd, z_pt),
  justified = c(partition_justified, z_pt > 3)
)

write_csv(partition_results, "data/processed/analysis5b_partition_harris_boyd.csv")
message("[OK] Analysis 5b saved.\n")


# ==============================================================================
# ANALYSIS 6: BOOTSTRAP CI FOR FP RATES
# ==============================================================================
# For each test: bootstrap (B=1000) the overall FP rate per RI source.
# Report 95% CI.

message("[ANALYSIS 6] Bootstrap CI for false-positive rates...")

set.seed(123)
B_fp <- 1000

boot_fp_ci <- function(data, test_name, lo, hi, B = 1000) {
  d <- data |> filter(test == test_name)
  vals <- d$result_num
  n <- length(vals)
  boot_rates <- replicate(B, {
    s <- sample(vals, n, replace = TRUE)
    mean(s < lo | s > hi) * 100
  })
  tibble(
    fp_mean = mean(boot_rates),
    fp_ci_lo = quantile(boot_rates, 0.025),
    fp_ci_hi = quantile(boot_rates, 0.975)
  )
}

# Manufacturer
fp_boot_mfr <- ri_manufacturer |>
  rowwise() |>
  mutate(res = list(boot_fp_ci(all_data_3sd, test, lower, upper, B_fp))) |>
  unnest(res) |>
  mutate(ri_source = "Manufacturer") |>
  select(test, ri_source, fp_mean, fp_ci_lo, fp_ci_hi)

# refineR (pooled per test — use the single-partition RI for PT/Fib, weighted for aPTT)
fp_boot_ref <- bind_rows(
  ri_refiner |> filter(test == "PT") |> slice(1) |>
    rowwise() |> mutate(res = list(boot_fp_ci(all_data_3sd, test, lower, upper, B_fp))) |> unnest(res),
  # aPTT: compute on full aPTT data using age-appropriate RI
  tibble(test = "aPTT") |>
    mutate(res = list({
      d <- all_data_3sd |> filter(test == "aPTT")
      ri_young <- ri_refiner |> filter(test == "aPTT", age_group == "1-12")
      ri_old   <- ri_refiner |> filter(test == "aPTT", age_group == "12-18")
      vals_y <- d |> filter(age_int < 12) |> pull(result_num)
      vals_o <- d |> filter(age_int >= 12) |> pull(result_num)
      boot_rates <- replicate(B_fp, {
        sy <- sample(vals_y, length(vals_y), replace = TRUE)
        so <- sample(vals_o, length(vals_o), replace = TRUE)
        n_out <- sum(sy < ri_young$lower | sy > ri_young$upper) +
                 sum(so < ri_old$lower | so > ri_old$upper)
        n_out / (length(sy) + length(so)) * 100
      })
      tibble(fp_mean = mean(boot_rates),
             fp_ci_lo = quantile(boot_rates, 0.025),
             fp_ci_hi = quantile(boot_rates, 0.975))
    })) |> unnest(res),
  ri_refiner |> filter(test == "Fibrinogen") |> slice(1) |>
    rowwise() |> mutate(res = list(boot_fp_ci(all_data_3sd, test, lower, upper, B_fp))) |> unnest(res)
) |>
  mutate(ri_source = "refineR") |>
  select(test, ri_source, fp_mean, fp_ci_lo, fp_ci_hi)

# Continuous (use per-age RI, compute on joined data)
boot_cont_fn <- function(tname) {
  d <- all_data_3sd |> filter(test == tname) |>
    left_join(cont_ri_all |> rename(age_int = age), by = c("test", "age_int")) |>
    filter(!is.na(lower))
  vals <- d$result_num; lo <- d$lower; hi <- d$upper
  n <- nrow(d)
  boot_rates <- replicate(B_fp, {
    idx <- sample(n, n, replace = TRUE)
    mean(vals[idx] < lo[idx] | vals[idx] > hi[idx]) * 100
  })
  tibble(test = tname, ri_source = "Continuous",
         fp_mean = mean(boot_rates),
         fp_ci_lo = quantile(boot_rates, 0.025),
         fp_ci_hi = quantile(boot_rates, 0.975))
}
fp_boot_cont <- bind_rows(
  boot_cont_fn("PT"), boot_cont_fn("aPTT"), boot_cont_fn("Fibrinogen")
)

# Guven
boot_guven_fn <- function(tname) {
  d <- all_data_3sd |> filter(test == tname) |>
    left_join(guven_ri |> rename(age_int = age), by = c("test", "age_int")) |>
    filter(!is.na(lower))
  vals <- d$result_num; lo <- d$lower; hi <- d$upper
  n <- nrow(d)
  boot_rates <- replicate(B_fp, {
    idx <- sample(n, n, replace = TRUE)
    mean(vals[idx] < lo[idx] | vals[idx] > hi[idx]) * 100
  })
  tibble(test = tname, ri_source = "Guven",
         fp_mean = mean(boot_rates),
         fp_ci_lo = quantile(boot_rates, 0.025),
         fp_ci_hi = quantile(boot_rates, 0.975))
}
fp_boot_guven <- bind_rows(boot_guven_fn("PT"), boot_guven_fn("aPTT"))

fp_boot_all <- bind_rows(fp_boot_mfr, fp_boot_ref, fp_boot_cont, fp_boot_guven) |>
  mutate(fp_label = sprintf("%.1f (%.1f-%.1f)", fp_mean, fp_ci_lo, fp_ci_hi))

message("\n=== TABLE: FP Rates with 95% Bootstrap CI ===")
fp_boot_all |>
  select(test, ri_source, fp_label) |>
  pivot_wider(names_from = ri_source, values_from = fp_label) |>
  print()

write_csv(fp_boot_all, "data/processed/analysis6_fp_bootstrap_ci.csv")
message("[OK] Analysis 6 saved.\n")


# ==============================================================================
# DONE
# ==============================================================================
message("============================================================")
message("[DONE] All analyses complete.")
message("  CSV outputs: data/processed/analysis*")
message("  Figure: figures/*/figure_fp_curves_by_age.*")
message("============================================================")
