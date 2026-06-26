# ==============================================================================
# RAW-DATA (NO 3SD FILTER) COMPARISON
# Repeats Analyses 1, 2, 3, 6 using all_data (raw, no global 3SD filter)
# instead of all_data_3sd. Outputs are written with the _RAW suffix only;
# NO existing files are overwritten.
# ==============================================================================
# Requires: source("src/common.R") to be run first.
# Output files:
#   data/processed/analysis1_fp_by_age_RAW.csv
#   data/processed/analysis1_fp_summary_RAW.csv
#   data/processed/analysis2_guven_transferability_RAW.csv
#   data/processed/analysis3_fp_curves_data_RAW.csv
#   data/processed/analysis6_fp_bootstrap_ci_RAW.csv
#   data/processed/COMPARISON_3sd_vs_raw.csv
#   data/processed/tables/COMPARISON_3sd_vs_raw.xlsx
#   figures/{TIFF_600DPI,PNG_300DPI}/figure_fp_curves_by_age_RAW.{tiff,png}
# ==============================================================================

if (!exists("PAL")) source("src/common.R")

suppressPackageStartupMessages({
  library(openxlsx)
})

# Use RAW data (no 3SD filter) -------------------------------------------------
all_data_raw <- all_data |> mutate(age_int = as.integer(as.character(age_int)))

message(sprintf("[RAW] Input sizes -- raw: %d, 3SD-filtered: %d (delta=%d)",
                nrow(all_data_raw), nrow(all_data_3sd),
                nrow(all_data_raw) - nrow(all_data_3sd)))


# ==============================================================================
# ANALYSIS 1 (RAW): MISCLASSIFICATION / FP RATES
# ==============================================================================

message("[ANALYSIS 1 RAW] FP rates by RI source (no 3SD)...")

ri_mfr_lookup <- ri_manufacturer |>
  crossing(age_int = 1:18) |>
  rename(mfr_lower = lower, mfr_upper = upper)

ri_ref_lookup <- bind_rows(
  ri_refiner |> filter(test == "PT") |> crossing(age_int = 1:18),
  ri_refiner |> filter(test == "Fibrinogen") |> crossing(age_int = 1:18),
  ri_refiner |> filter(test == "aPTT", age_group == "1-12") |> crossing(age_int = 1:11),
  ri_refiner |> filter(test == "aPTT", age_group == "12-18") |> crossing(age_int = 12:18)
) |>
  rename(ref_lower = lower, ref_upper = upper) |>
  select(test, age_int, ref_lower, ref_upper)

ri_cont_lookup <- cont_ri_all |>
  rename(age_int = age, cont_lower = lower, cont_upper = upper) |>
  select(test, age_int, cont_lower, cont_upper)

ri_guven_lookup <- guven_ri |>
  rename(age_int = age, guven_lower = lower, guven_upper = upper) |>
  select(test, age_int, guven_lower, guven_upper)

fp_detail_raw <- all_data_raw |>
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

fp_by_age_raw <- fp_detail_raw |>
  group_by(test, age_int) |>
  summarise(
    n = n(),
    fp_manufacturer = mean(out_mfr,  na.rm = TRUE) * 100,
    fp_refiner      = mean(out_ref,  na.rm = TRUE) * 100,
    fp_continuous   = mean(out_cont, na.rm = TRUE) * 100,
    fp_guven        = mean(out_guven, na.rm = TRUE) * 100,
    .groups = "drop"
  )

fp_summary_raw <- fp_by_age_raw |>
  group_by(test) |>
  summarise(
    N = sum(n),
    FP_Manufacturer = weighted.mean(fp_manufacturer, n, na.rm = TRUE),
    FP_refineR      = weighted.mean(fp_refiner,      n, na.rm = TRUE),
    FP_Continuous   = weighted.mean(fp_continuous,    n, na.rm = TRUE),
    FP_Guven        = weighted.mean(fp_guven,         n, na.rm = TRUE),
    .groups = "drop"
  )

message("\n=== FP Summary [RAW] ===")
print(fp_summary_raw, n = Inf)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write_csv(fp_by_age_raw,  "data/processed/analysis1_fp_by_age_RAW.csv")
write_csv(fp_summary_raw, "data/processed/analysis1_fp_summary_RAW.csv")
message("[OK] Analysis 1 RAW saved.\n")


# ==============================================================================
# ANALYSIS 2 (RAW): GUVEN TRANSFERABILITY (CLSI EP28-A3c + binomial test)
# ==============================================================================

message("[ANALYSIS 2 RAW] Guven transferability (no 3SD)...")

set.seed(42)
B <- 200
n_validate <- 20
threshold <- 2

validate_one_age_raw <- function(test_val, age_val, lo, hi) {
  d <- all_data_raw |> filter(test == test_val, age_int == age_val)
  n_total <- nrow(d)
  if (n_total < n_validate) {
    return(tibble(n_total = n_total, pass_rate = NA_real_,
                  n_outside = NA_integer_, prop_outside = NA_real_,
                  binom_p = NA_real_, binom_reject = NA_character_))
  }
  passes <- replicate(B, {
    samp <- slice_sample(d, n = n_validate, replace = TRUE)
    n_out <- sum(samp$result_num < lo | samp$result_num > hi)
    n_out <= threshold
  })
  n_out_full <- sum(d$result_num < lo | d$result_num > hi)
  bt <- binom.test(n_out_full, n_total, p = 0.05, alternative = "greater")
  tibble(n_total = n_total, pass_rate = mean(passes) * 100,
         n_outside = n_out_full, prop_outside = n_out_full / n_total * 100,
         binom_p = bt$p.value,
         binom_reject = ifelse(bt$p.value < 0.05, "Yes", "No"))
}

transfer_results_raw <- guven_ri |>
  mutate(res = pmap(list(test, age, lower, upper), validate_one_age_raw)) |>
  unnest(res) |>
  mutate(transferable = ifelse(is.na(pass_rate), NA_character_,
                               ifelse(pass_rate >= 90, "Yes", "No")))

transfer_summary_raw <- transfer_results_raw |>
  group_by(test) |>
  summarise(
    ages_tested = n(),
    ages_pass   = sum(transferable == "Yes", na.rm = TRUE),
    ages_fail   = sum(transferable == "No",  na.rm = TRUE),
    ages_na     = sum(is.na(transferable)),
    mean_pass_rate    = mean(pass_rate, na.rm = TRUE),
    ages_binom_reject = sum(binom_reject == "Yes", na.rm = TRUE),
    .groups = "drop"
  )

message("\n=== Guven Transferability Summary [RAW] ===")
print(transfer_summary_raw)

write_csv(transfer_results_raw, "data/processed/analysis2_guven_transferability_RAW.csv")
message("[OK] Analysis 2 RAW saved.\n")


# ==============================================================================
# ANALYSIS 3 (RAW): FP-CURVE FIGURE
# ==============================================================================

message("[ANALYSIS 3 RAW] FP curve figure (no 3SD)...")

fp_long_raw <- fp_by_age_raw |>
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

make_fp_panel_raw <- function(test_name, tag_label) {
  fp_long_raw |>
    filter(test == test_name, !is.na(fp_rate)) |>
    ggplot(aes(x = age_int, y = fp_rate, colour = ri_source, linetype = ri_source)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 5, linetype = "dotted", colour = PAL$base2) +
    annotate("text", x = 17, y = 5.8, label = "Expected 5%",
             colour = PAL$base2, size = 2.8, hjust = 1) +
    scale_colour_manual(values = c(
      "Manufacturer"        = PAL$base2,
      "Guven et al."        = PAL$accent3,
      "refineR"             = PAL$refiner,
      "Continuous (GAMLSS)" = PAL$contin
    )) +
    scale_linetype_manual(values = c(
      "Manufacturer"        = "dashed",
      "Guven et al."        = "dotdash",
      "refineR"             = "solid",
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

p_fp_pt_r   <- make_fp_panel_raw("PT",         "A")
p_fp_aptt_r <- make_fp_panel_raw("aPTT",       "B")
p_fp_fib_r  <- make_fp_panel_raw("Fibrinogen", "C") +
  theme(legend.position = "bottom")

fig_fp_curves_raw <- p_fp_pt_r / p_fp_aptt_r / p_fp_fib_r +
  plot_annotation(theme = theme_tufte_academic())

save_fig(fig_fp_curves_raw, "figure_fp_curves_by_age_RAW", width = 180, height = 280)
write_csv(fp_long_raw, "data/processed/analysis3_fp_curves_data_RAW.csv")
message("[OK] Analysis 3 RAW saved.\n")


# ==============================================================================
# ANALYSIS 6 (RAW): BOOTSTRAP CI FOR FP RATES
# ==============================================================================

message("[ANALYSIS 6 RAW] Bootstrap CI for FP rates (no 3SD)...")

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
    fp_mean  = mean(boot_rates),
    fp_ci_lo = quantile(boot_rates, 0.025),
    fp_ci_hi = quantile(boot_rates, 0.975)
  )
}

# Manufacturer
fp_boot_mfr_r <- ri_manufacturer |>
  rowwise() |>
  mutate(res = list(boot_fp_ci(all_data_raw, test, lower, upper, B_fp))) |>
  unnest(res) |>
  mutate(ri_source = "Manufacturer") |>
  select(test, ri_source, fp_mean, fp_ci_lo, fp_ci_hi)

# refineR (PT/Fib single pooled; aPTT two age groups joined to data)
fp_boot_ref_pt_fib <- ri_refiner |>
  filter(test %in% c("PT", "Fibrinogen")) |>
  rowwise() |>
  mutate(res = list(boot_fp_ci(all_data_raw, test, lower, upper, B_fp))) |>
  unnest(res) |>
  mutate(ri_source = "refineR") |>
  select(test, ri_source, fp_mean, fp_ci_lo, fp_ci_hi)

# aPTT custom: combine 1-12 and 12-18 partitioned
boot_fp_ci_aptt_part <- function(data) {
  d <- data |> filter(test == "aPTT")
  lo1 <- ri_refiner$lower[ri_refiner$test == "aPTT" & ri_refiner$age_group == "1-12"]
  hi1 <- ri_refiner$upper[ri_refiner$test == "aPTT" & ri_refiner$age_group == "1-12"]
  lo2 <- ri_refiner$lower[ri_refiner$test == "aPTT" & ri_refiner$age_group == "12-18"]
  hi2 <- ri_refiner$upper[ri_refiner$test == "aPTT" & ri_refiner$age_group == "12-18"]
  d <- d |> mutate(lo = ifelse(age_int < 12, lo1, lo2),
                   hi = ifelse(age_int < 12, hi1, hi2),
                   outside = result_num < lo | result_num > hi)
  vals <- d$outside
  n <- length(vals)
  boot_rates <- replicate(B_fp, mean(sample(vals, n, replace = TRUE)) * 100)
  tibble(test = "aPTT", ri_source = "refineR",
         fp_mean = mean(boot_rates),
         fp_ci_lo = quantile(boot_rates, 0.025),
         fp_ci_hi = quantile(boot_rates, 0.975))
}
fp_boot_ref_aptt <- boot_fp_ci_aptt_part(all_data_raw)
fp_boot_ref_r <- bind_rows(fp_boot_ref_pt_fib, fp_boot_ref_aptt)

# Continuous (per-age — outside flag computed per row then bootstrapped)
boot_fp_ci_continuous <- function(data, test_name) {
  d <- data |> filter(test == test_name) |>
    left_join(ri_cont_lookup, by = c("test", "age_int")) |>
    mutate(outside = result_num < cont_lower | result_num > cont_upper) |>
    filter(!is.na(outside))
  vals <- d$outside
  n <- length(vals)
  boot_rates <- replicate(B_fp, mean(sample(vals, n, replace = TRUE)) * 100)
  tibble(test = test_name, ri_source = "Continuous",
         fp_mean = mean(boot_rates),
         fp_ci_lo = quantile(boot_rates, 0.025),
         fp_ci_hi = quantile(boot_rates, 0.975))
}
fp_boot_cont_r <- bind_rows(
  boot_fp_ci_continuous(all_data_raw, "PT"),
  boot_fp_ci_continuous(all_data_raw, "aPTT"),
  boot_fp_ci_continuous(all_data_raw, "Fibrinogen")
)

# Guven (per-age — PT and aPTT only)
boot_guven_fn <- function(tname) {
  d <- all_data_raw |>
    filter(test == tname) |>
    left_join(guven_ri |> rename(age_int = age), by = c("test", "age_int")) |>
    mutate(outside = result_num < lower | result_num > upper) |>
    filter(!is.na(outside))
  vals <- d$outside
  n <- length(vals)
  boot_rates <- replicate(B_fp, mean(sample(vals, n, replace = TRUE)) * 100)
  tibble(test = tname, ri_source = "Guven",
         fp_mean = mean(boot_rates),
         fp_ci_lo = quantile(boot_rates, 0.025),
         fp_ci_hi = quantile(boot_rates, 0.975))
}
fp_boot_guven_r <- bind_rows(boot_guven_fn("PT"), boot_guven_fn("aPTT"))

fp_boot_all_raw <- bind_rows(fp_boot_mfr_r, fp_boot_ref_r, fp_boot_cont_r, fp_boot_guven_r) |>
  mutate(fp_label = sprintf("%.1f (%.1f-%.1f)", fp_mean, fp_ci_lo, fp_ci_hi))

message("\n=== Bootstrap FP CI [RAW] ===")
fp_boot_all_raw |>
  select(test, ri_source, fp_label) |>
  pivot_wider(names_from = ri_source, values_from = fp_label) |>
  print()

write_csv(fp_boot_all_raw, "data/processed/analysis6_fp_bootstrap_ci_RAW.csv")
message("[OK] Analysis 6 RAW saved.\n")


# ==============================================================================
# COMPARISON: 3SD-filtered vs RAW
# ==============================================================================

message("[COMPARE] Building side-by-side comparison...")

# --- FP summary (analysis 1) ---
fp_summary_3sd <- read_csv("data/processed/analysis1_fp_summary.csv", show_col_types = FALSE)
fp_compare <- fp_summary_3sd |>
  rename_with(~paste0(.x, "_3SD"), -test) |>
  left_join(
    fp_summary_raw |> rename_with(~paste0(.x, "_RAW"), -test),
    by = "test"
  ) |>
  mutate(
    delta_N           = N_RAW - N_3SD,
    delta_Manufacturer = FP_Manufacturer_RAW - FP_Manufacturer_3SD,
    delta_refineR      = FP_refineR_RAW      - FP_refineR_3SD,
    delta_Continuous   = FP_Continuous_RAW   - FP_Continuous_3SD,
    delta_Guven        = FP_Guven_RAW        - FP_Guven_3SD
  )

# --- Bootstrap CI (analysis 6) ---
boot_3sd <- read_csv("data/processed/analysis6_fp_bootstrap_ci.csv", show_col_types = FALSE)
boot_compare <- boot_3sd |>
  select(test, ri_source, fp_label_3SD = fp_label,
         fp_mean_3SD = fp_mean, fp_ci_lo_3SD = fp_ci_lo, fp_ci_hi_3SD = fp_ci_hi) |>
  left_join(
    fp_boot_all_raw |>
      select(test, ri_source, fp_label_RAW = fp_label,
             fp_mean_RAW = fp_mean, fp_ci_lo_RAW = fp_ci_lo, fp_ci_hi_RAW = fp_ci_hi),
    by = c("test", "ri_source")
  ) |>
  mutate(delta_mean = fp_mean_RAW - fp_mean_3SD)

# --- Guven transferability (analysis 2) ---
guv_3sd <- read_csv("data/processed/analysis2_guven_transferability.csv", show_col_types = FALSE)
guv_compare <- guv_3sd |>
  select(test, age,
         n_total_3SD     = n_total,
         pass_rate_3SD   = pass_rate,
         prop_out_3SD    = prop_outside,
         binom_p_3SD     = binom_p,
         transferable_3SD = transferable) |>
  left_join(
    transfer_results_raw |>
      select(test, age,
             n_total_RAW      = n_total,
             pass_rate_RAW    = pass_rate,
             prop_out_RAW     = prop_outside,
             binom_p_RAW      = binom_p,
             transferable_RAW = transferable),
    by = c("test", "age")
  ) |>
  mutate(delta_n = n_total_RAW - n_total_3SD,
         delta_pass_rate = pass_rate_RAW - pass_rate_3SD)

# --- Save CSVs ---
write_csv(fp_compare,   "data/processed/COMPARISON_fp_summary_3sd_vs_raw.csv")
write_csv(boot_compare, "data/processed/COMPARISON_fp_bootstrap_3sd_vs_raw.csv")
write_csv(guv_compare,  "data/processed/COMPARISON_guven_3sd_vs_raw.csv")

# --- Build Excel file with 4 sheets ---
out_xlsx <- "data/processed/tables/COMPARISON_3sd_vs_raw.xlsx"
dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)

wb <- createWorkbook()
header_style <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                            halign = "center", border = "bottom")

add_sheet <- function(name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = header_style)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}

# Sheet 0: input sizes summary
sizes_df <- tibble(
  test = c("PT", "aPTT", "Fibrinogen", "TOTAL"),
  n_RAW = c(sum(all_data_raw$test == "PT"),
            sum(all_data_raw$test == "aPTT"),
            sum(all_data_raw$test == "Fibrinogen"),
            nrow(all_data_raw)),
  n_3SD = c(sum(all_data_3sd$test == "PT"),
            sum(all_data_3sd$test == "aPTT"),
            sum(all_data_3sd$test == "Fibrinogen"),
            nrow(all_data_3sd))
) |>
  mutate(delta = n_RAW - n_3SD,
         pct_removed = round(delta / n_RAW * 100, 2))

add_sheet("Sample_Sizes",       sizes_df)
add_sheet("FP_Summary",         fp_compare)
add_sheet("FP_Bootstrap_CI",    boot_compare)
add_sheet("Guven_Transfer",     guv_compare)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message(sprintf("[OK] Comparison Excel saved: %s", out_xlsx))

message("\n============================================================")
message("[DONE] RAW analyses + comparison complete.")
message("============================================================")
