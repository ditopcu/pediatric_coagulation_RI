# ==============================================================================
# INFANT ANALYSIS - STEP C
# Classify 30-365 day infant observations against the current pediatric (1-18 y)
# refineR RIs. NO outlier filtering, NO CI computation.
# Output: data/processed/infant_C_classification.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(janitor); library(readr); library(tidyr)
})

LO_Y <- 30 / 365     # 30 days  ~ 0.0822 years
HI_Y <- 1.0          # 365 days = 1.0  years

# --- Load raw legacy files (no 3SD filter, no department filter) -------------
load_legacy <- function(path, test_name) {
  read_excel(path) |>
    clean_names() |>
    select(any_of(c("sample_id", "sex", "age_year", "result_num"))) |>
    mutate(test = test_name)
}

infant <- bind_rows(
  load_legacy("data/coa_results/PT-2 2025.10.17.xlsx",   "PT"),
  load_legacy("data/coa_results/aPTT 2025.10.17.xlsx",   "aPTT"),
  load_legacy("data/coa_results/Fib-2 2025.10.17.xlsx",  "Fibrinogen")
) |>
  filter(age_year >= LO_Y, age_year < HI_Y) |>
  mutate(age_days = round(age_year * 365))

cat(sprintf("\n[INFANT] Subset 30-365 days: %d rows\n", nrow(infant)))
cat("Per test counts:\n")
print(infant |> count(test))

# --- Pediatric (1-18 y) refineR RI lookup ------------------------------------
# PT and Fibrinogen: single pooled RI for 1-18 y
# aPTT: 1-12 y partition (closest to infants); also report 12-18 y for reference
peds_ri <- tribble(
  ~test,        ~ri_label,         ~lower, ~upper,
  "PT",         "1-18 y pooled",    8.56,   10.90,
  "aPTT",       "1-12 y partition", 24.30,  35.30,
  "Fibrinogen", "1-18 y pooled",    1.89,   3.79
)

# --- Classify each observation -----------------------------------------------
infant_class <- infant |>
  left_join(peds_ri, by = "test") |>
  mutate(
    classification = case_when(
      result_num <  lower ~ "Below",
      result_num >  upper ~ "Above",
      TRUE                ~ "Within"
    )
  )

# --- Summary table per test --------------------------------------------------
summary_tbl <- infant_class |>
  group_by(test, ri_label, lower, upper) |>
  summarise(
    n_infant  = n(),
    n_below   = sum(classification == "Below"),
    n_within  = sum(classification == "Within"),
    n_above   = sum(classification == "Above"),
    pct_below = round(n_below  / n_infant * 100, 1),
    pct_within= round(n_within / n_infant * 100, 1),
    pct_above = round(n_above  / n_infant * 100, 1),
    pct_outside = round((n_below + n_above) / n_infant * 100, 1),
    .groups = "drop"
  )

cat("\n=== INFANT CLASSIFICATION VS PEDIATRIC refineR RI ===\n")
print(summary_tbl)

# --- Median (IQR) of infant raw result_num for quick context ----------------
quick_desc <- infant |>
  group_by(test) |>
  summarise(
    n = n(),
    median = round(median(result_num), 2),
    q25 = round(quantile(result_num, 0.25), 2),
    q75 = round(quantile(result_num, 0.75), 2),
    min = round(min(result_num), 2),
    max = round(max(result_num), 2),
    .groups = "drop"
  )

cat("\n=== INFANT result_num summary (raw, no outlier filter) ===\n")
print(quick_desc)

# --- Save ---------------------------------------------------------------------
write_csv(summary_tbl, "data/processed/infant_C_classification.csv")
write_csv(quick_desc,  "data/processed/infant_C_descriptive_quick.csv")
cat("\n[OK] Saved:\n")
cat("  data/processed/infant_C_classification.csv\n")
cat("  data/processed/infant_C_descriptive_quick.csv\n")
