# ==============================================================================
# INFANT ANALYSIS - STEP A
# Descriptive statistics for infants (30-365 days), pooled across sex.
# No outlier filtering. No CI.
# Output: data/processed/infant_A_descriptive.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(janitor); library(readr)
  if (!requireNamespace("moments", quietly = TRUE))
    install.packages("moments", repos = "https://cran.r-project.org")
  library(moments)
})

LO_Y <- 30 / 365
HI_Y <- 1.0

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

# --- Descriptive table (pooled across sex) -----------------------------------
desc <- infant |>
  group_by(test) |>
  summarise(
    n        = n(),
    n_female = sum(sex == "K", na.rm = TRUE),
    n_male   = sum(sex == "E", na.rm = TRUE),
    pct_female = round(n_female / n * 100, 1),
    age_days_median = round(median(age_days), 0),
    age_days_min    = round(min(age_days), 0),
    age_days_max    = round(max(age_days), 0),
    mean   = round(mean(result_num), 2),
    sd     = round(sd(result_num),   2),
    median = round(median(result_num), 2),
    q25    = round(quantile(result_num, 0.25), 2),
    q75    = round(quantile(result_num, 0.75), 2),
    min    = round(min(result_num), 2),
    max    = round(max(result_num), 2),
    skewness = round(moments::skewness(result_num), 2),
    kurtosis = round(moments::kurtosis(result_num), 2),
    .groups = "drop"
  ) |>
  mutate(test = factor(test, levels = c("PT", "aPTT", "Fibrinogen"))) |>
  arrange(test)

cat("\n=== INFANT DESCRIPTIVE STATISTICS (30-365 days, pooled) ===\n")
print(desc, width = 200)

write_csv(desc, "data/processed/infant_A_descriptive.csv")
cat("\n[OK] Saved: data/processed/infant_A_descriptive.csv\n")
