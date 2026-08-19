# ==============================================================================
# INFANT ANALYSIS - STEP B
# refineR RI estimation for 30-365 day infants (PT, aPTT, Fibrinogen).
# NBootstrap = 200. No outlier filtering (refineR handles pathological tail).
# Saves RDS + CSV + comparison Excel + comparison figure.
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(janitor); library(readr); library(tidyr)
  library(ggplot2); library(patchwork); library(openxlsx)
  library(refineR)
})

if (!exists("PAL")) source("src/common.R")

LO_Y <- 30 / 365
HI_Y <- 1.0
NBOOT <- 200

# --- Load infant subset -------------------------------------------------------
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
  filter(age_year >= LO_Y, age_year < HI_Y)

cat(sprintf("\n[INFANT B] Loaded %d rows (PT=%d, aPTT=%d, Fib=%d)\n",
            nrow(infant),
            sum(infant$test == "PT"),
            sum(infant$test == "aPTT"),
            sum(infant$test == "Fibrinogen")))

# --- refineR per test ---------------------------------------------------------
run_refiner <- function(test_name) {
  d <- infant |> filter(test == test_name)
  cat(sprintf("\n[REFINER] %s (n=%d) -- NBootstrap=%d ...\n",
              test_name, nrow(d), NBOOT))
  t0 <- proc.time()
  res <- refineR::findRI(d$result_num, NBootstrap = NBOOT)
  dt <- (proc.time() - t0)["elapsed"]
  cat(sprintf("[REFINER] %s done in %.1f s\n", test_name, dt))
  res
}

set.seed(42)
fit_pt   <- run_refiner("PT")
fit_aptt <- run_refiner("aPTT")
fit_fib  <- run_refiner("Fibrinogen")

# --- Save RDS ----------------------------------------------------------------
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(fit_pt,   "data/processed/infant_refineR_PT.RDS")
saveRDS(fit_aptt, "data/processed/infant_refineR_aPTT.RDS")
saveRDS(fit_fib,  "data/processed/infant_refineR_Fibrinogen.RDS")
cat("[OK] RDS files saved.\n")

# --- Extract RI + 90% CI -----------------------------------------------------
extract_ri <- function(fit, test_name, n_sub) {
  out <- getRI(fit, RIperc = c(0.025, 0.50, 0.975), CIprop = 0.90,
               pi.type = "tolerance")
  # getRI returns a data frame with columns Percentile, PointEst, CILow, CIHigh
  out <- as.data.frame(out)
  tibble(
    test = test_name,
    n = n_sub,
    lower = out$PointEst[out$Percentile == 0.025],
    lower_ci_lo = out$CILow [out$Percentile == 0.025],
    lower_ci_hi = out$CIHigh[out$Percentile == 0.025],
    median = out$PointEst[out$Percentile == 0.500],
    upper = out$PointEst[out$Percentile == 0.975],
    upper_ci_lo = out$CILow [out$Percentile == 0.975],
    upper_ci_hi = out$CIHigh[out$Percentile == 0.975]
  )
}

infant_ri <- bind_rows(
  extract_ri(fit_pt,   "PT",         sum(infant$test == "PT")),
  extract_ri(fit_aptt, "aPTT",       sum(infant$test == "aPTT")),
  extract_ri(fit_fib,  "Fibrinogen", sum(infant$test == "Fibrinogen"))
)

# Round for display
infant_ri_disp <- infant_ri |>
  mutate(across(c(lower, lower_ci_lo, lower_ci_hi,
                  median, upper, upper_ci_lo, upper_ci_hi),
                ~round(.x, 2)))

cat("\n=== INFANT refineR RI (30-365 days, NBootstrap=200) ===\n")
print(infant_ri_disp)

# --- Build pediatric vs infant comparison table ------------------------------
# Pediatric pooled RI from common.R (ri_refiner)
peds_ri <- bind_rows(
  ri_refiner |> filter(test == "PT") |> mutate(age_group = "1-18 y"),
  ri_refiner |> filter(test == "aPTT", age_group == "1-12") |>
    mutate(age_group = "1-12 y (pediatric reference partition)"),
  ri_refiner |> filter(test == "Fibrinogen") |> mutate(age_group = "1-18 y")
) |>
  select(test, age_group, peds_lower = lower, peds_upper = upper) |>
  mutate(across(c(peds_lower, peds_upper), ~round(.x, 2)))

compare_tbl <- infant_ri_disp |>
  select(test, n_infant = n,
         infant_lower = lower, infant_lower_ci_lo = lower_ci_lo, infant_lower_ci_hi = lower_ci_hi,
         infant_upper = upper, infant_upper_ci_lo = upper_ci_lo, infant_upper_ci_hi = upper_ci_hi) |>
  left_join(peds_ri, by = "test") |>
  select(test, age_group_pediatric = age_group,
         peds_lower, peds_upper,
         n_infant,
         infant_lower, infant_lower_ci_lo, infant_lower_ci_hi,
         infant_upper, infant_upper_ci_lo, infant_upper_ci_hi)

cat("\n=== PEDIATRIC vs INFANT COMPARISON ===\n")
print(compare_tbl, width = 250)

write_csv(infant_ri_disp, "data/processed/infant_B_refineR_RI.csv")
write_csv(compare_tbl,    "data/processed/infant_B_pediatric_vs_infant.csv")
cat("[OK] CSVs saved.\n")

# --- Comparison figure: RI bands side-by-side --------------------------------
fig_data <- bind_rows(
  compare_tbl |> select(test, group = "Pediatric (current study)",
                        lower = peds_lower, upper = peds_upper) |>
    mutate(group = "Pediatric (current study)",
           lower_ci_lo = NA_real_, lower_ci_hi = NA_real_,
           upper_ci_lo = NA_real_, upper_ci_hi = NA_real_),
  compare_tbl |> select(test, lower = infant_lower, upper = infant_upper,
                        lower_ci_lo = infant_lower_ci_lo, lower_ci_hi = infant_lower_ci_hi,
                        upper_ci_lo = infant_upper_ci_lo, upper_ci_hi = infant_upper_ci_hi) |>
    mutate(group = "Infant (30-365 days)")
) |>
  mutate(group = factor(group, levels = c("Pediatric (current study)",
                                          "Infant (30-365 days)")),
         test = factor(test, levels = c("PT", "aPTT", "Fibrinogen")))

make_panel <- function(test_name, y_label, tag_label) {
  df <- fig_data |> filter(test == test_name)
  ggplot(df, aes(x = group, y = lower)) +
    geom_linerange(aes(ymin = lower, ymax = upper, colour = group),
                   linewidth = 6, alpha = 0.55) +
    geom_point(aes(y = lower, colour = group), size = 2.6) +
    geom_point(aes(y = upper, colour = group), size = 2.6) +
    geom_errorbar(aes(ymin = lower_ci_lo, ymax = lower_ci_hi, colour = group),
                  width = 0.18, linewidth = 0.6, na.rm = TRUE) +
    geom_errorbar(aes(ymin = upper_ci_lo, ymax = upper_ci_hi, colour = group),
                  width = 0.18, linewidth = 0.6, na.rm = TRUE) +
    scale_colour_manual(values = c(
      "Pediatric (current study)" = PAL$refiner,
      "Infant (30-365 days)"      = PAL$accent3
    )) +
    labs(x = NULL, y = y_label, tag = tag_label) +
    theme_tufte_academic() +
    theme(plot.tag = element_text(size = 14, face = "bold"),
          legend.position = "none",
          axis.text.x = element_text(size = 9))
}

p_pt   <- make_panel("PT",         "PT, seconds",     "A")
p_aptt <- make_panel("aPTT",       "aPTT, seconds",   "B")
p_fib  <- make_panel("Fibrinogen", "Fibrinogen, g/L", "C")

fig <- (p_pt | p_aptt | p_fib) +
  plot_annotation(theme = theme_tufte_academic())

save_fig(fig, "figure_supplement_infant_RI", width = 220, height = 100)
cat("[OK] Figure saved.\n")

# --- Supplement Excel --------------------------------------------------------
out_xlsx <- "data/processed/tables/SUPPLEMENT_infant_30_365d.xlsx"
dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                  halign = "center", border = "bottom")
add_sheet <- function(nm, df) {
  addWorksheet(wb, nm); writeData(wb, nm, df, headerStyle = hs)
  freezePane(wb, nm, firstRow = TRUE)
  setColWidths(wb, nm, cols = seq_len(ncol(df)), widths = "auto")
}

desc <- read_csv("data/processed/infant_A_descriptive.csv", show_col_types = FALSE)
class_tbl <- read_csv("data/processed/infant_C_classification.csv", show_col_types = FALSE)

add_sheet("Descriptive",          desc)
add_sheet("Classification_vs_RI", class_tbl)
add_sheet("Infant_refineR_RI",    infant_ri_disp)
add_sheet("Pediatric_vs_Infant",  compare_tbl)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("[OK] Supplement Excel saved: %s\n", out_xlsx))

cat("\n=========================================================\n")
cat("[DONE] Infant analysis (steps A + B + C) complete.\n")
cat("=========================================================\n")
