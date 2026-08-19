# ==============================================================================
# ANALYSIS 4 - DEPARTMENT-BASED EXCLUSION SENSITIVITY
# Scenarios:
#   S1 baseline   : all departments (current manuscript)
#   S2 no ER      : exclude COCUK ACIL
#   S3 outpatient : only departments matching "POL." (polyclinics)
#   S4 ER only    : only COCUK ACIL
# refineR run with NBootstrap = 1 (point estimates only; no CI for now).
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(janitor); library(readr); library(tidyr)
  library(openxlsx); library(refineR)
})

if (!exists("PAL")) source("src/common.R")

main_files <- c(
  PT         = "data/coa_results/tce 2026 1-18 pt.xlsx",
  aPTT       = "data/coa_results/tce 2026 1-18 aptt.xlsx",
  Fibrinogen = "data/coa_results/tce 2026 1-18 fib.xlsx"
)

load_main <- function(path, tname) {
  suppressMessages(read_excel(path)) |>
    clean_names() |>
    mutate(test = tname,
           age_int = floor(age_year))
}

raw <- bind_rows(
  load_main(main_files["PT"],         "PT"),
  load_main(main_files["aPTT"],       "aPTT"),
  load_main(main_files["Fibrinogen"], "Fibrinogen")
)

# --- Classify department -----------------------------------------------------
classify_dept <- function(d) {
  case_when(
    grepl("COCUK ACIL", d, ignore.case = TRUE)              ~ "ER",
    grepl("\\bPOL\\.?\\b|POLIKLINIK", d, ignore.case = TRUE) ~ "Outpatient",
    grepl("SERVIS|KLINIK|MT1|YBU|YOGUN", d, ignore.case = TRUE) ~ "Inpatient",
    TRUE                                                    ~ "Other"
  )
}

raw <- raw |> mutate(dept_cat = classify_dept(department))

cat("\n=== DEPARTMENT CATEGORY DISTRIBUTION (per test) ===\n")
dept_dist <- raw |>
  group_by(test, dept_cat) |>
  summarise(n = n(), .groups = "drop") |>
  pivot_wider(names_from = dept_cat, values_from = n, values_fill = 0) |>
  mutate(Total = ER + Outpatient + Inpatient + Other)
print(dept_dist)

# --- Define scenarios --------------------------------------------------------
scenarios <- list(
  "S1_baseline"   = function(df) df,
  "S2_no_ER"      = function(df) filter(df, dept_cat != "ER"),
  "S3_outpatient" = function(df) filter(df, dept_cat == "Outpatient"),
  "S4_ER_only"    = function(df) filter(df, dept_cat == "ER")
)

# --- refineR helper (NBootstrap = 1, point estimate only) -------------------
refiner_point <- function(values, test_name, scenario, age_grp = NA) {
  n <- length(values)
  if (n < 100) {
    return(tibble(test = test_name, scenario = scenario, age_group = age_grp,
                  n = n, lower = NA_real_, median = NA_real_, upper = NA_real_,
                  note = "n<100 skipped"))
  }
  fit <- try(refineR::findRI(values, NBootstrap = 1), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(tibble(test = test_name, scenario = scenario, age_group = age_grp,
                  n = n, lower = NA_real_, median = NA_real_, upper = NA_real_,
                  note = "refineR failed"))
  }
  ri <- as.data.frame(getRI(fit, RIperc = c(0.025, 0.50, 0.975)))
  tibble(test = test_name, scenario = scenario, age_group = age_grp,
         n = n,
         lower  = round(ri$PointEst[ri$Percentile == 0.025], 2),
         median = round(ri$PointEst[ri$Percentile == 0.500], 2),
         upper  = round(ri$PointEst[ri$Percentile == 0.975], 2),
         note   = "")
}

# --- Run refineR per scenario per test --------------------------------------
results <- list()
set.seed(42)

for (sname in names(scenarios)) {
  cat(sprintf("\n[Scenario %s] running ...\n", sname))
  d_scn <- scenarios[[sname]](raw)

  # PT
  d <- d_scn |> filter(test == "PT")
  results[[paste(sname, "PT", sep="_")]] <-
    refiner_point(d$result_num, "PT", sname, "1-18 y")
  cat(sprintf("  PT n=%d\n", nrow(d)))

  # aPTT 1-12y partition
  d <- d_scn |> filter(test == "aPTT", age_int < 12)
  results[[paste(sname, "aPTT_1_12", sep="_")]] <-
    refiner_point(d$result_num, "aPTT", sname, "1-12 y")
  cat(sprintf("  aPTT 1-12y n=%d\n", nrow(d)))

  # aPTT 12-18y partition
  d <- d_scn |> filter(test == "aPTT", age_int >= 12)
  results[[paste(sname, "aPTT_12_18", sep="_")]] <-
    refiner_point(d$result_num, "aPTT", sname, "12-18 y")
  cat(sprintf("  aPTT 12-18y n=%d\n", nrow(d)))

  # Fibrinogen 1-18y
  d <- d_scn |> filter(test == "Fibrinogen")
  results[[paste(sname, "Fib", sep="_")]] <-
    refiner_point(d$result_num, "Fibrinogen", sname, "1-18 y")
  cat(sprintf("  Fib n=%d\n", nrow(d)))
}

a4 <- bind_rows(results)

cat("\n=== ANALYSIS 4 RESULTS ===\n")
print(a4 |> arrange(test, age_group, scenario), n = Inf)

# --- FP rates per scenario (using baseline manuscript RI, not the scenario RI)
# This shows whether different sub-populations would have different FP rates
# under the published RI.
ri_lookup <- bind_rows(
  ri_refiner |> filter(test == "PT")        |> mutate(group_filter = "PT_all"),
  ri_refiner |> filter(test == "aPTT", age_group == "1-12")  |> mutate(group_filter = "aPTT_1_12"),
  ri_refiner |> filter(test == "aPTT", age_group == "12-18") |> mutate(group_filter = "aPTT_12_18"),
  ri_refiner |> filter(test == "Fibrinogen")|> mutate(group_filter = "Fib_all")
) |>
  select(group_filter, baseline_lower = lower, baseline_upper = upper)

fp_by_scenario <- bind_rows(
  raw |> mutate(group_filter = "PT_all", grp = "PT")        |> filter(test == "PT"),
  raw |> filter(test == "aPTT", age_int < 12)  |> mutate(group_filter = "aPTT_1_12",  grp = "aPTT_1_12"),
  raw |> filter(test == "aPTT", age_int >= 12) |> mutate(group_filter = "aPTT_12_18", grp = "aPTT_12_18"),
  raw |> filter(test == "Fibrinogen")          |> mutate(group_filter = "Fib_all",    grp = "Fib")
) |>
  left_join(ri_lookup, by = "group_filter") |>
  mutate(outside = result_num < baseline_lower | result_num > baseline_upper)

fp_summary <- bind_rows(
  fp_by_scenario |> mutate(scenario = "S1_baseline"),
  fp_by_scenario |> filter(dept_cat != "ER")        |> mutate(scenario = "S2_no_ER"),
  fp_by_scenario |> filter(dept_cat == "Outpatient")|> mutate(scenario = "S3_outpatient"),
  fp_by_scenario |> filter(dept_cat == "ER")        |> mutate(scenario = "S4_ER_only")
) |>
  group_by(scenario, grp) |>
  summarise(n = n(), fp_pct = round(mean(outside) * 100, 2), .groups = "drop") |>
  pivot_wider(names_from = grp, values_from = c(n, fp_pct))

cat("\n=== FP RATE AT BASELINE RI ACROSS SCENARIOS ===\n")
print(fp_summary, width = 220)

# --- Save outputs ------------------------------------------------------------
write_csv(a4,         "data/processed/analysis4_sensitivity_RI.csv")
write_csv(fp_summary, "data/processed/analysis4_sensitivity_FP.csv")
write_csv(dept_dist,  "data/processed/analysis4_department_distribution.csv")

# Excel supplement
out_xlsx <- "data/processed/tables/SUPPLEMENT_analysis4_sensitivity.xlsx"
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                  halign = "center", border = "bottom")
add_sheet <- function(nm, df) {
  addWorksheet(wb, nm); writeData(wb, nm, df, headerStyle = hs)
  freezePane(wb, nm, firstRow = TRUE)
  setColWidths(wb, nm, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet("Dept_Distribution", dept_dist)
add_sheet("RI_by_Scenario",    a4)
add_sheet("FP_by_Scenario",    fp_summary)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("\n[OK] Saved: %s\n", out_xlsx))
