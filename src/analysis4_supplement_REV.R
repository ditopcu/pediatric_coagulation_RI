# ==============================================================================
# PATIENT-ORIGIN SUPPLEMENT (reviewer R2-1): outpatient vs inpatient contrast
# ==============================================================================
# Builds the supplementary table requested by Reviewer 2 (were outpatient vs
# inpatient analyses performed?). Reuses the department classification of
# analysis4_sensitivity.R and adds the inpatient-only scenario that the
# original scenario set (baseline / no-ER / outpatient-only / ER-only) lacked.
#
# refineR is run with NBootstrap = 1 (point estimates only), matching
# analysis4_sensitivity.R. Main manuscript RIs are NOT changed.
#
# Outputs (new files only):
#   data/processed/REVIEW_A4_patient_origin_RI_REV.csv
#   data/processed/REVIEW_A4_patient_origin_flagging_REV.csv
#   data/processed/tables/SUPPLEMENT_patient_origin_REV.xlsx
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

# Same classification as analysis4_sensitivity.R
classify_dept <- function(d) {
  case_when(
    grepl("COCUK ACIL", d, ignore.case = TRUE)                  ~ "ER",
    grepl("\\bPOL\\.?\\b|POLIKLINIK", d, ignore.case = TRUE)    ~ "Outpatient",
    grepl("SERVIS|KLINIK|MT1|YBU|YOGUN", d, ignore.case = TRUE) ~ "Inpatient",
    TRUE                                                        ~ "Other"
  )
}

raw <- raw |> mutate(dept_cat = classify_dept(department))

# --- Category counts per analyte ---------------------------------------------
origin_counts <- raw |>
  group_by(test, dept_cat) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(test) |>
  mutate(pct = round(n / sum(n) * 100, 1)) |>
  ungroup() |>
  arrange(test, dept_cat)

cat("\n=== PATIENT ORIGIN COUNTS ===\n")
print(origin_counts, n = Inf)

# --- Scenarios for the reviewer contrast -------------------------------------
scenarios <- list(
  "S1_baseline"        = function(df) df,
  "S3_outpatient_only" = function(df) filter(df, dept_cat == "Outpatient"),
  "S5_inpatient_only"  = function(df) filter(df, dept_cat == "Inpatient")
)

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

results <- list()
set.seed(42)

for (sname in names(scenarios)) {
  cat(sprintf("\n[Scenario %s] running ...\n", sname))
  d_scn <- scenarios[[sname]](raw)

  d <- d_scn |> filter(test == "PT")
  results[[paste(sname, "PT", sep = "_")]] <-
    refiner_point(d$result_num, "PT", sname, "1-18 y")

  d <- d_scn |> filter(test == "aPTT", age_int < 12)
  results[[paste(sname, "aPTT_1_12", sep = "_")]] <-
    refiner_point(d$result_num, "aPTT", sname, "1-12 y")

  d <- d_scn |> filter(test == "aPTT", age_int >= 12)
  results[[paste(sname, "aPTT_12_18", sep = "_")]] <-
    refiner_point(d$result_num, "aPTT", sname, "12-18 y")

  d <- d_scn |> filter(test == "Fibrinogen")
  results[[paste(sname, "Fib", sep = "_")]] <-
    refiner_point(d$result_num, "Fibrinogen", sname, "1-18 y")
}

ri_contrast <- bind_rows(results) |>
  arrange(test, age_group, scenario)

cat("\n=== RI CONTRAST (baseline vs outpatient-only vs inpatient-only) ===\n")
print(ri_contrast, n = Inf)

# --- Flagging at the published (baseline) RI, by patient origin --------------
ri_lookup <- bind_rows(
  ri_refiner |> filter(test == "PT")                          |> mutate(grp = "PT 1-18 y"),
  ri_refiner |> filter(test == "aPTT", age_group == "1-12")   |> mutate(grp = "aPTT 1-12 y"),
  ri_refiner |> filter(test == "aPTT", age_group == "12-18")  |> mutate(grp = "aPTT 12-18 y"),
  ri_refiner |> filter(test == "Fibrinogen")                  |> mutate(grp = "Fibrinogen 1-18 y")
) |>
  select(grp, baseline_lower = lower, baseline_upper = upper)

flag_base <- bind_rows(
  raw |> filter(test == "PT")                  |> mutate(grp = "PT 1-18 y"),
  raw |> filter(test == "aPTT", age_int < 12)  |> mutate(grp = "aPTT 1-12 y"),
  raw |> filter(test == "aPTT", age_int >= 12) |> mutate(grp = "aPTT 12-18 y"),
  raw |> filter(test == "Fibrinogen")          |> mutate(grp = "Fibrinogen 1-18 y")
) |>
  left_join(ri_lookup, by = "grp") |>
  mutate(outside = result_num < baseline_lower | result_num > baseline_upper)

flagging_by_origin <- bind_rows(
  flag_base |> mutate(origin = dept_cat),
  flag_base |> mutate(origin = "All")
) |>
  group_by(grp, origin) |>
  summarise(n = n(),
            flag_pct = round(mean(outside) * 100, 2),
            .groups = "drop") |>
  arrange(grp, origin)

cat("\n=== FLAGGING AT PUBLISHED RI, BY PATIENT ORIGIN ===\n")
print(flagging_by_origin, n = Inf)

# --- Save --------------------------------------------------------------------
write_csv(ri_contrast,        "data/processed/REVIEW_A4_patient_origin_RI_REV.csv")
write_csv(flagging_by_origin, "data/processed/REVIEW_A4_patient_origin_flagging_REV.csv")

out_xlsx <- "data/processed/tables/SUPPLEMENT_patient_origin_REV.xlsx"
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                  halign = "center", border = "bottom")
add_sheet <- function(nm, df) {
  addWorksheet(wb, nm); writeData(wb, nm, df, headerStyle = hs)
  freezePane(wb, nm, firstRow = TRUE)
  setColWidths(wb, nm, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet("Origin_Counts",      origin_counts)
add_sheet("RI_Contrast",        ri_contrast)
add_sheet("Flagging_by_Origin", flagging_by_origin)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("\n[OK] Saved: %s\n", out_xlsx))
