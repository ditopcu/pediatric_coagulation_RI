# ==============================================================================
# BINARY PATIENT ORIGIN: OUTPATIENT vs INPATIENT (reviewer R2-1)
# ==============================================================================
# Reviewer item 1 asks a binary question -- was an outpatient versus inpatient
# analysis performed. analysis4_sensitivity.R and analysis4_supplement_REV.R
# answer it with a four-way split (ER / Outpatient / Inpatient / Other); this
# script produces the two-group contrast the reviewer asked for.
#
# Classification of the `department` string, applied in this order:
#   1. contains KONSULTASYON              -> Inpatient
#      The field records the consulting service, not where the patient was.
#      Neither the analysis sets nor the raw LIS archive carry a patient
#      location column, so this cannot be resolved from the data.
#   2. contains ACIL                      -> Outpatient
#      Emergency attendances are treated as ambulatory.
#   3. POL. / POLIKLI... / MUAYENE        -> Outpatient
#      POLIKLI also catches the Turkish suffixed form POLIKLINIGI and the
#      string truncated at 50 characters in the source file.
#   4. SERVIS / SERV / KLINIK / KLINIGI / MT1-4 floor code -> Inpatient
#      SERV and KLINIGI catch the abbreviated and suffixed ward names that a
#      plain SERVIS / KLINIK pattern misses.
#   5. remainder                          -> Outpatient
#      Numbered specialty examination units that carry no clinic or ward word
#      at all (KULAK BURUN BOGAZ HASTALIKLARI -4, SASILIK - 1, DERMATOLOJI - 7).
#      Written out in full to the department map so the bucket is auditable.
#   COCUK NOROLOJI (SK) is held out as Unclassified: the SK suffix is not
#   resolved, and it is 3 rows in PT and aPTT, 0 in fibrinogen.
#
# refineR is run with NBootstrap = 1 (point estimates only), matching
# analysis4_sensitivity.R and analysis4_supplement_REV.R. Main manuscript RIs
# are NOT changed.
#
# Outputs (new files only):
#   data/processed/REVIEW_A4_binary_origin_RI_REV.csv
#   data/processed/REVIEW_A4_binary_origin_flagging_REV.csv
#   data/processed/REVIEW_A4_binary_origin_deptmap_REV.csv
#   data/processed/tables/SUPPLEMENT_binary_origin_REV.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(janitor); library(readr); library(tidyr)
  library(tibble); library(openxlsx); library(refineR)
})

# Reference limits come from src/common.R, which is the single source of truth
# for `ri_refiner`, `ri_manufacturer` and `ri_direct`. Sourcing it also loads
# the analysis data and the GAMLSS continuous curves from cont_out/, so that
# directory has to hold a complete run before this script will start.
source("src/common.R")

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

# --- Binary classification ---------------------------------------------------
classify_origin <- function(d) {
  u <- toupper(d)
  case_when(
    grepl("\\(SK\\)", u, perl = TRUE)                       ~ "Unclassified",
    grepl("KONSULTASYON", u, perl = TRUE)                   ~ "Inpatient",
    grepl("\\bACIL\\b", u, perl = TRUE)                     ~ "Outpatient",
    grepl("\\bPOL\\.?\\b|POLIKLI|MUAYENE", u, perl = TRUE)  ~ "Outpatient",
    grepl("SERVIS|\\bSERV\\b|\\bKLINIK|KLINIGI|MT[1-4]-", u, perl = TRUE) ~ "Inpatient",
    TRUE                                                    ~ "Outpatient"
  )
}

# The rule that assigned each row, kept for the audit trail.
which_rule <- function(d) {
  u <- toupper(d)
  case_when(
    grepl("\\(SK\\)", u, perl = TRUE)                       ~ "R0_unresolved_SK",
    grepl("KONSULTASYON", u, perl = TRUE)                   ~ "R1_konsultasyon",
    grepl("\\bACIL\\b", u, perl = TRUE)                     ~ "R2_acil",
    grepl("\\bPOL\\.?\\b|POLIKLI|MUAYENE", u, perl = TRUE)  ~ "R3_poliklinik_muayene",
    grepl("SERVIS|\\bSERV\\b|\\bKLINIK|KLINIGI|MT[1-4]-", u, perl = TRUE) ~ "R4_servis_klinik",
    TRUE                                                    ~ "R5_numarali_brans_birimi"
  )
}

raw <- raw |>
  mutate(origin = classify_origin(department),
         origin_rule = which_rule(department))

cat("\n=== BINARY ORIGIN DISTRIBUTION (per test) ===\n")
origin_dist <- raw |>
  group_by(test, origin) |>
  summarise(n = n(), .groups = "drop") |>
  pivot_wider(names_from = origin, values_from = n, values_fill = 0)
print(origin_dist)

cat("\n=== ROWS PER CLASSIFICATION RULE ===\n")
rule_dist <- raw |>
  group_by(test, origin_rule, origin) |>
  summarise(n = n(), .groups = "drop") |>
  arrange(test, origin_rule)
print(rule_dist, n = Inf)

# Department -> bucket map, so no department is silently absorbed.
dept_map <- raw |>
  group_by(test, department, origin_rule, origin) |>
  summarise(n = n(), .groups = "drop") |>
  arrange(test, origin, desc(n))

cat("\n=== RULE 5 (no clinic or ward word in the name) -- full contents ===\n")
print(dept_map |> filter(origin_rule == "R5_numarali_brans_birimi") |>
        select(test, department, n), n = Inf)

# --- Scenarios ---------------------------------------------------------------
scenarios <- list(
  "S1_all"          = function(df) df,
  "S2_outpatient"   = function(df) filter(df, origin == "Outpatient"),
  "S3_inpatient"    = function(df) filter(df, origin == "Inpatient")
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
  cat(sprintf("  PT n=%d\n", nrow(d)))

  d <- d_scn |> filter(test == "aPTT", age_int < 12)
  results[[paste(sname, "aPTT_1_12", sep = "_")]] <-
    refiner_point(d$result_num, "aPTT", sname, "1-12 y")
  cat(sprintf("  aPTT 1-12y n=%d\n", nrow(d)))

  d <- d_scn |> filter(test == "aPTT", age_int >= 12)
  results[[paste(sname, "aPTT_12_18", sep = "_")]] <-
    refiner_point(d$result_num, "aPTT", sname, "12-18 y")
  cat(sprintf("  aPTT 12-18y n=%d\n", nrow(d)))

  d <- d_scn |> filter(test == "Fibrinogen")
  results[[paste(sname, "Fib", sep = "_")]] <-
    refiner_point(d$result_num, "Fibrinogen", sname, "1-18 y")
  cat(sprintf("  Fib n=%d\n", nrow(d)))
}

ri_contrast <- bind_rows(results) |> arrange(test, age_group, scenario)

cat("\n=== RI CONTRAST (all vs outpatient vs inpatient) ===\n")
print(ri_contrast, n = Inf)

# --- Flagging at the published refineR RI, by origin -------------------------
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
  flag_base |> mutate(origin_lbl = origin),
  flag_base |> mutate(origin_lbl = "All")
) |>
  group_by(grp, origin_lbl) |>
  summarise(n = n(),
            flag_pct = round(mean(outside) * 100, 2),
            .groups = "drop") |>
  arrange(grp, origin_lbl)

cat("\n=== FLAGGING AT PUBLISHED RI, BY BINARY ORIGIN ===\n")
print(flagging_by_origin, n = Inf)

# --- Save --------------------------------------------------------------------
write_csv(ri_contrast,        "data/processed/REVIEW_A4_binary_origin_RI_REV.csv")
write_csv(flagging_by_origin, "data/processed/REVIEW_A4_binary_origin_flagging_REV.csv")
write_csv(dept_map,           "data/processed/REVIEW_A4_binary_origin_deptmap_REV.csv")

out_xlsx <- "data/processed/tables/SUPPLEMENT_binary_origin_REV.xlsx"
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                  halign = "center", border = "bottom")
add_sheet <- function(nm, df) {
  addWorksheet(wb, nm); writeData(wb, nm, df, headerStyle = hs)
  freezePane(wb, nm, firstRow = TRUE)
  setColWidths(wb, nm, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet("Origin_Distribution", origin_dist)
add_sheet("Rule_Distribution",   rule_dist)
add_sheet("RI_Contrast",         ri_contrast)
add_sheet("Flagging_by_Origin",  flagging_by_origin)
add_sheet("Department_Map",      dept_map)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("\n[OK] Saved: %s\n", out_xlsx))
