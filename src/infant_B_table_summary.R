# ==============================================================================
# SUPPLEMENT TABLE -- Reference intervals with sample size and 90% CI
# Builds a clean publication-ready table combining infant (30-365 d) and
# pediatric (1-18 y / 1-12 y / 12-18 y) refineR reference intervals.
# Output:
#   data/processed/tables/PUB_table_S_infant_pediatric_RI.csv
#   appended sheet "RI_Table_Summary" in SUPPLEMENT_infant_30_365d.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(openxlsx)
})

if (!exists("PAL")) source("src/common.R")

# --- Pediatric refineR with CI (from PUB_table_2) --------------------------
peds <- read_csv("data/processed/tables/PUB_table_2_reference_intervals.csv",
                 show_col_types = FALSE) |>
  filter(method == "refineR") |>
  select(test, age_group, n, lower, lower_ci_lo, lower_ci_hi,
         upper, upper_ci_lo, upper_ci_hi) |>
  mutate(group = "Pediatric")

# --- Infant refineR with CI ------------------------------------------------
infant <- read_csv("data/processed/infant_B_refineR_RI.csv",
                   show_col_types = FALSE) |>
  mutate(age_group = "30-365 d", group = "Infant") |>
  select(test, age_group, n, lower, lower_ci_lo, lower_ci_hi,
         upper, upper_ci_lo, upper_ci_hi, group)

# --- Combine and order -----------------------------------------------------
fmt <- function(x, ci_lo, ci_hi, d = 2) {
  sprintf("%.*f (%.*f-%.*f)", d, x, d, ci_lo, d, ci_hi)
}

combined <- bind_rows(infant, peds) |>
  mutate(test = factor(test, levels = c("PT", "aPTT", "Fibrinogen")),
         age_order = case_when(
           age_group == "30-365 d" ~ 1,
           age_group == "1-12 y"   ~ 2,
           age_group == "12-18 y"  ~ 3,
           age_group == "1-18 y"   ~ 4,
           TRUE                    ~ 5
         )) |>
  arrange(test, age_order) |>
  mutate(
    lower_str = fmt(lower, lower_ci_lo, lower_ci_hi),
    upper_str = fmt(upper, upper_ci_lo, upper_ci_hi),
    test_label = case_when(
      test == "PT"         ~ "PT, s",
      test == "aPTT"       ~ "aPTT, s",
      test == "Fibrinogen" ~ "Fibrinogen, g/L"
    )
  ) |>
  select(`Test` = test_label,
         `Group` = group,
         `Age group` = age_group,
         `n` = n,
         `Lower limit, 2.5th (90% CI)` = lower_str,
         `Upper limit, 97.5th (90% CI)` = upper_str)

cat("\n=== SUPPLEMENT TABLE: Infant + Pediatric refineR RI ===\n")
print(combined)

# --- Save CSV --------------------------------------------------------------
csv_path <- "data/processed/tables/PUB_table_S_infant_pediatric_RI.csv"
write_csv(combined, csv_path)
cat(sprintf("\n[OK] CSV: %s\n", csv_path))

# --- Append "RI_Table_Summary" sheet to existing supplement Excel ----------
xlsx_path <- "data/processed/tables/SUPPLEMENT_infant_30_365d.xlsx"
wb <- loadWorkbook(xlsx_path)

# Remove sheet if it exists (so we cleanly replace)
if ("RI_Table_Summary" %in% names(wb)) removeWorksheet(wb, "RI_Table_Summary")

addWorksheet(wb, "RI_Table_Summary")
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                  halign = "center", border = "bottom")
writeData(wb, "RI_Table_Summary", combined, headerStyle = hs)
freezePane(wb, "RI_Table_Summary", firstRow = TRUE)
setColWidths(wb, "RI_Table_Summary",
             cols = seq_len(ncol(combined)), widths = "auto")

saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(sprintf("[OK] Excel sheet appended: %s -> 'RI_Table_Summary'\n", xlsx_path))
