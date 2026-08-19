# ==============================================================================
# EXPORT MAIN PUBLICATION TABLES TO EXCEL
# Combines PUB_table_1..6 CSVs into a single .xlsx with six sheets.
# Output: data/processed/tables/PUB_main_tables.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(openxlsx)
})

tbl_dir <- "data/processed/tables"
out_path <- file.path(tbl_dir, "PUB_main_tables.xlsx")

sheets <- list(
  "T1_Study_Population"      = "PUB_table_1_study_population.csv",
  "T2_Reference_Intervals"   = "PUB_table_2_reference_intervals.csv",
  "T3_False_Positive_Rates"  = "PUB_table_3_false_positive_rates.csv",
  "T4_Guven_Transferability" = "PUB_table_4_guven_transferability.csv",
  "T5_International_Comp"    = "PUB_table_5_international_comparison.csv",
  "T6_Harris_Boyd"           = "PUB_table_6_harris_boyd.csv"
)

wb <- createWorkbook()

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9E1F2",
  halign = "center",
  border = "bottom"
)

for (sheet_name in names(sheets)) {
  csv_path <- file.path(tbl_dir, sheets[[sheet_name]])
  df <- read_csv(csv_path, show_col_types = FALSE)

  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df, headerStyle = header_style)
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")

  message(sprintf("[OK] %s: %d rows x %d cols", sheet_name, nrow(df), ncol(df)))
}

saveWorkbook(wb, out_path, overwrite = TRUE)
message(sprintf("\nSaved: %s", out_path))
