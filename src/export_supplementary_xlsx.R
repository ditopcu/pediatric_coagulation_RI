# ==============================================================================
# EXPORT SUPPLEMENTARY TABLES TO EXCEL
# Combines PUB_table_S1, S2, S3 CSVs into a single .xlsx with three sheets.
# Output: data/processed/tables/PUB_supplementary_tables.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(openxlsx)
})

tbl_dir <- "data/processed/tables"
out_path <- file.path(tbl_dir, "PUB_supplementary_tables.xlsx")

# Sheet numbering matches manuscript references:
#   S1 = sample sizes + descriptive statistics
#   S2 = analytical performance (within-run CV%)
#   S3 = sex-specific refineR reference limits
#   S4 = continuous GAMLSS reference limit values
sheets <- list(
  "S1_Descriptive"           = "PUB_table_S1_descriptive_all.csv",
  "S2_Analytical_Performance"= "TEZ_tablo_4_5_analytical_performance.csv",
  "S3_Sex_refineR_RI"        = "PUB_table_S2_refineR_all_by_age_sex.csv",
  "S4_Continuous_RI"         = "PUB_table_S3_continuous_RI_values.csv"
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
