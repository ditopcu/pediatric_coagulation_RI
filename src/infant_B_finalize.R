# ==============================================================================
# INFANT ANALYSIS - STEP B (FINALIZE)
# Resumes from saved RDS files (refineR fits) and produces:
#   - Comparison figure (pediatric vs infant RI bands)
#   - Supplement Excel (Descriptive, Classification, Infant_RI, Pediatric_vs_Infant)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(patchwork); library(openxlsx)
  library(refineR)
})

if (!exists("PAL")) source("src/common.R")

# --- Load saved infant RI CSV (already extracted) ---------------------------
infant_ri_disp <- read_csv("data/processed/infant_B_refineR_RI.csv",
                           show_col_types = FALSE)
compare_tbl    <- read_csv("data/processed/infant_B_pediatric_vs_infant.csv",
                           show_col_types = FALSE)

cat("\n=== Infant refineR RI (loaded) ===\n"); print(infant_ri_disp)
cat("\n=== Pediatric vs Infant (loaded) ===\n"); print(compare_tbl)

# --- Build long-format data for figure -------------------------------------
peds_rows <- compare_tbl |>
  transmute(test, group = "Pediatric (current study)",
            lower = peds_lower, upper = peds_upper,
            lower_ci_lo = NA_real_, lower_ci_hi = NA_real_,
            upper_ci_lo = NA_real_, upper_ci_hi = NA_real_)

infant_rows <- compare_tbl |>
  transmute(test, group = "Infant (30-365 days)",
            lower = infant_lower, upper = infant_upper,
            lower_ci_lo = infant_lower_ci_lo, lower_ci_hi = infant_lower_ci_hi,
            upper_ci_lo = infant_upper_ci_lo, upper_ci_hi = infant_upper_ci_hi)

fig_data <- bind_rows(peds_rows, infant_rows) |>
  mutate(group = factor(group, levels = c("Pediatric (current study)",
                                          "Infant (30-365 days)")),
         test  = factor(test,  levels = c("PT", "aPTT", "Fibrinogen")))

make_panel <- function(test_name, y_label, tag_label) {
  df <- fig_data |> filter(test == test_name)
  ggplot(df, aes(x = group, colour = group)) +
    geom_linerange(aes(ymin = lower, ymax = upper),
                   linewidth = 6, alpha = 0.55) +
    geom_point(aes(y = lower), size = 2.6) +
    geom_point(aes(y = upper), size = 2.6) +
    geom_errorbar(aes(ymin = lower_ci_lo, ymax = lower_ci_hi),
                  width = 0.18, linewidth = 0.6, na.rm = TRUE) +
    geom_errorbar(aes(ymin = upper_ci_lo, ymax = upper_ci_hi),
                  width = 0.18, linewidth = 0.6, na.rm = TRUE) +
    scale_colour_manual(values = c(
      "Pediatric (current study)" = PAL$refiner,
      "Infant (30-365 days)"      = PAL$accent3
    )) +
    labs(x = NULL, y = y_label, tag = tag_label) +
    theme_tufte_academic() +
    theme(plot.tag = element_text(size = 14, face = "bold"),
          legend.position = "none",
          axis.text.x = element_text(size = 8, angle = 12, hjust = 0.8))
}

p_pt   <- make_panel("PT",         "PT, seconds",     "A")
p_aptt <- make_panel("aPTT",       "aPTT, seconds",   "B")
p_fib  <- make_panel("Fibrinogen", "Fibrinogen, g/L", "C")

fig <- (p_pt | p_aptt | p_fib) +
  plot_annotation(theme = theme_tufte_academic())

save_fig(fig, "figure_supplement_infant_RI", width = 220, height = 100)
cat("[OK] Figure saved.\n")

# --- Supplement Excel ------------------------------------------------------
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

desc      <- read_csv("data/processed/infant_A_descriptive.csv",      show_col_types = FALSE)
class_tbl <- read_csv("data/processed/infant_C_classification.csv",  show_col_types = FALSE)

add_sheet("Descriptive",          desc)
add_sheet("Classification_vs_RI", class_tbl)
add_sheet("Infant_refineR_RI",    infant_ri_disp)
add_sheet("Pediatric_vs_Infant",  compare_tbl)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("[OK] Supplement Excel saved: %s\n", out_xlsx))

cat("\n=========================================================\n")
cat("[DONE] Infant analysis (steps A + B + C) complete.\n")
cat("=========================================================\n")
