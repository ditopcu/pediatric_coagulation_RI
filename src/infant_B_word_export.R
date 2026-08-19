# ==============================================================================
# Build a single Word supplement document containing:
#   - Supplementary Figure S2 (PNG) + caption
#   - Supplementary Table     + caption
# Output: article_update/Supplement_infant_RI.docx
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr)
  library(officer); library(flextable)
})

png_path <- "figures/PNG_300DPI/figure_supplement_infant_RI.png"
csv_path <- "data/processed/tables/PUB_table_S_infant_pediatric_RI.csv"
out_docx <- "article_update/Supplement_infant_RI.docx"

stopifnot(file.exists(png_path), file.exists(csv_path))

tbl <- read_csv(csv_path, show_col_types = FALSE)

# --- Captions ----------------------------------------------------------------
# Opening commentary -- frames the analysis and pre-empts the "why not in main
# manuscript?" reviewer question by explaining the limitations of the infant subset.

intro_para_1 <- paste0(
  "This supplementary material presents an exploratory analysis of coagulation ",
  "reference intervals in infants aged 30 to 365 days, complementing the primary ",
  "1-18-year reference intervals reported in the main manuscript. Although our ",
  "hospital-based dataset contained measurements from this younger age group, ",
  "we did not derive primary infant reference intervals because (i) the available ",
  "infant fibrinogen sample (n = 217) is below the conventional threshold for ",
  "indirect estimation and yields a wide upper-limit confidence interval; ",
  "(ii) coagulation undergoes rapid postnatal maturation that would not be ",
  "adequately captured by a single 30-365 d interval, while partitioning into ",
  "developmentally homogeneous sub-bins (e.g., 1-3, 4-6, 7-12 months) would have ",
  "yielded individual sub-bins too small for reliable indirect estimation; ",
  "(iii) the infant patient distribution is heavily skewed toward symptomatic ",
  "presentations (predominantly pediatric emergency and gastroenterology services), ",
  "which weakens the healthy-referral assumption that underlies refineR more ",
  "severely than in older age groups; and (iv) several aPTT measurements within ",
  "the infant subset showed extreme values (minimum 1.37 s, maximum 270 s) ",
  "consistent with preanalytical contamination or critical pathology. We therefore ",
  "present the analysis below as a supportive, hypothesis-generating observation ",
  "rather than as a definitive infant reference interval."
)

intro_para_2 <- paste0(
  "Despite these caveats, the infant cohort consistently showed longer PT and aPTT ",
  "upper reference limits and lower fibrinogen limits than the pediatric cohort ",
  "(Figure S2, Table below). Applying the pediatric 1-18 y refineR reference interval ",
  "to the infant subset, 19.5% of PT, 21.1% of aPTT, and 14.3% of fibrinogen results ",
  "fell outside the pediatric reference range - well above the 5% expected under a ",
  "correctly transferred interval. The direction of the deviations is consistent ",
  "with developmental hemostasis: vitamin-K-dependent and contact-pathway factor ",
  "activities rise toward adult values over the first year of life, while fibrinogen ",
  "increases somewhat more slowly. These observations underscore that the pediatric ",
  "refineR reference intervals derived in children aged 1 year and older should not ",
  "be applied to infants under one year, and that age-specific infant reference ",
  "intervals derived from a dedicated healthy infant cohort remain warranted."
)

# Short figure caption (sample sizes/numbers live in the Table -- not repeated here)
fig_title <- "Supplementary Figure S2."

fig_caption <- paste0(
  "refineR reference interval estimates for PT (A), aPTT (B), and fibrinogen (C) ",
  "in infants (30-365 days, blue) and the pediatric cohort (1-18 y, orange). ",
  "For aPTT, the pediatric cohort is shown by the manuscript's Harris-Boyd partitions ",
  "(1-12 y and 12-18 y). Vertical bars span the 2.5th-97.5th percentile interval; ",
  "whiskers represent 90% bootstrap confidence intervals at each reference limit ",
  "(NBootstrap = 200). Sample sizes and numerical values are provided in the ",
  "Supplementary Table below."
)

tbl_title <- "Supplementary Table."

tbl_caption <- paste0(
  "Two-sided refineR reference intervals for PT, aPTT, and fibrinogen, separately for ",
  "the infant (30-365 days) and pediatric (1-18 y) cohorts. Lower (2.5th percentile) ",
  "and upper (97.5th percentile) limits are reported with 90% bootstrap confidence ",
  "intervals in parentheses (NBootstrap = 200; tolerance intervals). For aPTT, the ",
  "pediatric cohort is partitioned at 12 years per the manuscript's Harris-Boyd analysis; ",
  "PT and fibrinogen are pooled across 1-18 y. No outlier filter was applied to the ",
  "infant data."
)

# --- Build flextable ---------------------------------------------------------
ft <- flextable(tbl) |>
  bold(part = "header") |>
  bg(part = "header", bg = "#D9E1F2") |>
  align(align = "center", part = "all") |>
  align(j = 1, align = "left", part = "all") |>
  border_outer(border = fp_border(color = "#444444", width = 0.75)) |>
  border_inner_h(border = fp_border(color = "#BBBBBB", width = 0.4)) |>
  fontsize(size = 10, part = "all") |>
  set_table_properties(layout = "autofit", width = 1)

# --- Compose Word document --------------------------------------------------
doc <- read_docx() |>
  body_add_par("Supplementary Materials - Infant (30-365 days) Reference Intervals",
               style = "heading 1") |>
  body_add_par("", style = "Normal") |>
  body_add_fpar(fpar(ftext(intro_para_1, prop = fp_text(font.size = 11)))) |>
  body_add_par("", style = "Normal") |>
  body_add_fpar(fpar(ftext(intro_para_2, prop = fp_text(font.size = 11)))) |>
  body_add_par("", style = "Normal") |>
  body_add_img(src = png_path, width = 6.5, height = 2.4) |>
  body_add_fpar(fpar(
    ftext(fig_title, prop = fp_text(bold = TRUE, font.size = 10)),
    ftext(" ", prop = fp_text(font.size = 10)),
    ftext(fig_caption, prop = fp_text(font.size = 10))
  )) |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft) |>
  body_add_fpar(fpar(
    ftext(tbl_title, prop = fp_text(bold = TRUE, font.size = 10)),
    ftext(" ", prop = fp_text(font.size = 10)),
    ftext(tbl_caption, prop = fp_text(font.size = 10))
  ))

dir.create(dirname(out_docx), recursive = TRUE, showWarnings = FALSE)
print(doc, target = out_docx)
cat(sprintf("\n[OK] Word document saved: %s\n", out_docx))
