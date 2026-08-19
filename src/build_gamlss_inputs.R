# ==============================================================================
# BUILD GAMLSS PIPELINE INPUTS FROM THE CANONICAL FILES
# ==============================================================================
# Writes one semicolon-separated file per analyte in the format the Ammer
# pipeline expects: PID; Age; Value.
#
#   PID   = patient_id   (the pipeline uses this to recognise samples from the
#                         same subject; ProbNPFinal = ProbNP / Freq)
#   Age   = floor(age_year)   integer years, per CLAUDE.md
#   Value = result_num
#
# Run from the project root:  Rscript src/build_gamlss_inputs.R
# ==============================================================================

library(readxl)
library(dplyr)

OUT <- "data/coa_results/gamlss_input"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

analytes <- tibble::tribble(
  ~name,  ~src,                                       ~expect_n,
  "PT",   "data/coa_results/tce 2026 1-18 pt.xlsx",    17000,
  "aPTT", "data/coa_results/tce 2026 1-18 aptt.xlsx",  17001,
  "Fib",  "data/coa_results/tce 2026 1-18 fib.xlsx",    3166
)

ok <- TRUE

for (i in seq_len(nrow(analytes))) {
  a <- analytes[i, ]
  d <- read_excel(a$src) |> janitor::clean_names()

  stopifnot(all(c("patient_id", "age_year", "result_num") %in% names(d)))

  out <- d |>
    transmute(PID = patient_id,
              Age = as.integer(floor(age_year)),
              Value = result_num)

  # Checks -- abort loudly rather than write a silently wrong input.
  n_na <- sum(is.na(out$PID) | is.na(out$Age) | is.na(out$Value))
  n_ok <- nrow(out) == a$expect_n

  path <- file.path(OUT, sprintf("gamlss_input_%s.csv", a$name))
  write.table(out, file = path, sep = ";", row.names = FALSE,
              quote = FALSE, fileEncoding = "UTF-8")

  cat(sprintf("\n=== %s ===\n", a$name))
  cat(sprintf("  source        : %s\n", a$src))
  cat(sprintf("  rows written  : %d (expected %d) %s\n",
              nrow(out), a$expect_n, if (n_ok) "OK" else "MISMATCH"))
  cat(sprintf("  NA in any col : %d %s\n", n_na, if (n_na == 0) "OK" else "PROBLEM"))
  cat(sprintf("  Age range     : %d - %d (integer)\n", min(out$Age), max(out$Age)))
  cat(sprintf("  Value range   : %.3f - %.3f\n", min(out$Value), max(out$Value)))
  cat(sprintf("  unique PID    : %d\n", dplyr::n_distinct(out$PID)))
  cat(sprintf("  samples/subject: max %d\n",
              max(table(out$PID))))
  cat(sprintf("  written to    : %s\n", path))

  if (!n_ok || n_na > 0) ok <- FALSE
}

cat("\n==============================================\n")
cat(sprintf("ALL INPUTS VALID: %s\n", ok))
cat("==============================================\n")
if (!ok) quit(status = 1)
