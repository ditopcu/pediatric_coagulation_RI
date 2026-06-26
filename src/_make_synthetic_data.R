# ==============================================================================
# Generate a small SYNTHETIC dataset for smoke-testing the pipeline.
# Schema mirrors the real Excel files but values are sampled from broad
# clinically plausible distributions. NOT a substitute for the real data:
# numerical results from synthetic-driven pipelines will not match the
# published manuscript. Use only to verify that the code paths run.
# Output: data/synthetic/{tce_2026_1_18_pt, _aptt, _fib}.xlsx
#         data/synthetic/{PT_2, aPTT, Fib_2}_2025_10_17.xlsx (with infant rows)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(writexl)
})

set.seed(2026)

out_dir <- "data/synthetic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

DEPTS <- c("COCUK ACIL",
           "COCUK CERRAHI POL. - 1 (MH 2.KAT)",
           "COCUK GASTROENTROLOJI POL. - 1 (MH 2.KAT)",
           "COCUK KARDIYOLOJI SERVISI (MT1-7C)",
           "COCUK KLINIK SERVIS (MT1-8A)",
           "COCUK ONKOLOJI POL. - 1 (MH 2.KAT)",
           "COCUK SAGLIGI VE HAST POL. -15 (KD 2.KAT)")

# --- Helper: sample one analyte for ages 1..18 with mild age dependence ----
gen_main <- function(n_pt = 1700, n_aptt = 1700, n_fib = 320) {
  base <- function(n) {
    tibble(
      patient_id = sprintf("PID%07d", sample(10000:99999, n, replace = TRUE)),
      sample_id  = sprintf("ISH%08d", sample(10000000:99999999, n, replace = TRUE)),
      department = sample(DEPTS, n, replace = TRUE,
                          prob = c(0.36, 0.06, 0.06, 0.04, 0.04, 0.04, 0.40)),
      sex        = sample(c("K", "E"), n, replace = TRUE, prob = c(0.435, 0.565)),
      age_year   = round(runif(n, 1, 17.99), 2)
    ) |>
      mutate(age_month = floor(age_year * 12),
             age_group = as.character(floor(age_year)))
  }

  # PT around 9.7 (s)
  pt <- base(n_pt) |>
    mutate(test_name = "PT",
           result_num = round(pmax(7.5, rnorm(n(), 9.7 + 0.005 * age_year, 0.55)), 2),
           raw_result = sprintf("%.2f", result_num),
           result_chr = raw_result,
           ek_sonuc_2 = NA_character_)

  # aPTT mean drops slightly with age, partition at 12
  aptt <- base(n_aptt) |>
    mutate(test_name = "aPTT",
           result_num = round(pmax(20, rnorm(n(), 29.8 - 0.10 * age_year, 2.4)), 2),
           raw_result = sprintf("%.2f", result_num),
           result_chr = raw_result,
           ek_sonuc_2 = NA_character_)

  # Fibrinogen
  fib <- base(n_fib) |>
    mutate(test_name = "Fibrinogen",
           result_num = round(pmax(0.8, rnorm(n(), 2.85, 0.42)), 2),
           raw_result = sprintf("%.2f", result_num),
           result_chr = raw_result,
           ek_sonuc_2 = NA_character_)

  list(PT = pt, aPTT = aptt, Fib = fib)
}

# --- Generate legacy file with infants (30-365 days) included --------------
gen_legacy <- function(test_name, n_infant, n_peds, mean_inf, sd_inf,
                       mean_peds, sd_peds) {
  # infant rows: age_year in [30/365, 1)
  infant <- tibble(
    sample_id = sprintf("ISH%08d", sample(10000000:99999999, n_infant, replace = TRUE)),
    sex       = sample(c("K", "E"), n_infant, replace = TRUE, prob = c(0.435, 0.565)),
    age_year  = round(runif(n_infant, 30/365, 0.999), 4),
    age_group = "<1",
    result_num= round(pmax(0.5, rnorm(n_infant, mean_inf, sd_inf)), 2)
  )
  peds <- tibble(
    sample_id = sprintf("ISH%08d", sample(10000000:99999999, n_peds, replace = TRUE)),
    sex       = sample(c("K", "E"), n_peds, replace = TRUE, prob = c(0.435, 0.565)),
    age_year  = round(runif(n_peds, 1, 17.99), 2),
    age_group = as.character(floor(age_year)),
    result_num= round(pmax(0.8, rnorm(n_peds, mean_peds, sd_peds)), 2)
  )
  bind_rows(infant, peds)
}

main <- gen_main()
write_xlsx(main$PT,   file.path(out_dir, "tce_2026_1_18_pt.xlsx"))
write_xlsx(main$aPTT, file.path(out_dir, "tce_2026_1_18_aptt.xlsx"))
write_xlsx(main$Fib,  file.path(out_dir, "tce_2026_1_18_fib.xlsx"))

pt_legacy   <- gen_legacy("PT",         140, 1700, 10.5, 0.7, 9.7, 0.55)
aptt_legacy <- gen_legacy("aPTT",       140, 1700, 32.5, 3.1, 29.5, 2.4)
fib_legacy  <- gen_legacy("Fibrinogen",  22,  320, 2.40, 0.38, 2.85, 0.42)
write_xlsx(pt_legacy,   file.path(out_dir, "PT_2_2025_10_17.xlsx"))
write_xlsx(aptt_legacy, file.path(out_dir, "aPTT_2025_10_17.xlsx"))
write_xlsx(fib_legacy,  file.path(out_dir, "Fib_2_2025_10_17.xlsx"))

cat(sprintf("[OK] Synthetic dataset written under %s\n", out_dir))
cat("Files:\n")
print(list.files(out_dir, full.names = TRUE))
cat("\nNote: synthetic results will NOT match published numbers.\n")
cat("Use only to verify that scripts execute end-to-end.\n")
