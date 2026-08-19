# ==============================================================================
# CONFIDENCE INTERVALS FOR SUPPLEMENTARY TABLES S6 AND S7
# ==============================================================================
# The two reviewer-facing sensitivity tables were first produced with refineR
# point estimates only (NBootstrap = 1). This script recomputes every cell with
# NBootstrap = 200, the same setting as the published reference intervals, so the
# tables carry 90% confidence intervals on each limit. It also adds 95% bootstrap
# confidence intervals on the out-of-RI flagging rates, by the same procedure as
# Table 3 (B = 1000, percentile method, seed 123).
#
# Subset definitions are copied from the two scripts that produced the point
# estimates, and the resulting subset sizes are asserted against the values those
# scripts reported, so a silent divergence fails loudly:
#   src/analysis4_binary_origin_REV.R          (S6, patient origin)
#   src/analysis_crosstest_criteria_REV.R      (S7, cross-test exclusion)
#
# Fits run sequentially; refineR parallelises internally. 16 fits, ~7 min each.
#
# Outputs (new files only):
#   data/processed/REVIEW_S6_patient_origin_CI.csv
#   data/processed/REVIEW_S7_crosstest_CI.csv
#   data/processed/REVIEW_S6_S7_flagging_CI.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(janitor); library(readr); library(tidyr)
  library(tibble); library(refineR)
})

NBOOT_RI   <- 200    # matches the published reference intervals
NBOOT_FLAG <- 1000   # matches Table 3


# Reference limits come from src/common.R, which is the single source of truth
# for `ri_refiner`, `ri_manufacturer` and `ri_direct`. Sourcing it also loads
# the analysis data and the GAMLSS continuous curves from cont_out/, so that
# directory has to hold a complete run before this script will start.
source("src/common.R")

# Expand an age-partitioned limit table into one row per integer age. The
# partition labels follow common.R: "1-12" is age_year < 12, so age_int 1-11.
expand_by_age <- function(ri, ages = 1:18) {
  ri |>
    rowwise() |>
    mutate(age_int = list(switch(age_group,
                                 "1-18"  = ages,
                                 "1-12"  = ages[ages <  12],
                                 "12-18" = ages[ages >= 12]))) |>
    ungroup() |>
    unnest(age_int) |>
    select(test, age_int, lower, upper)
}

main_files <- c(
  PT         = "data/coa_results/tce 2026 1-18 pt.xlsx",
  aPTT       = "data/coa_results/tce 2026 1-18 aptt.xlsx",
  Fibrinogen = "data/coa_results/tce 2026 1-18 fib.xlsx"
)
load_main <- function(path, tname) {
  suppressMessages(read_excel(path)) |>
    clean_names() |>
    mutate(test = tname, age_int = as.integer(floor(age_year)))
}
raw <- bind_rows(
  load_main(main_files["PT"],         "PT"),
  load_main(main_files["aPTT"],       "aPTT"),
  load_main(main_files["Fibrinogen"], "Fibrinogen")
)

# --- S6: binary patient origin (rules from analysis4_binary_origin_REV.R) ----
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
raw <- raw |> mutate(origin = classify_origin(department))

# --- S7: criterion-based cross-test exclusion (from analysis_crosstest_criteria_REV.R)
crit_manufacturer <- ri_manufacturer |>
  crossing(age_int = 1:18) |>
  select(test, age_int, lower, upper)
crit_direct <- expand_by_age(ri_direct)

keep_after_crosstest <- function(dat, crit, target) {
  ab <- dat |>
    left_join(crit, by = c("test", "age_int")) |>
    mutate(abn = result_num < lower | result_num > upper) |>
    select(sample_id, test, abn)
  others <- setdiff(c("PT", "aPTT", "Fibrinogen"), target)
  d_t <- dat |> filter(test == target)
  flags <- d_t |> select(sample_id)
  for (o in others) {
    f_o <- ab |> filter(test == o) |> select(sample_id, !!o := abn)
    flags <- flags |> left_join(f_o, by = "sample_id")
  }
  flags <- flags |>
    mutate(across(all_of(others), ~ coalesce(.x, FALSE))) |>
    mutate(drop_any = if_any(all_of(others), identity))
  d_t |> filter(!flags$drop_any)
}

# --- Assemble every subset that needs a fit ----------------------------------
part_rows <- function(d, target) {
  if (target == "aPTT") {
    list(list(age_group = "1-12 y",  x = d$result_num[d$age_int <  12]),
         list(age_group = "12-18 y", x = d$result_num[d$age_int >= 12]))
  } else {
    list(list(age_group = "1-18 y", x = d$result_num))
  }
}

jobs <- list()
add_jobs <- function(jobs, d, target, table_id, variant) {
  for (pr in part_rows(d, target)) {
    jobs[[length(jobs) + 1]] <- list(table = table_id, test = target,
                                     variant = variant, age_group = pr$age_group,
                                     x = pr$x)
  }
  jobs
}

# The full-cohort rows of both tables are the published reference intervals
# themselves -- same data, same NBootstrap -- so they are taken from
# ri_comparison rather than refitted. Only the subsets are fitted here.
for (t in c("PT", "aPTT", "Fibrinogen")) {
  d_all <- raw |> filter(test == t)
  jobs <- add_jobs(jobs, d_all |> filter(origin == "Outpatient"), t, "S6", "Outpatient")
  jobs <- add_jobs(jobs, d_all |> filter(origin == "Inpatient"),  t, "S6", "Inpatient")

  jobs <- add_jobs(jobs, keep_after_crosstest(raw, crit_manufacturer, t), t,
                   "S7", "Manufacturer RI")
  jobs <- add_jobs(jobs, keep_after_crosstest(raw, crit_direct, t), t,
                   "S7", "Direct RI")
}

cat("\n=== subset sizes ===\n")
sizes <- tibble(table = vapply(jobs, `[[`, "", "table"),
                test = vapply(jobs, `[[`, "", "test"),
                variant = vapply(jobs, `[[`, "", "variant"),
                age_group = vapply(jobs, `[[`, "", "age_group"),
                n = vapply(jobs, function(j) length(j$x), integer(1)))
print(sizes, n = Inf)

# Assert against the sizes the point-estimate scripts reported.
expect <- tribble(
  ~table, ~test, ~variant, ~age_group, ~n,
  "S6", "PT", "Outpatient", "1-18 y", 15127L,
  "S6", "PT", "Inpatient",  "1-18 y",  1870L,
  "S6", "aPTT", "Outpatient", "1-12 y", 8975L,
  "S6", "aPTT", "Inpatient",  "1-12 y", 1113L,
  "S6", "Fibrinogen", "Outpatient", "1-18 y", 2618L,
  "S6", "Fibrinogen", "Inpatient",  "1-18 y",  548L,
  "S7", "PT", "Manufacturer RI", "1-18 y", 11777L,
  "S7", "PT", "Direct RI",       "1-18 y", 15169L,
  "S7", "Fibrinogen", "Manufacturer RI", "1-18 y", 2134L,
  "S7", "Fibrinogen", "Direct RI",       "1-18 y", 2691L
)
chk <- expect |> left_join(sizes, by = c("table", "test", "variant", "age_group"),
                           suffix = c("_expected", "_actual"))
if (any(chk$n_expected != chk$n_actual)) {
  print(chk)
  stop("subset sizes do not match the point-estimate scripts -- definitions have diverged")
}
cat("[OK] all asserted subset sizes match the point-estimate scripts.\n")

# --- Fit ---------------------------------------------------------------------
# refineR parallelises internally over all available cores. Wrapping it in
# future workers makes it try to start its own cluster inside a worker that is
# limited to one core, which it refuses to do. So the fits run sequentially and
# refineR is left to use the machine itself.
LOG <- "data/processed/supp_S6_S7_CI_log.txt"
say <- function(...) {
  txt <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", ...)
  cat(txt, "\n"); cat(txt, "\n", file = LOG, append = TRUE)
}
cat("", file = LOG)
say(sprintf("[FIT] %d subsets, NBootstrap = %d, sequential, refineR uses %d cores",
            length(jobs), NBOOT_RI, parallel::detectCores(logical = TRUE)))

res <- lapply(seq_along(jobs), function(i) {
  j <- jobs[[i]]
  say(sprintf("  [%2d/%2d] %s %s %s %s  n=%d ...",
              i, length(jobs), j$table, j$test, j$variant, j$age_group, length(j$x)))
  t0 <- Sys.time()
  fit <- tryCatch(refineR::findRI(j$x, NBootstrap = NBOOT_RI),
                  error = function(e) e)
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (inherits(fit, "error")) {
    return(data.frame(table = j$table, test = j$test, variant = j$variant,
                      age_group = j$age_group, n = length(j$x),
                      lower = NA_real_, lower_ci_lo = NA_real_, lower_ci_hi = NA_real_,
                      upper = NA_real_, upper_ci_lo = NA_real_, upper_ci_hi = NA_real_,
                      secs = secs,
                      note = paste("ERROR:", conditionMessage(fit))))
  }
  r <- as.data.frame(refineR::getRI(fit, RIperc = c(0.025, 0.975), CIprop = 0.90))
  lo <- r[r$Percentile == 0.025, ]; hi <- r[r$Percentile == 0.975, ]
  data.frame(table = j$table, test = j$test, variant = j$variant,
             age_group = j$age_group, n = length(j$x),
             lower = round(lo$PointEst, 2),
             lower_ci_lo = round(lo$CILow, 2), lower_ci_hi = round(lo$CIHigh, 2),
             upper = round(hi$PointEst, 2),
             upper_ci_lo = round(hi$CILow, 2), upper_ci_hi = round(hi$CIHigh, 2),
             secs = secs,
             note = "")
})

ri_out <- as_tibble(bind_rows(res))
say("[FIT] done")

nfail <- sum(grepl("^ERROR", ri_out$note))
say(sprintf("[FIT] failures: %d of %d;  total fit time %.1f min",
            nfail, nrow(ri_out), sum(ri_out$secs) / 60))
if (nfail > 0) {
  cat("\n=== FAILURES ===\n")
  print(as.data.frame(ri_out[grepl("^ERROR", ri_out$note),
                             c("table", "test", "variant", "age_group", "n", "note")]))
  stop("refineR failed on ", nfail, " subsets -- not writing partial output")
}

cat("\n=== REFERENCE INTERVALS WITH 90% CI ===\n")
print(as.data.frame(ri_out))

write_csv(ri_out |> filter(table == "S6"),
          "data/processed/REVIEW_S6_patient_origin_CI.csv")
write_csv(ri_out |> filter(table == "S7"),
          "data/processed/REVIEW_S7_crosstest_CI.csv")

# --- Flagging rates with 95% bootstrap CI ------------------------------------
lut <- bind_rows(
  ri_refiner |> filter(test == "PT")                         |> crossing(age_int = 1:18),
  ri_refiner |> filter(test == "Fibrinogen")                 |> crossing(age_int = 1:18),
  ri_refiner |> filter(test == "aPTT", age_group == "1-12")  |> crossing(age_int = 1:11),
  ri_refiner |> filter(test == "aPTT", age_group == "12-18") |> crossing(age_int = 12:18)
) |> select(test, age_int, pub_lower = lower, pub_upper = upper)

boot_ci <- function(ind, B = NBOOT_FLAG) {
  n <- length(ind)
  if (n == 0) return(c(NA, NA, NA))
  est <- mean(ind) * 100
  bs <- replicate(B, mean(sample(ind, n, replace = TRUE)) * 100)
  c(est, unname(quantile(bs, 0.025)), unname(quantile(bs, 0.975)))
}

flag_base <- raw |>
  left_join(lut, by = c("test", "age_int")) |>
  mutate(grp = case_when(
    test == "PT"                   ~ "PT 1-18 y",
    test == "aPTT" & age_int <  12 ~ "aPTT 1-12 y",
    test == "aPTT" & age_int >= 12 ~ "aPTT 12-18 y",
    test == "Fibrinogen"           ~ "Fibrinogen 1-18 y"),
    outside = result_num < pub_lower | result_num > pub_upper)

set.seed(123)
fl <- list()
for (g in unique(flag_base$grp)) {
  sub <- flag_base |> filter(grp == g)
  for (o in c("All", "Outpatient", "Inpatient")) {
    ind <- if (o == "All") sub$outside else sub$outside[sub$origin == o]
    v <- boot_ci(ind)
    fl[[length(fl) + 1]] <- tibble(grp = g, origin = o, n = length(ind),
                                   flag_pct = round(v[1], 2),
                                   ci_lo = round(v[2], 2), ci_hi = round(v[3], 2))
  }
}
flag_out <- bind_rows(fl)
cat("\n=== FLAGGING AT PUBLISHED RI, BY ORIGIN, WITH 95% CI ===\n")
print(as.data.frame(flag_out))
write_csv(flag_out, "data/processed/REVIEW_S6_S7_flagging_CI.csv")

cat("\n[OK] All outputs written.\n")
