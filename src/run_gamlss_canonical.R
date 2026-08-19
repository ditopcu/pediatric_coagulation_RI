# ==============================================================================
# GAMLSS PIPELINE RUN FROM THE CANONICAL FILES (NBootstrap = 100)
# ==============================================================================
# Runs the Ammer pipeline end to end -- point estimate AND bootstrap -- for PT,
# aPTT and fibrinogen, using the inputs built by src/build_gamlss_inputs.R from
# the three canonical xlsx files.
#
# Expects cont_out/ to be empty of PT/aPTT/Fib results, so both the point
# estimates and all bootstrap replicates are computed fresh. Mixing replicate
# vintages produces confidence intervals drawn from different runs, because
# common.R scans indices 1:200 behind a file.exists guard.
#
# Run from the project root:  Rscript src/run_gamlss_canonical.R
# Log: data/processed/gamlss_canonical_run_log.txt
# ==============================================================================

library(refineR)
library(gamlss)
library(parallel)
library(future)
library(future.apply)

options("future.globals.maxSize" = 1500 * 1024^2)

# NOTE: src/pipeline/mode_deneme.R is deliberately NOT sourced. It redefines
# defineAgeGroupsWithTol with a debugging version that hard-codes covarName,
# colnameID and colnameValue inside the body, reads a global `inputData`, and
# carries a leftover print() call. The original runner scripts source it last,
# so they pick that version up; algoRICurves.R holds the real one.
source("src/pipeline/utils.R")
source("src/pipeline/algoRICurves.R")

# --- Workaround for a latent free-variable bug in the protected pipeline ---
# Both defineAgeGroupsWithTol() (point estimate) and defineAgeGroupsWithTolBS()
# (bootstrap) build an initial age group (0,0). When the input ages start at 1
# that group is empty, so the local `ids` is never assigned, and the expansion
# step then reads `ids` before assignment in removeDuplIDs(). Under future
# multisession the variable cannot be resolved in the workers ("object 'ids'
# not found"). A global `ids` is NOT exported by future's globals scanner (it
# classifies `ids` as local because the function also assigns it), so the fix
# must live inside the function body: inject a local `ids <- character(0)` at
# the top of the function IN MEMORY. src/pipeline/ is not modified.
# With an empty `ids`, the duplicate-ID screen is a no-op for the empty group,
# which is the behaviourally neutral choice.
patch_ids_init <- function(fn) {
  b <- as.list(body(fn))
  body(fn) <- as.call(append(b, list(quote(ids <- character(0))), after = 1))
  fn
}
defineAgeGroupsWithTol   <- patch_ids_init(defineAgeGroupsWithTol)
defineAgeGroupsWithTolBS <- patch_ids_init(defineAgeGroupsWithTolBS)

outputDir  <- "cont_out/"
NBootstrap <- 100
log_file   <- "data/processed/gamlss_canonical_run_log.txt"

log_msg <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
  message(txt)
  cat(txt, "\n", file = log_file, append = TRUE)
}

log_msg("=== GAMLSS canonical run START ===")
log_msg("[PATCH] local ids initialization injected into ",
        "defineAgeGroupsWithTol and defineAgeGroupsWithTolBS.")
log_msg("[NOTE] mode_deneme.R not sourced (it overrides ",
        "defineAgeGroupsWithTol with a debugging version).")
log_msg("R: ", R.version.string,
        " | cores: ", detectCores(logical = TRUE),
        " | NBootstrap: ", NBootstrap)

analytes <- list(
  list(csv = "data/coa_results/gamlss_input/gamlss_input_PT.csv",
       filename = "PT",   expect_n = 17000),
  list(csv = "data/coa_results/gamlss_input/gamlss_input_aPTT.csv",
       filename = "aPTT", expect_n = 17001),
  list(csv = "data/coa_results/gamlss_input/gamlss_input_Fib.csv",
       filename = "Fib",  expect_n = 3166)
)

# Guard: refuse to start if results from an earlier run are still present.
stale <- list.files(outputDir, pattern = "^(PT|aPTT|Fib)_.*\\.RData$")
if (length(stale) > 0) {
  log_msg("FATAL: cont_out/ still holds ", length(stale),
          " PT/aPTT/Fib result files. Archive them first so point estimates ",
          "and bootstrap replicates come from one run. ABORTING.")
  quit(status = 1)
}
log_msg("cont_out/ is clear of PT/aPTT/Fib results -- starting fresh.")

for (a in analytes) tryCatch({

  log_msg("[", a$filename, "] reading ", a$csv)
  inputData <- read.table(file = a$csv, header = TRUE, sep = ";",
                          stringsAsFactors = FALSE)
  log_msg("[", a$filename, "] input rows: ", nrow(inputData),
          " (expected ", a$expect_n, ")")
  if (nrow(inputData) != a$expect_n) {
    stop("row count does not match the canonical file")
  }

  seed <- 123
  set.seed(seed)

  t0 <- Sys.time()
  log_msg("[", a$filename, "] point estimate START")
  pointEstGamlss <- runPipelinePointEst(inputData = inputData,
                                        covarName = "Age",
                                        colnameID = "PID",
                                        colnameValue = "Value",
                                        outputDir = outputDir,
                                        filename = a$filename,
                                        NCores = NULL)
  log_msg("[", a$filename, "] point estimate DONE in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

  t1 <- Sys.time()
  log_msg("[", a$filename, "] bootstrap START (", NBootstrap, " iterations)")
  bootstrapRes <- runPipelineCIs(inputData = inputData,
                                 pointEstGamlss = pointEstGamlss,
                                 NBootstrap = NBootstrap,
                                 covarName = "Age",
                                 colnameID = "PID",
                                 colnameValue = "Value",
                                 outputDir = outputDir,
                                 filename = a$filename,
                                 NCores = NULL)

  n_ok <- sum(vapply(bootstrapRes, function(x) !is.null(x$Models), logical(1)))
  log_msg("[", a$filename, "] bootstrap DONE: ", length(bootstrapRes),
          " iterations (", n_ok, " with a successful GAMLSS fit) in ",
          round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), " min")
  log_msg("[", a$filename, "] TOTAL ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

}, error = function(e) {
  log_msg("[", a$filename, "] ERROR: ", conditionMessage(e),
          " -- continuing with the next analyte.")
})

# Final inventory so a partial run is obvious rather than silent.
log_msg("--- replicate inventory ---")
for (nm in c("PT", "aPTT", "Fib")) {
  n <- length(list.files(outputDir,
                         pattern = paste0("^", nm, "_gamlssModel_Est_[0-9]+\\.RData$")))
  log_msg("  ", nm, ": ", n, " bootstrap replicate files")
}
log_msg("=== GAMLSS canonical run END ===")
