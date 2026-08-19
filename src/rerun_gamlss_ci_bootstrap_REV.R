# ==============================================================================
# RERUN GAMLSS BOOTSTRAP CIs WITH NBootstrap = 200 (reviewer revision, P1)
# ==============================================================================
# The published continuous RI 90% CI bands were based on NBootstrap = 5
# (indefensible). This driver re-runs the bootstrap stage of the Ammer pipeline
# with NBootstrap = 200 for PT, aPTT, and Fibrinogen.
#
# Key properties:
# - Point estimates are NOT recomputed: runPipelinePointEst() loads the cached
#   <Test>_gamlss_PointEst.RData from cont_out/ (pipeline caching).
# - Bootstrap iterations are ADDITIVE: existing <Test>_gamlssModel_Est_1..5
#   files are reused; only the missing iterations 6..200 are computed and
#   saved as new files in cont_out/. No existing file is modified or deleted.
# - Input files are exactly those referenced by the per-analyte runner scripts
#   (src/RUNNER_Pipeline {pt,aptt,fib}.R). An input-consistency check against
#   the data stored inside the cached point-estimate object aborts the run on
#   mismatch instead of bootstrapping from the wrong input.
#
# Run from the project root:  Rscript src/rerun_gamlss_ci_bootstrap_REV.R
# Timing log: data/processed/REVIEW_gamlss_rerun_log_REV.txt
# ==============================================================================

library(refineR)
library(gamlss)
library(parallel)
library(future)
library(future.apply)

options("future.globals.maxSize" = 1500 * 1024^2)

source("src/pipeline/utils.R")
source("src/pipeline/algoRICurves.R")

# --- Workaround for a latent free-variable bug in the protected pipeline ---
# defineAgeGroupsWithTolBS() builds an initial age group (0,0); when the input
# ages start at 1 this group is empty, its local `ids` is never assigned, and
# the expansion step then reads `ids` before assignment (used in the
# removeDuplIDs() calls). Under future multisession the variable cannot be
# resolved in the workers ("object 'ids' not found"). A global `ids` is NOT
# exported by future's globals scanner (the scanner classifies `ids` as local
# because the function also assigns it), so the fix must live inside the
# function body: inject a local `ids <- character(0)` at the top of the
# function IN MEMORY. The protected pipeline source file is NOT modified.
# With an empty `ids`, the duplicate-ID screen is a no-op for the empty
# group, which is the behaviorally neutral choice (PIDs are character).
patch_ids_init <- function(fn) {
  b <- as.list(body(fn))
  body(fn) <- as.call(append(b, list(quote(ids <- character(0))), after = 1))
  fn
}
defineAgeGroupsWithTolBS <- patch_ids_init(defineAgeGroupsWithTolBS)
message("[PATCH] defineAgeGroupsWithTolBS: local ids initialization injected.")

outputDir  <- "cont_out/"
NBootstrap <- 200
log_file   <- "data/processed/REVIEW_gamlss_rerun_log_REV.txt"

log_msg <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
  message(txt)
  cat(txt, "\n", file = log_file, append = TRUE)
}

log_msg("R version: ", R.version.string,
        " | cores: ", detectCores(logical = TRUE),
        " | NBootstrap target: ", NBootstrap)

analytes <- list(
  list(csv = "data/coa_results/CSV cont PT 2025.10.17.csv",   filename = "PT"),
  list(csv = "data/coa_results/CSV cont aPTT 2025.10.17.csv", filename = "aPTT"),
  list(csv = "data/coa_results/CSV cont fib 2025.10.17.csv",  filename = "Fib")
)

for (a in analytes) tryCatch({

  log_msg("[", a$filename, "] reading input: ", a$csv)
  inputData <- read.table(file = a$csv, header = TRUE, sep = ";",
                          stringsAsFactors = FALSE)
  log_msg("[", a$filename, "] input rows: ", nrow(inputData))

  # Load cached point estimate (runPipelinePointEst returns it without
  # refitting because the .RData files already exist in cont_out/).
  pointEstGamlss <- runPipelinePointEst(inputData = inputData,
                                        covarName = "Age",
                                        colnameID = "PID",
                                        colnameValue = "Value",
                                        outputDir = outputDir,
                                        filename = a$filename,
                                        NCores = NULL)

  # Input-consistency check: the point-estimate object stores the data it was
  # fitted on; the bootstrap must resample the same data set.
  stored_data <- tryCatch(pointEstGamlss$.user$data, error = function(e) NULL)
  if (!is.null(stored_data) && !is.null(nrow(stored_data))) {
    if (nrow(stored_data) != nrow(inputData)) {
      log_msg("[", a$filename, "] FATAL: input rows (", nrow(inputData),
              ") != rows stored in point estimate (", nrow(stored_data),
              "). Wrong input file? ABORTING this analyte.")
      stop("Input-consistency check failed for ", a$filename)
    }
    log_msg("[", a$filename, "] input-consistency check PASSED (",
            nrow(stored_data), " rows).")
  } else {
    log_msg("[", a$filename, "] input-consistency check SKIPPED: point ",
            "estimate object does not expose fitted data; proceeding with ",
            "the runner-referenced input file.")
  }

  seed <- 123
  set.seed(seed)

  t0 <- Sys.time()
  log_msg("[", a$filename, "] bootstrap start (existing iterations reused, ",
          "missing ones computed).")

  bootstrapRes <- runPipelineCIs(inputData = inputData,
                                 pointEstGamlss = pointEstGamlss,
                                 NBootstrap = NBootstrap,
                                 covarName = "Age",
                                 colnameID = "PID",
                                 colnameValue = "Value",
                                 outputDir = outputDir,
                                 filename = a$filename,
                                 NCores = NULL)

  elapsed <- difftime(Sys.time(), t0, units = "mins")
  n_ok <- sum(vapply(bootstrapRes, function(x) !is.null(x$Models), logical(1)))
  log_msg("[", a$filename, "] bootstrap done: ", length(bootstrapRes),
          " iterations (", n_ok, " with successful GAMLSS fit) in ",
          round(as.numeric(elapsed), 1), " min.")
}, error = function(e) {
  log_msg("[", a$filename, "] ERROR: ", conditionMessage(e),
          " -- continuing with next analyte.")
})

log_msg("ALL DONE.")
