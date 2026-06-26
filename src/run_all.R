# ==============================================================================
# RUN ALL: Execute the full pipeline in order
# Pediatric Coagulation Reference Intervals
# ==============================================================================
# Usage: setwd("path/to/ISH_Tuba_COA") then source("src/run_all.R")
#
# Pipeline:
#   1. common.R   — Load libraries, design system, data, RI definitions
#   2. pub.R      — Generate all publication figures
#   3. analysis.R — Run manuscript analyses (FP rates, transferability, etc.)
#   4. tables.R   — Generate all thesis + publication tables as CSV
#
# Prerequisites:
#   - refineR RDS files in data/processed/
#   - GAMLSS pipeline outputs in cont_out/
#   - Raw data in data/coa_results/tce 2026 1-18 *.xlsx
# ==============================================================================

message("==============================================================")
message("[START] Pediatric Coagulation RI Pipeline")
message("==============================================================\n")

t0 <- Sys.time()

# Step 1: Common setup
message("--- Step 1/4: Common setup ---")
source("src/common.R")

# Step 2: Publication figures
message("\n--- Step 2/4: Publication figures ---")
source("src/pub.R")

# Step 3: Manuscript analyses
message("\n--- Step 3/4: Manuscript analyses ---")
source("src/analysis.R")

# Step 4: Tables
message("\n--- Step 4/4: Tables ---")
source("src/tables.R")

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)

message("\n==============================================================")
message(sprintf("[DONE] Full pipeline completed in %s minutes.", elapsed))
message("  Figures: figures/TIFF_600DPI/ and figures/PNG_300DPI/")
message("  Tables:  data/processed/tables/")
message("  Analyses: data/processed/analysis*.csv")
message("==============================================================")
