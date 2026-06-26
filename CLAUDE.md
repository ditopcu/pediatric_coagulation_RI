# Pediatric Coagulation Reference Intervals

## Context [REQUIRED]

Establishing age- and sex-specific reference intervals (RIs) for PT, aPTT, and fibrinogen in children aged 1-18 years, using retrospective hospital data from a Turkish tertiary center (Roche Cobas t511/t711 analyzers). The work combines indirect RI estimation (refineR), continuous RI modeling (GAMLSS via the Ammer pipeline), and direct method validation per CLSI EP28-A3c. The goal is a peer-reviewed publication that prioritizes clinical impact over methodological complexity.

## Stack [REQUIRED]

- Language: R 4.4.x / 4.5.x on Windows (RStudio / Positron / VS Code)
- Key packages: tidyverse, readxl, refineR, gamlss, patchwork, ggrepel, cowplot, janitor, moments
- Pipeline: Ammer et al. GAMLSS pipeline (RUNNER_Pipeline.R, algoRICurves.R, plotRICurve.R, utils.R)
- Figures: Tufte/Annesley academic design system (theme_tufte_academic in common.R)
- Output: TIFF 600 DPI + PNG 300 DPI dual export

## Folder structure

```
ISH_Tuba_COA/
├── data/
│   ├── coa_results/                # Raw Excel data (DO NOT MODIFY)
│   │   ├── tce 2026 1-18 pt.xlsx      (17000 rows, 12 cols incl. department)
│   │   ├── tce 2026 1-18 aptt.xlsx    (17001 rows)
│   │   └── tce 2026 1-18 fib.xlsx     (3166 rows)
│   └── processed/
│       ├── tables/              # 18 CSV tables (TEZ_ + PUB_)
│       ├── analysis*.csv        # 7 analysis output CSVs
│       ├── ptz_ref_*.RDS        # refineR intermediate results
│       ├── aptt_ref_*.RDS
│       └── fib_ref_*.RDS
├── src/
│   ├── common.R                 # Shared: libraries, design system, data, RI definitions, GAMLSS extraction
│   ├── pub.R                    # Publication figures (6 figures, sources common.R)
│   ├── analysis.R               # Manuscript analyses 1-6 (sources common.R)
│   ├── tables.R                 # Table generation — 18 CSVs (sources common.R)
│   ├── run_all.R                # Master runner: common -> pub -> analysis -> tables (~1 min)
│   ├── coa hesap 3.R            # Main refineR analysis (age/sex partitioning)
│   ├── plots final.R            # Original RI comparison plots (legacy)
│   ├── RUNNER_Pipeline aptt.R   # Pipeline runner for aPTT
│   ├── RUNNER_Pipeline pt.R     # Pipeline runner for PT
│   ├── RUNNER_Pipeline fib.R    # Pipeline runner for Fibrinogen
│   └── pipeline/                # Ammer pipeline source (DO NOT MODIFY)
│       ├── RUNNER_Pipeline.R
│       ├── algoRICurves.R
│       ├── plotRICurve.R
│       ├── utils.R
│       └── mode_deneme.R
├── cont_out/                    # GAMLSS pipeline outputs (DO NOT MODIFY)
│   ├── PT_gamlss_PointEst.RData / PT_gamlssModel_Est_1..5.RData
│   ├── aPTT_gamlss_PointEst.RData / aPTT_gamlssModel_Est_1..5.RData
│   └── Fib_gamlss_PointEst.RData / Fib_gamlssModel_Est_1..5.RData
├── figures/
│   ├── TIFF_600DPI/             # 6 TIFF figures
│   └── PNG_300DPI/              # 6 PNG figures
├── docs/
│   ├── figure_captions.md       # Figure titles and captions (not embedded in figures)
│   ├── RESULTS.md               # Comprehensive results summary
│   └── manuscript/
├── _archive/                    # Old versions (never deleted, just moved)
├── CLAUDE.md                    # This file
└── STATUS.md                    # Current progress
```

## Working rules [REQUIRED]

- Code and code comments in English only, no emojis in R code or messages
- Conversation in Turkish unless specified otherwise
- Never modify files in `data/coa_results/`, `cont_out/`, or `src/pipeline/`
- Before writing or overwriting any file, confirm path and intent
- Do not refactor or rename code that was not part of the request
- Preserve existing variable and function names unless explicitly asked
- Never make assumptions about data structure -- inspect first, then act
- No unsolicited commentary on style, performance, or architecture
- All figure outputs: dual save as TIFF 600 DPI + PNG 300 DPI
- Figure titles/subtitles NOT embedded -- use docs/figure_captions.md
- Working directory is always the project root
- At the end of every session, update `STATUS.md` with what was done, what changed, and what is next

## Domain rules

- Follow CLSI EP28-A3c for direct RI validation and IFCC/EFLM for indirect methods
- Never include real patient identifiers in outputs, examples, or logs
- refineR RIs use 2.5th and 97.5th percentiles with 90% CI (bootstrap N=200)
- GAMLSS continuous RIs use BCCG family with power transformation (ptrans)
- **Critical**: GAMLSS `centiles.pred()` requires `pp`, `ptrans`, `families`, `fam`, `covarName` in .GlobalEnv
- **Critical**: `getRICurve()` returns 7-row RIMat (all percentiles) -- always index by rowname (`"Perc_0.025"`)
- **Critical**: `estimateCIs()` CILow/CIHigh are indexed 1..length(RIperc). Always pass all 7 default percentiles and use indices 1/4/7 for 2.5th/50th/97.5th CIs
- Data column convention: `sample_id`, `sex` (K=Female, E=Male), `age_year`, `result_num`
- `age_int = floor(age_year)` for integer age grouping

## Color system (semantic, consistent across ALL figures)

```
refineR (indirect)  = PAL$refiner  = "#E67E22" (orange)
Direct              = PAL$direct   = "#2980B9" (blue)
Continuous (GAMLSS) = PAL$contin   = "#27AE60" (green)
Manufacturer        = PAL$base2    = "#BDC3C7" (gray)
Guven et al.        = PAL$accent3  = "#8E44AD" (purple)
Female              = PAL$female   = "#C0392B" (red)
Male                = PAL$male     = "#2980B9" (blue)
```

## Data notes

- Excel files use Turkish decimal separator (comma) -- readxl handles this automatically
- aPTT dataset contains a likely data entry error: Value = 0.000477 (ProbNP = 0)
- Age 18 has very few observations (1 in aPTT weighted GAMLSS data) -- edge estimates unreliable
- 3SD filter is applied per test globally (not per age group) in current scripts
- Pipeline input requires tab-delimited file with columns: Age (integer years), Value, PID
- ICU patients already excluded at Excel stage; COCUK ACIL (39%) still in data
- Guven et al. (2026, Scand J Clin Lab Invest) is the key competing study -- same Cobas t511, Turkish pediatric, PT/INR/aPTT (no fibrinogen)

## Open questions / decisions

- [x] pp/ptrans fix verified against thesis values? -- YES, PASS
- [x] Include fibrinogen continuous RI? -- YES
- [x] Analysis 1-3, 5, 5b, 6 completed
- [x] Analysis 5: International comparison table completed (Weidhofer, Luo, Guven)
- [ ] Analysis 4: Exclusion sensitivity -- department column available, pending decision
- [ ] Target journal finalized?

## Current state

See `STATUS.md` for the latest progress and next steps.

## How to run

```r
# Set working directory to project root
setwd("path/to/ISH_Tuba_COA")

# Run full pipeline (~1 min): figures + analyses + tables
source("src/run_all.R")

# Or run individually:
source("src/common.R")     # data + RI definitions
source("src/pub.R")        # 6 figures
source("src/analysis.R")   # 6 analyses + 1 figure
source("src/tables.R")     # 18 tables

# Main refineR analysis (generates RDS files in data/processed/)
source("src/coa hesap 3.R")
```

## Out of scope

- No production deployment -- research code only
- Do not optimize for performance unless asked -- clarity first
- Do not add multi-algorithm comparison (explicitly rejected in favor of clinical framing)
- Do not create new statistical methods -- use established refineR + GAMLSS framework
- No web apps, Shiny dashboards, or interactive outputs
