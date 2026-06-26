# Results Summary

> Pediatric Coagulation Reference Intervals: PT, aPTT, Fibrinogen
> Turkish tertiary center, ages 1-18 years, Roche Cobas t511/t711
> Methods: refineR (indirect), GAMLSS (continuous), Direct (a posteriori)

---

## Data

| Parameter | N (raw) | N (3SD filtered) | Age range | Analyzer |
|-----------|---------|-------------------|-----------|----------|
| PT | 17,000 | 16,907 | 1-18 years | Cobas t511/t711 |
| aPTT | 17,001 | 16,900 | 1-18 years | Cobas t511/t711 |
| Fibrinogen | 3,166 | 3,111 | 1-18 years | Cobas t511/t711 |
| **Total** | **37,167** | **36,918** | | |

Exclusion: ICU patients removed at data extraction. 3SD global filter applied per test.

---

## Figures

### Figure 1 — Distribution of coagulation test results by age and sex
**File:** `sekil_4_1_distribution`
3-panel boxplot + jitter (PT, aPTT, Fibrinogen) after 3SD exclusion. Shows raw data spread by integer age, colored by sex (Female=red, Male=blue). Demonstrates minimal sex-specific differences across all three tests.

### Figure 2 — Sex-specific partitioning CI (PT and aPTT)
**File:** `sekil_4_2_4_3_partitioning`
refineR 90% CI errorbar plots for 2.5th and 97.5th percentiles by age and sex. CI overlap between sexes at most ages supports combined (non-sex-specific) RI for both PT and aPTT.

### Figure 3 — Indirect vs Direct RI comparison
**File:** `sekil_4_4_indirect_vs_direct`
Crossbar plot comparing refineR (indirect) and Direct (a posteriori) RIs with 90% CI for all test/age-group combinations. Dashed lines show manufacturer limits. Demonstrates good agreement between methods.

### Figure 4 — Continuous reference intervals (GAMLSS)
**File:** `sekil_4_5_4_6_4_7_continuous_RI` and `figure_1_continuous_RI_publication`
3-panel continuous RI curves (2.5th, 50th, 97.5th percentiles) with 90% CI bands for PT, aPTT, and Fibrinogen. Generated from GAMLSS pipeline (BCCG family, power transformation). Age-dependent dynamics clearly visible for aPTT (decreasing with age) and Fibrinogen (U-shape).

**Key values (age 1 / age 17):**

| Test | 2.5th | 50th | 97.5th |
|------|-------|------|--------|
| PT | 8.4-8.5 | 9.5-9.7 | 10.9-11.2 |
| aPTT | 24.2→22.8 | 26.1→24.6 | 35.5→33.5 |
| Fibrinogen | 1.6→1.9 | 2.5→2.7 | 3.5→3.9 |

### Figure 5 — Indirect vs Direct with manufacturer overlay
**File:** `figure_2_indirect_vs_direct_publication`
Extended version of Figure 3 with manufacturer reference limits as dashed vertical lines. Highlights how manufacturer RI for aPTT (23.6-30.6s) is narrower than both indirect and direct estimates.

### Figure 6 — Three methods combined
**File:** `figure_3_combined_three_methods`
3-panel overlay: continuous GAMLSS curves as background, refineR segments (orange solid), Direct segments (blue dashed), manufacturer limits (gray dotted) for all three tests. Single figure showing all RI sources simultaneously.

### Figure 7 — Age-specific false-positive rates
**File:** `figure_fp_curves_by_age`
3-panel line plot showing false-positive rate (%) by age for each RI source (manufacturer, refineR, continuous GAMLSS, Guven et al.). Dotted line at expected 5%. Key finding: manufacturer aPTT RI causes up to 40% FP at some ages.

### Figure 8 — RI comparison (from thesis data)
**File:** `figure_RI_all_comparison_updated`
Crossbar comparison from `RI Tuba Tez Final All.xlsx`, showing refineR vs Direct across all age subgroups.

---

## Tables

### Table 1 — Reference intervals by method

| Test | Age group | Method | N | Lower (90% CI) | Upper (90% CI) |
|------|-----------|--------|---|-----------------|-----------------|
| PT | 1-18 y | refineR | 17,000 | 8.56 (8.54-8.58) | 10.90 (10.70-10.90) |
| PT | 1-18 y | Direct | 305 | 8.69 (8.55-8.84) | 10.80 (10.60-11.20) |
| aPTT | 1-12 y | refineR | 10,091 | 24.30 (23.30-24.50) | 35.30 (34.40-35.40) |
| aPTT | 1-12 y | Direct | 170 | 22.96 (22.40-23.40) | 34.74 (34.20-35.60) |
| aPTT | 12-18 y | refineR | 6,910 | 23.30 (22.70-23.70) | 33.90 (33.10-34.10) |
| aPTT | 12-18 y | Direct | 150 | 24.10 (23.50-24.30) | 32.69 (31.60-33.40) |
| Fibrinogen | 1-18 y | refineR | 3,166 | 1.89 (1.84-1.92) | 3.79 (3.55-4.06) |
| Fibrinogen | 1-18 y | Direct | 198 | 2.10 (2.07-2.16) | 3.56 (3.37-3.89) |

### Table 2 — Manufacturer reference intervals

| Test | Lower | Upper | Unit |
|------|-------|-------|------|
| PT | 8.4 | 10.6 | seconds |
| aPTT | 23.6 | 30.6 | seconds |
| Fibrinogen | 1.93 | 4.12 | g/L |

### Table 3 — False-positive rates by RI source (with 95% bootstrap CI)

| Test | Manufacturer | refineR | Continuous (GAMLSS) | Guven et al. |
|------|-------------|---------|-------------------|-------------|
| PT | 17.0% (16.4-17.5) | 13.1% (12.7-13.6) | **11.4% (10.9-11.9)** | 29.3% (28.7-30.0) |
| aPTT | **31.9% (31.2-32.6)** | 10.9% (10.5-11.4) | **10.0% (9.5-10.5)** | 14.9% (14.4-15.4) |
| Fibrinogen | 14.3% (13.1-15.6) | 19.3% (18.1-20.7) | 17.9% (16.6-19.3) | -- |

**Key finding:** Manufacturer aPTT RI produces 31.9% false-positives vs 10.0% with continuous GAMLSS (3.2x higher). All CIs non-overlapping — difference is statistically significant.

**File:** `data/processed/analysis1_fp_summary.csv`, `analysis1_fp_by_age.csv`, `analysis6_fp_bootstrap_ci.csv`

### Table 4 — Guven et al. RI transferability (CLSI EP28-A3c)

| Test | Ages tested | Ages passed | Ages failed | Mean pass rate |
|------|-------------|-------------|-------------|----------------|
| PT | 17 | **0** | 17 | 5.6% |
| aPTT | 17 | **0** | 17 | 43.1% |

Method: 20-sample bootstrap validation (B=200). Pass criterion: <=2/20 (10%) outside RI at >=90% of iterations. **Result: Guven et al. RIs are not transferable to our population for either test.**

**File:** `data/processed/analysis2_guven_transferability.csv`

### Table 5 — Age partitioning justification (Harris-Boyd test)

| Test | Cutpoint | n (<12y) | n (>=12y) | Mean diff | z-score | Justified? |
|------|----------|----------|-----------|-----------|---------|-----------|
| aPTT | 12 years | 10,013 | 6,887 | 0.95 s | **18.16** | **YES** |
| PT | 12 years | 8,696 | 8,211 | 0.06 s | 5.02 | YES* |

*PT z-score exceeds threshold but clinical difference is minimal (0.06s). Partition used only for aPTT.

**File:** `data/processed/analysis5b_partition_harris_boyd.csv`

### Table 6 — International comparison

**File:** `data/processed/analysis5_international_comparison.csv`

Selected comparisons (1-18 year range, PT and aPTT in seconds):

| Study | Analyzer | Method | PT lower | PT upper | aPTT lower | aPTT upper |
|-------|----------|--------|----------|----------|------------|------------|
| **This study (refineR)** | Cobas t511/t711 | Indirect | 8.56 | 10.90 | 23.3-24.3 | 33.9-35.3 |
| **This study (Direct)** | Cobas t511/t711 | Direct | 8.69 | 10.80 | 23.0-24.1 | 32.7-34.7 |
| Guven et al. 2026 | Cobas t511 | Indirect | 7.8-8.2 | 10.0-10.7 | 24.3-25.4 | 32.4-36.2 |
| Luo et al. 2026 | Sysmex CN-6000 | Direct | 10.1-10.4 | 12.4-13.0 | 23.4-24.7 | 32.6-32.9 |
| Weidhofer et al. 2018 | STA-Compact/BCS-XP | Indirect | (% units) | (% units) | ~26 | ~42 |
| Manufacturer (Roche) | Cobas t511/t711 | Insert | 8.4 | 10.6 | 23.6 | 30.6 |

Note: PT values are analyzer-dependent (Cobas vs Sysmex differ by ~2s). aPTT shows more inter-study variation, especially in upper limits.

---

## Output file inventory

### Figures (TIFF 600 DPI + PNG 300 DPI)

| Filename | Description |
|----------|-------------|
| `sekil_4_1_distribution` | Boxplot + jitter, 3 tests by age/sex |
| `sekil_4_2_4_3_partitioning` | Sex-specific partition CI (PT + aPTT) |
| `sekil_4_4_indirect_vs_direct` | Indirect vs Direct crossbar |
| `sekil_4_5_4_6_4_7_continuous_RI` | GAMLSS continuous RI (3 panels) |
| `figure_1_continuous_RI_publication` | Publication Figure 1 |
| `figure_2_indirect_vs_direct_publication` | Publication Figure 2 |
| `figure_3_combined_three_methods` | Three methods overlay |
| `figure_fp_curves_by_age` | Age-specific FP rate curves |
| `figure_RI_all_comparison_updated` | RI comparison from thesis data |

### CSV tables

| Filename | Description |
|----------|-------------|
| `analysis1_fp_summary.csv` | Overall FP rates per test/source |
| `analysis1_fp_by_age.csv` | FP rates by age/test/source |
| `analysis2_guven_transferability.csv` | CLSI EP28-A3c validation results |
| `analysis3_fp_curves_data.csv` | FP curve data (long format) |
| `analysis5_international_comparison.csv` | International RI comparison table |
| `analysis5b_partition_harris_boyd.csv` | Harris-Boyd test results |
| `analysis6_fp_bootstrap_ci.csv` | FP rates with 95% bootstrap CI |

---

## How to reproduce

```r
setwd("path/to/ISH_Tuba_COA")
source("src/run_all.R")  # runs common.R -> pub.R -> analysis.R (~1.5 min)
```

Individual scripts can also be run separately:
```r
source("src/common.R")    # load data + RI definitions
source("src/pub.R")       # generate figures only
source("src/analysis.R")  # run analyses only
```
