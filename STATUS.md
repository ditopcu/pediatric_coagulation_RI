# STATUS.md -- Pediatric Coagulation RI Project

> **Last updated:** 2026-06-26
> **Updated by:** Claude Code (VS Code session)

## Session 2026-06-26

- **GitHub repository preparation** (target: <https://github.com/ditopcu/pediatric_coagulation_RI>).
  - Renamed raw-data folder `data/coa_tuba/` -> `data/coa_results/`. Replaced
    all 27 references across `src/`, `src/_archive/`, root config, docs, and
    audit logs in a single in-place pass.
  - Archived 11 dev/debug scripts to `src/_archive/` with a per-file purpose
    table (`src/_archive/README.md`). Active `src/` now contains only the
    analysis pipeline (21 scripts).
  - Added `.gitignore` (raw data, working .txt extractions, IDE/OS scratch,
    memory/.claude/), `LICENSE` (MIT), and a top-level `README.md` describing
    reproducibility model, environment, layout, data schema, and run order.
  - Added a small synthetic dataset under `data/synthetic/` plus generator
    `src/_make_synthetic_data.R` for end-to-end smoke testing (numerical results
    will NOT match the manuscript -- pipeline integrity check only).
  - Captured pinned R session via `docs/sessionInfo.txt`
    (`src/_capture_sessioninfo.R` to refresh).
  - Marked `src/tables.R` Supplementary Table S2 (analytical performance) values
    as illustrative placeholders -- authoritative QC is reported directly in the
    manuscript supplement.

## Session 2026-06-23

- **Manuscript switched from 3SD-filtered to RAW data** for FP rate and Guven
  transferability analyses (Analyses 1, 2, 3, 6). 3SD filter only kept for
  Figure 1 visualisation. New script: `src/analysis_raw_comparison.R`. Side-by-
  side comparison exported in `tables/COMPARISON_3sd_vs_raw.xlsx`. Numerical
  differences small (~0.4-1.5 pp higher FP rates with RAW); main conclusions
  unchanged.
- Manuscript + Tables.docx fully synchronised with RAW results (Abstract,
  Results, Discussion, Tables 1/3/4). Verification report:
  `article_update/FINAL_CHANGES.txt`. Remaining manual edits noted there
  (age mean/SD, Harris-Boyd z).

### Infant (30-365 day) supplementary analysis -- COMPLETE
- New scripts:
  - `src/infant_A_descriptive.R` -- pooled descriptive statistics
  - `src/infant_B_refineR.R` + `src/infant_B_finalize.R` -- refineR RI
    (NBootstrap = 200, 90% CI; tolerance interval)
  - `src/infant_C_classification.R` -- Below/Within/Above against
    pediatric refineR RI
- Source data: `data/coa_results/{PT-2,aPTT,Fib-2} 2025.10.17.xlsx`
  (filter: `age_year >= 30/365 & age_year < 1.0`).
- Sample sizes: PT 1,367 | aPTT 1,365 | Fibrinogen 217.
- Infant refineR RI (vs pediatric 1-18 y / 1-12 y for aPTT):
  - PT: 8.36 (8.19-8.43) to 11.5 (11.3-11.7)  vs  8.56-10.9
  - aPTT: 23.8 (23.2-24.3) to 40.8 (39.3-41.6)  vs  24.3-35.3
  - Fib: 1.61 (1.52-1.78) to 3.34 (2.85-3.45)  vs  1.89-3.79
- Classification (infant results vs pediatric RI):
  - PT 21.7% outside | aPTT 26.1% outside | Fib 19.8% outside
  - Pattern matches developmental hemostasis (longer PT/aPTT, lower Fib in infants).
- Deliverables:
  - Excel: `data/processed/tables/SUPPLEMENT_infant_30_365d.xlsx`
    (Descriptive, Classification_vs_RI, Infant_refineR_RI, Pediatric_vs_Infant)
  - Figure: `figures/{TIFF_600DPI,PNG_300DPI}/figure_supplement_infant_RI.{tiff,png}`
  - Captions: `docs/figure_captions.md` (Supplementary Figure S2 and tables).
- RDS fits saved at `data/processed/infant_refineR_{PT,aPTT,Fibrinogen}.RDS`.

### Data provenance note
- `data_test/` files contain a `department` column with full granularity
  (~25 outpatient clinics + ER + inpatient services). Used data_test/ Sheet 1
  -- NOT data_test/ Sayfa6 (which is a pivot summary).
- ID overlap between `data_test/` and `data/coa_results/tce 2026 1-18 *.xlsx`
  is PT 94.4% / aPTT 88.3% / Fib 96.8%. Datasets are not identical snapshots;
  decision: keep current manuscript dataset for main analyses; use data_test/
  only for joining the `department` column (Analysis 4).

## Session 2026-05-27

- Built Excel deliverables for tables:
  - `data/processed/tables/PUB_main_tables.xlsx` -- 6 sheets (T1-T6)
  - `data/processed/tables/PUB_supplementary_tables.xlsx` -- 3 sheets (S1-S3)
  - Generators: `src/export_main_xlsx.R`, `src/export_supplementary_xlsx.R`
- **Analysis 2 enhancement**: added exact binomial proportion test (one-sided,
  H0: p <= 0.05) per age group, applied to the full age-group sample (not the
  20-subject subset). Complements the EP28-A3c bootstrap pass rate. New columns
  in `analysis2_guven_transferability.csv` / `PUB_table_4`:
  `n_outside`, `prop_outside`, `binom_p`, `binom_reject`.
  - Result: 17/17 ages reject H0 for both PT and aPTT -- consistent with
    bootstrap finding that Guven RIs are not transferable.

---

## What was done this session

### Bug fixes
1. **BUG-001 (CRITICAL)**: GAMLSS curves showed wrong values. Root cause: (a) `ptrans`/`families`/`fam` not in .GlobalEnv; (b) RIMat 7-row indexing bug (read row 1/2/3 as 2.5th/10th/25th instead of 2.5th/50th/97.5th). Fixed with rowname-based access. Verified against thesis values -- PASS.
2. **BUG-002**: Emoji in R messages replaced with ASCII.
3. **BUG-003**: geom_text clipping fixed (x=17.2, coord_cartesian clip=off).

### Code organization
4. **Refactored into 5 scripts:**
   - `src/common.R` -- shared libraries, design system (PAL, theme), data loading, all RI definitions, GAMLSS continuous RI extraction
   - `src/pub.R` -- 6 publication figures (no embedded titles, semantic colors)
   - `src/analysis.R` -- 6 analyses + 1 figure (FP curves)
   - `src/tables.R` -- 18 CSV tables (thesis + publication + supplementary)
   - `src/run_all.R` -- master runner (common -> pub -> analysis -> tables, ~1 min)

### Analyses completed
5. **Analysis 1**: False-positive rates by RI source (4 sources x 3 tests)
6. **Analysis 2**: Guven et al. RI transferability (CLSI EP28-A3c, 20-sample bootstrap)
7. **Analysis 3**: Age-specific FP rate curves (figure)
8. **Analysis 5**: International comparison table (6 studies)
9. **Analysis 5b**: Harris-Boyd age partitioning test (aPTT z=18.16, justified)
10. **Analysis 6**: Bootstrap 95% CI for FP rates (B=1000)

### Figure design
11. Continuous RI colors: red -> green (PAL$contin) for semantic consistency
12. All embedded titles/subtitles removed -> `docs/figure_captions.md`
13. Redundant figures removed (sekil_4_4, figure_RI_all)
14. Figure 5: formal legend box added (4 methods)
15. Fibrinogen continuous RI added to all relevant figures

### Tables
16. 18 CSV tables generated: 9 thesis (TEZ_tablo_4_1-4_9), 6 publication (PUB_table_1-6), 3 supplementary (PUB_table_S1-S3)

## Key results

### False-Positive Rates (95% Bootstrap CI)

| Test | Manufacturer | refineR | Continuous (GAMLSS) | Guven et al. |
|------|-------------|---------|-------------------|-------------|
| PT | 17.0% (16.4-17.5) | 13.1% (12.7-13.6) | **11.4% (10.9-11.9)** | 29.3% (28.7-30.0) |
| aPTT | **31.9% (31.2-32.6)** | 10.9% (10.5-11.4) | **10.0% (9.5-10.5)** | 14.9% (14.4-15.4) |
| Fibrinogen | 14.3% (13.1-15.6) | 19.3% (18.1-20.7) | 17.9% (16.6-19.3) | -- |

### Guven Transferability
- PT: 0/17 ages pass (mean pass rate 5.6%) -- **not transferable**
- aPTT: 0/17 ages pass (mean pass rate 43.1%) -- **not transferable**

### Harris-Boyd Partitioning
- aPTT at age 12: z = 18.16 (threshold = 3) -- **partition justified**

## Final output inventory

**Figures (6 TIFF + 6 PNG):**
figure_1_distribution, figure_2_partitioning, figure_3_indirect_vs_direct,
figure_4_continuous_RI, figure_5_combined_three_methods, figure_fp_curves_by_age

**Tables (18 CSV in data/processed/tables/):**
TEZ_tablo_4_1_sample_sizes, TEZ_tablo_4_2_descriptive_PT, TEZ_tablo_4_3_descriptive_aPTT,
TEZ_tablo_4_4_descriptive_Fib, TEZ_tablo_4_5_analytical_performance,
TEZ_tablo_4_6_refineR_PT_by_age_sex, TEZ_tablo_4_7_refineR_aPTT_by_age_sex,
TEZ_tablo_4_8_refineR_Fib_by_age, TEZ_tablo_4_9_indirect_vs_direct,
PUB_table_1_study_population, PUB_table_2_reference_intervals,
PUB_table_3_false_positive_rates, PUB_table_4_guven_transferability,
PUB_table_5_international_comparison, PUB_table_6_harris_boyd,
PUB_table_S1_descriptive_all, PUB_table_S2_refineR_all_by_age_sex,
PUB_table_S3_continuous_RI_values

**Analysis CSV (7 in data/processed/):**
analysis1_fp_by_age, analysis1_fp_summary, analysis2_guven_transferability,
analysis3_fp_curves_data, analysis5_international_comparison,
analysis5b_partition_harris_boyd, analysis6_fp_bootstrap_ci

## CI index fix (this session, continued)

19. **CI index bug fixed**: `estimateCIs()` with `RIperc=c(0.025,0.5,0.975)` produced 3-element CILow/CIHigh but 7-row RIMat, causing CI values to map to wrong percentiles (p50 CI showed p10 CI, p97.5 CI showed p25 CI). Fix: pass all 7 default percentiles, use indices 1/4/7 for 2.5th/50th/97.5th. Applied in `common.R`, `pub.R`, `tables.R`. Verified: PT age 1 p97.5 CI = [10.82, 10.99] (was [9.06, 9.13]).

20. **Figure 2 partitioning redesign:**
    - Facet order reversed: Upper RL (97.5th) on top, Lower RL (2.5th) on bottom
    - Overall RI bands added: orange horizontal line + shaded 90% CI band (from ri_comparison)
    - aPTT: two separate bands for 1-12 y and 12-18 y partitions
    - "All" point added at right edge (combined estimate, orange diamond)
    - Figure caption updated in docs/figure_captions.md

## Known bugs / blockers

None. All previously known bugs resolved.

## Next steps

1. **Exclusion sensitivity analysis** (Analysis 4): Department column available -- pending decision
2. **Target journal selection**
3. **Manuscript drafting**

## Reference documents

- `docs/RESULTS.md` -- comprehensive results summary with all tables and figures
- `docs/figure_captions.md` -- figure titles and captions
- Thesis: `Tuba_Cakmak_Ercan_TEZ.pdf`
- Guven et al. 2026, Scand J Clin Lab Invest
- Weidhofer et al. 2018, Clin Chim Acta
- Luo et al. 2026, Health Sci Rep

---

> **Convention**: Claude Code updates this file at the end of every session.
