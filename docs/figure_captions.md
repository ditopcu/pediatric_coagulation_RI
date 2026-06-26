# Figure Captions

## Figure 1 — figure_1_distribution
**Title:** Distribution of coagulation test results by age and sex after outlier exclusion.
**Caption:** Boxplot with jitter overlay for PT (A), aPTT (B), and fibrinogen (C) in children aged 1-18 years after exclusion of values beyond +/-3 SD. Female (red) and male (blue) results shown side by side at each integer age. Horizontal lines within boxes indicate medians; box edges indicate 25th and 75th percentiles. Individual data points shown with transparency. Roche Cobas t511/t711; N = 36,918.

## Figure 2 — figure_2_partitioning
**Title:** Sex-specific reference limit estimates with 90% confidence intervals by age.
**Caption:** refineR-estimated upper (97.5th, top row) and lower (2.5th, bottom row) reference limits for PT (A) and aPTT (B) by integer age and sex. Error bars represent 90% bootstrap confidence intervals (N = 200). Orange horizontal lines and shaded bands indicate the overall (all ages combined) reference limit point estimate and 90% CI, respectively. For aPTT, separate bands are shown for the 1-12 y and 12-18 y age partitions. The rightmost point ("All") represents the combined estimate across all ages. Overlap of sex-specific confidence intervals at most ages supports combined (non-sex-specific) reference intervals. Female (red), Male (blue).

## Figure 3 — figure_3_indirect_vs_direct
**Title:** Comparison of indirect and direct reference intervals with manufacturer limits.
**Caption:** Crossbar plot comparing indirect (refineR, orange) and direct (a posteriori, blue) reference intervals for PT, aPTT, and fibrinogen by age group. Error bars represent 90% confidence intervals for lower and upper reference limits. Dashed vertical lines indicate manufacturer-supplied reference limits (Roche Cobas t511/t711). For aPTT, separate age-partitioned intervals (1-12 y and 12-18 y) are shown based on Harris-Boyd partitioning analysis.

## Figure 4 — figure_4_continuous_RI
**Title:** Age-dependent continuous reference intervals estimated by GAMLSS.
**Caption:** Continuous reference interval curves for PT (A), aPTT (B), and fibrinogen (C) in children aged 1-18 years. Solid lines represent the 2.5th, 50th (median, thick), and 97.5th percentiles estimated using GAMLSS (BCCG family with power transformation). Shaded bands represent 90% bootstrap confidence intervals (5 iterations). Direct labels identify each percentile curve. Roche Cobas t511/t711.

## Figure 5 — figure_5_combined_three_methods
**Title:** Overlay comparison of three reference interval estimation methods.
**Caption:** Combined display of continuous GAMLSS curves (green), indirect refineR reference limits (orange, solid segments), direct a posteriori limits (blue, dashed segments), and manufacturer limits (gray, dotted lines) for PT (A), aPTT (B), and fibrinogen (C). Each panel shows all four RI sources on the same scale to facilitate direct comparison. For aPTT and fibrinogen, age-dependent variation of the continuous RI contrasts with the single-interval approaches.

## Figure 6 — figure_fp_curves_by_age
**Title:** Age-specific false-positive rates by reference interval source.
**Caption:** Proportion of 3SD-filtered results classified as outside reference limits at each integer age for PT (A), aPTT (B), and fibrinogen (C). Lines represent different RI sources: manufacturer (gray, dashed), Guven et al. (purple, dot-dash), indirect refineR (orange, solid), and continuous GAMLSS (green, solid). Dotted horizontal line indicates the expected 5% false-positive rate. Bootstrap 95% confidence intervals for overall rates are reported in Table 3.

## Supplementary Figure S2 — figure_supplement_infant_RI
**Title:** Reference interval bands across age groups: infants (30-365 days) versus the pediatric (1-18 y) intervals.
**Caption:** refineR reference interval estimates for PT (A), aPTT (B), and fibrinogen (C). Each panel shows age groups in ascending order on the x-axis: for PT and fibrinogen, infants (30-365 days) and pediatric (1-18 y); for aPTT, infants (30-365 days), pediatric 1-12 y partition, and pediatric 12-18 y partition. Vertical bars span the 2.5th-97.5th percentile interval; vertical whiskers represent 90% bootstrap confidence intervals at each reference limit (N = 200 for both pediatric and infant estimates). The infant group is shown in blue and pediatric groups in orange. Sample sizes: PT infant n = 1,367 / pediatric n = 17,000; aPTT infant n = 1,365 / 1-12 y n = 10,091 / 12-18 y n = 6,910; fibrinogen infant n = 217 / pediatric n = 3,166. Infant intervals differ markedly from the pediatric intervals for all three analytes — PT and aPTT upper limits are higher and fibrinogen limits are lower in infants, consistent with developmental hemostasis.

## Supplementary Tables — SUPPLEMENT_infant_30_365d.xlsx
**Sheet "Descriptive":** Pooled descriptive statistics (n, sex distribution, age in days, mean, SD, median, IQR, range, skewness, kurtosis) for infants aged 30-365 days; no outlier filter applied.
**Sheet "Classification_vs_RI":** Proportion of infant results falling Below, Within, or Above the pediatric (1-18 y for PT and fibrinogen; 1-12 y partition for aPTT) refineR reference interval. Used to demonstrate non-transferability of the pediatric RI to infants.
**Sheet "Infant_refineR_RI":** Lower (2.5th) and upper (97.5th) refineR reference limits with 90% bootstrap confidence intervals for each analyte in the infant subset (NBootstrap = 200, tolerance interval).
**Sheet "Pediatric_vs_Infant":** Side-by-side comparison table of pediatric (1-18 y; aPTT 1-12 y partition) and infant (30-365 days) refineR reference intervals, with 90% CIs for the infant limits.

## Supplementary Table — PUB_table_S_infant_pediatric_RI / RI_Table_Summary
**Title:** Reference intervals with sample sizes and 90% confidence intervals -- infants (30-365 days) and pediatric (1-18 y) cohorts.
**Caption:** Two-sided refineR reference intervals for PT, aPTT, and fibrinogen estimated separately for infants aged 30 to 365 days and pediatric cohorts aged 1 to 18 years. Lower (2.5th percentile) and upper (97.5th percentile) reference limits are reported with their 90% bootstrap confidence intervals in parentheses (NBootstrap = 200; tolerance intervals). For aPTT, the pediatric cohort is shown by the manuscript's Harris-Boyd partition (1-12 y and 12-18 y); PT and fibrinogen are pooled across 1-18 y. The infant cohort uses no outlier filter. The same data underlie Supplementary Figure S2.
**File:** `data/processed/tables/PUB_table_S_infant_pediatric_RI.csv` (standalone CSV) and `SUPPLEMENT_infant_30_365d.xlsx` (sheet `RI_Table_Summary`).
