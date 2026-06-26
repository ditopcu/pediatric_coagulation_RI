# ==============================================================================
# BUILD SUPPLEMENTARY METHODS WORD DOCUMENT
# Covers full pipeline parameters for all statistical methods used in the
# manuscript. Referenced in main text near "Supplementary Methods".
# Output: article_update/Supplementary_Methods.docx
# ==============================================================================

suppressPackageStartupMessages({
  library(officer)
})

out_path <- "article_update/Supplementary_Methods.docx"

# Convenience helpers ---------------------------------------------------------
H1 <- function(text) function(d) body_add_par(d, text, style = "heading 1")
H2 <- function(text) function(d) body_add_par(d, text, style = "heading 2")
P  <- function(text, size = 11) function(d) {
  body_add_fpar(d, fpar(ftext(text, prop = fp_text(font.size = size))))
}
BL <- function() function(d) body_add_par(d, "", style = "Normal")
BUL <- function(items, size = 10) function(d) {
  for (it in items) {
    d <- body_add_fpar(d, fpar(ftext(paste0("- ", it),
                                     prop = fp_text(font.size = size))))
  }
  d
}

# Content ---------------------------------------------------------------------
sections <- list(
  H1("Supplementary Methods"),
  P("This supplement provides the full parameter specifications, software versions, and procedural details for the statistical analyses summarised in the main Methods section. The same dataset and identical preprocessing were used throughout, except where explicitly noted."),
  BL(),

  H2("1. Data preparation and quality screening"),
  P("Coagulation measurements (PT, aPTT, fibrinogen) acquired between 2018 and 2024 on Roche Cobas t511 and t711 analysers were extracted from the laboratory information system. Patients aged 1 to less than 18 years at the time of testing were retained. Only the first valid measurement per analyte per patient was kept to avoid pseudo-replication."),
  P("No global outlier filter was applied for the indirect estimation, transferability, or false-positive rate analyses. A three-standard-deviation (3SD) filter (per analyte, applied to log-transformed values) was used only for the descriptive distribution figure (Figure 1) to improve visual clarity; this 3SD subset was not used to derive any reference interval or test statistic."),
  BL(),

  H2("2. Indirect reference intervals (refineR)"),
  P("Two-sided reference intervals were estimated using the refineR package (R, CRAN). Key parameters:"),
  BUL(c(
    "NBootstrap = 200 (bootstrap iterations for the 90% tolerance confidence intervals around each reference limit).",
    "RIperc = c(0.025, 0.50, 0.975) (lower 2.5th, median, upper 97.5th percentiles).",
    "CIprop = 0.90, pi.type = \"tolerance\" (90% tolerance interval around each percentile, as defined in refineR).",
    "Random seed = 42 for reproducibility.",
    "For aPTT, age partitioning at 12 years was applied prior to refineR estimation (justified by the Harris-Boyd test below). For PT and fibrinogen, a single 1-18 y pooled interval was estimated. For the supplementary infant analysis, refineR was run separately on the 30-365 day subset with the same parameters."
  )),
  BL(),

  H2("3. Age and sex partitioning (Harris-Boyd test)"),
  P("Partitioning decisions were based on the Harris-Boyd standardised z-statistic comparing two adjacent groups:"),
  P("z = |mean1 - mean2| / sqrt(sd1^2/n1 + sd2^2/n2)"),
  P("A threshold of z > 3 was used as the conventional cut-off for partition justification. Sex-specific partitioning was assessed separately for the 2.5th and 97.5th percentile estimates by inspecting overlap of the corresponding 90% confidence intervals across sexes within each integer age (Supplementary Figure 2 / Supplementary Table S3); CIs overlapped at the great majority of ages for both PT and aPTT, and a combined-sex RI was therefore reported. For age, only aPTT exceeded the z > 3 threshold at a candidate cut-point of 12 years; PT and fibrinogen did not, and were reported as single intervals."),
  BL(),

  H2("4. Direct (a posteriori) reference intervals"),
  P("A direct comparator cohort was defined a posteriori by restricting the dataset to outpatient samples from clinics consistent with apparently healthy referral (well-child / preoperative evaluation pathways) and excluding all inpatient, emergency department, and chronic-disease subspecialty visits. Within this subset, lower and upper reference limits were derived using the non-parametric percentile method per CLSI EP28-A3c (2.5th and 97.5th percentiles with 90% bootstrap confidence intervals, B = 1000)."),
  BL(),

  H2("5. Continuous reference intervals (GAMLSS pipeline)"),
  P("Continuous, age-dependent reference intervals were estimated using the GAMLSS-based indirect pipeline of Ammer et al. (algoRICurves.R), implemented in R. Per-analyte runner scripts (RUNNER_Pipeline_{PT,aPTT,Fib}.R) were used with the following parameters:"),
  BUL(c(
    "Input format: tab-delimited file with columns Age (integer years), Value (analyte result), PID (sample identifier).",
    "Covariate: Age (covarName = \"Age\"), modelled with a Box-Cox power transformation; the transformation parameter (pp) was selected automatically by the pipeline and propagated to .GlobalEnv as ptrans(x, p) = if (p == 0) log(x) else x^p.",
    "Candidate distribution families: BCCG, BCCGo, BCT, BCTo, BCPEo, LOGNO; final family selected by the pipeline's model-search routine.",
    "Smoothing: penalised B-splines (pb()) applied to the mu, sigma, and nu parameters of the selected family as functions of ptrans(Age, pp); tau (where applicable) was held constant.",
    "Sample weighting: probability of non-pathological membership (probNP) computed per observation by the pipeline's KDE-based outlier-attenuation routine; weights were passed to gamlss() via the weights argument.",
    "Point-estimate model: stored as <Test>_gamlss_PointEst.RData in cont_out/.",
    "Bootstrap confidence intervals: NBootstrap = 5 iterations (default for figures); each bootstrap fit stored as <Test>_gamlssModel_Est_{1..5}.RData.",
    "Centile prediction: centiles.pred() called with the seven percentiles c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975) at each integer age 1-18. The 2.5th, 50th, and 97.5th percentile curves with 90% bootstrap CIs are reported (Figure 4; Supplementary Table S4).",
    "Reproducibility: pipeline source files in src/pipeline/ are unchanged from the Ammer reference implementation; per-analyte runner scripts in src/ document the call parameters."
  )),
  BL(),

  H2("6. Transferability of the Guven et al. reference intervals"),
  P("Transferability of the published per-age reference intervals of Guven et al. (2026; Cobas t511) was assessed within each integer age in our dataset using a bootstrap analogue of the CLSI EP28-A3c verification protocol:"),
  BUL(c(
    "Per age and analyte, B = 200 bootstrap samples of size 20 were drawn with replacement from the local data (slice_sample(replace = TRUE)).",
    "For each bootstrap sample the number of values falling outside the Guven et al. reference limits for that age was counted; the sample was classified as \"pass\" if at most two of twenty values were outside (the EP28-A3c decision rule).",
    "The pass rate per age was the proportion of bootstrap iterations meeting the pass criterion. A pass rate of >= 90% was used as the threshold for transferability.",
    "In addition, an exact one-sided binomial test (binom.test, alternative = \"greater\", p0 = 0.05) was applied to the full local data at each age to test whether the observed proportion of out-of-RI results significantly exceeded the nominal 5% rate expected under a correctly transferred interval. Departures from this null hypothesis with p < 0.05 were reported as binomial rejections.",
    "Random seed = 42."
  )),
  BL(),

  H2("7. Bootstrap confidence intervals for out-of-RI flagging rates"),
  P("For Table 3 and Figure 6, 95% bootstrap confidence intervals around the overall (age-pooled) out-of-RI flagging rate were generated as follows:"),
  BUL(c(
    "For each test x RI source (manufacturer, refineR, continuous GAMLSS, Guven et al.), B = 1000 bootstrap resamples of size n (with replacement, n = number of observations for that test) were drawn from the per-observation out-of-RI indicator.",
    "The 95% CI was derived from the 2.5th and 97.5th percentiles of the bootstrap distribution.",
    "Random seed = 123."
  )),
  BL(),

  H2("8. Software"),
  BUL(c(
    "R 4.4.x / 4.5.x on Windows.",
    "Core packages: tidyverse 2.0.0, readxl, refineR (CRAN), gamlss 5.5-0 (with gamlss.dist, gamlss.data, nlme), patchwork, ggrepel, cowplot, janitor, moments, openxlsx, officer, flextable.",
    "Ammer GAMLSS pipeline source: src/pipeline/ (RUNNER_Pipeline.R, algoRICurves.R, plotRICurve.R, utils.R), unmodified.",
    "All analysis scripts are in src/, with src/run_all.R as the master runner. Per-analyte refineR fits are persisted as .RDS files under data/processed/ and GAMLSS model objects under cont_out/ for full reproducibility."
  )),
  BL(),

  P("End of Supplementary Methods.", size = 9)
)

# Render ----------------------------------------------------------------------
doc <- read_docx()
for (f in sections) doc <- f(doc)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
print(doc, target = out_path)
cat(sprintf("[OK] Supplementary Methods saved: %s\n", out_path))
