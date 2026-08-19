# ==============================================================================
# ANALYSIS (REVIEWER RESPONSE): Distribution shape of raw results
# Pediatric Coagulation Reference Intervals
# ==============================================================================
# Reviewer request: "distribution curves of the results" + comment on PT
# right-skew. This standalone script produces:
#   1. figures/TIFF_600DPI/figure_supplement_distribution_REV.tiff
#      figures/PNG_300DPI/figure_supplement_distribution_REV.png
#      Three stacked panels (A: PT, B: aPTT, C: Fibrinogen):
#      raw-data histogram (density scale) + refineR-estimated
#      non-pathological density (solid orange) + refineR RI limits (dashed).
#   2. data/processed/REVIEW_distribution_shape_REV.csv
#      Raw-data shape statistics per analyte.
#
# Run from project root:  Rscript src/analysis_distribution_REV.R
#
# Method notes (no assumptions -- verified against installed refineR 2.0.0):
# - RDS objects are class RWDRI (refineR). Fitted non-pathological model is a
#   Box-Cox transformed normal: z = BoxCox(x - Shift, Lambda) ~ N(Mu, Sigma),
#   with P = estimated non-pathological fraction of the input data.
# - refineR:::plot.RWDRI draws the estimated curve as expected bin counts:
#     countsPred = pCorr * N * P *
#                  (pnorm(BoxCox(breakR - Shift, Lambda), Mu, Sigma) -
#                   pnorm(BoxCox(breakL - Shift, Lambda), Mu, Sigma))
#   where pCorr renormalizes the fitted normal over the observed data range.
#   The continuous density-scale equivalent used here is:
#     f(x) = P * pCorr * dnorm(BoxCox(x - Shift, Lambda), Mu, Sigma)
#                      * (x - Shift)^(Lambda - 1)
#   ((x - Shift)^(Lambda - 1) is the Jacobian of the Box-Cox transform; it
#   equals 1/(x - Shift) for Lambda = 0, matching BoxCox(x, 0) = log(x)).
# - Overlay scaling: the histogram is on the density scale of ALL raw data
#   (integrates to 1), and the model curve is scaled by the object's own
#   estimated non-pathological fraction P (x$P), exactly as plot.RWDRI scales
#   its countsPred by N * P. No ad-hoc rescaling is applied.
# - VERIFICATION GATE: the 2.5th / 97.5th percentiles of the computed density
#   (numerical CDF inversion) must reproduce getRI(obj)$PointEst within 1%
#   per analyte; the curve is only drawn for analytes that PASS. getRI()
#   internally uses refineR:::cdfTruncatedBoxCox (truncation at the Box-Cox
#   support bound); our density is additionally truncated to the observed data
#   range like plot.RWDRI (pCorr), a negligible difference at these tails.
# ==============================================================================

source("src/common.R")

library(refineR)
library(moments)

message("[REV-DIST] refineR version: ", as.character(packageVersion("refineR")))

# ------------------------------------------------------------------------------
# 1. Load refineR fit objects (RWDRI)
# ------------------------------------------------------------------------------

fit_pt   <- readRDS("data/processed/ptz_ref_1_18.RDS")    # PT, 1-18 y
fit_aptt <- readRDS("data/processed/aptt_ref_1_18.RDS")   # aPTT, pooled 1-18 y
fit_fib  <- readRDS("data/processed/fib_ref_1_18.RDS")    # Fibrinogen, 1-18 y

stopifnot(inherits(fit_pt, "RWDRI"), inherits(fit_aptt, "RWDRI"),
          inherits(fit_fib, "RWDRI"))

# ------------------------------------------------------------------------------
# 2. Non-pathological density (replicates refineR:::plot.RWDRI computation)
# ------------------------------------------------------------------------------

# Density of the estimated non-pathological subpopulation on the original
# scale, scaled by the estimated non-pathological fraction P so that it is
# directly comparable to a density-scale histogram of ALL raw data.
np_density <- function(obj, x) {
  lambda <- obj$Lambda; mu <- obj$Mu; sigma <- obj$Sigma
  shift  <- obj$Shift;  P  <- obj$P
  # Truncation correction over the observed data range, exactly as plot.RWDRI:
  pC <- suppressWarnings(refineR::BoxCox(
    c(max(min(obj$Data - shift), 1e-20), min(max(obj$Data - shift), 1e+20)),
    lambda = lambda))
  pC    <- pnorm(pC, mean = mu, sd = sigma)
  pCorr <- 1 / (pC[2] - pC[1])
  z <- suppressWarnings(refineR::BoxCox(x - shift, lambda = lambda))
  d <- P * pCorr * dnorm(z, mean = mu, sd = sigma) * (x - shift)^(lambda - 1)
  d[!is.finite(d) | d < 0] <- 0
  d
}

# Numerical 2.5th / 97.5th percentiles of the computed density
# (trapezoid CDF over a fine grid spanning the observed data range).
np_quantiles <- function(obj, probs = c(0.025, 0.975), n_grid = 40001) {
  lo <- max(min(obj$Data), obj$Shift + 1e-9)
  hi <- max(obj$Data)
  x  <- seq(lo, hi, length.out = n_grid)
  d  <- np_density(obj, x)
  dx  <- diff(x)
  cdf <- c(0, cumsum((d[-1] + d[-length(d)]) / 2 * dx))
  cdf <- cdf / cdf[length(cdf)]
  keep <- c(TRUE, diff(cdf) > 0)   # drop flat tail segments for inversion
  approx(x = cdf[keep], y = x[keep], xout = probs)$y
}

# ------------------------------------------------------------------------------
# 3. Verification gate (printed PASS/FAIL per analyte)
# ------------------------------------------------------------------------------

verify_fit <- function(obj, label, tol = 0.01) {
  ri  <- getRI(obj)                       # defaults: 2.5th and 97.5th
  q   <- np_quantiles(obj)
  rel <- abs(q - ri$PointEst) / ri$PointEst
  ok  <- all(rel <= tol)
  cat(sprintf(paste0("[GATE] %-10s q2.5 = %8.4f vs getRI %8.4f (%5.3f%%) | ",
                     "q97.5 = %8.4f vs getRI %8.4f (%5.3f%%) -> %s\n"),
              label, q[1], ri$PointEst[1], 100 * rel[1],
              q[2], ri$PointEst[2], 100 * rel[2],
              ifelse(ok, "PASS", "FAIL")))
  ok
}

cat("\n[REV-DIST] Verification gate: computed density percentiles vs getRI()\n")
gate <- c(
  PT         = verify_fit(fit_pt,   "PT"),
  aPTT       = verify_fit(fit_aptt, "aPTT"),
  Fibrinogen = verify_fit(fit_fib,  "Fibrinogen")
)
if (!all(gate)) {
  cat("[REV-DIST] WARNING: gate FAILED for:",
      paste(names(gate)[!gate], collapse = ", "),
      "- falling back to histogram + RI limit lines only for those panels.\n")
}

# ------------------------------------------------------------------------------
# 4. Raw-data shape statistics CSV
# ------------------------------------------------------------------------------
# Study refineR limits (from ri_refiner in common.R):
#   PT 8.56-10.90 (1-18 y); Fibrinogen 1.89-3.79 (1-18 y);
#   aPTT partitioned: 24.30-35.30 (1-12 y), 23.30-33.90 (12-18 y).
# For aPTT the age-appropriate partition limits are applied per row
# (age_year < 12 vs >= 12) before aggregating.

raw <- all_data |> filter(!is.na(result_num))
n_na <- nrow(all_data) - nrow(raw)
cat(sprintf("\n[REV-DIST] Raw rows: %d (%d NA result_num excluded)\n",
            nrow(raw), n_na))

shape_stats <- raw |>
  mutate(
    lower_lim = case_when(
      test == "PT"                        ~  8.56,
      test == "Fibrinogen"                ~  1.89,
      test == "aPTT" & age_year <  12     ~ 24.30,
      test == "aPTT" & age_year >= 12     ~ 23.30
    ),
    upper_lim = case_when(
      test == "PT"                        ~ 10.90,
      test == "Fibrinogen"                ~  3.79,
      test == "aPTT" & age_year <  12     ~ 35.30,
      test == "aPTT" & age_year >= 12     ~ 33.90
    )
  ) |>
  group_by(test) |>
  summarise(
    n        = n(),
    mean     = round(mean(result_num), 3),
    sd       = round(sd(result_num), 3),
    median   = round(median(result_num), 3),
    skewness = round(moments::skewness(result_num), 3),
    kurtosis = round(moments::kurtosis(result_num), 3),  # Pearson; normal = 3
    pct_below_lower_refineR_limit = round(100 * mean(result_num < lower_lim), 2),
    pct_above_upper_refineR_limit = round(100 * mean(result_num > upper_lim), 2),
    .groups = "drop"
  ) |>
  arrange(match(test, c("PT", "aPTT", "Fibrinogen")))

write_csv(shape_stats, "data/processed/REVIEW_distribution_shape_REV.csv")
cat("\n[REV-DIST] Shape statistics (raw data):\n")
print(as.data.frame(shape_stats), row.names = FALSE)
cat("[OK] Written: data/processed/REVIEW_distribution_shape_REV.csv\n\n")

# ------------------------------------------------------------------------------
# 5. Figure: three stacked panels (A: PT, B: aPTT, C: Fibrinogen)
# ------------------------------------------------------------------------------
# X-axis cap: upper bound at the 99.9th percentile of the raw data per analyte
# (keeps the pathological right tail visible while the bulk stays readable);
# lower bound at the 0.1th percentile (avoids stretching the axis to isolated
# extreme low values, e.g. the aPTT entry-error value 0.000477).
# Histogram bins are computed over the full raw data with binwidth =
# display range / 60, so the density scale remains normalized to ALL raw data.

make_dist_panel <- function(test_name, obj, xlab_txt, tag,
                            ri_lines, draw_curve, aptt_ann = NULL) {
  v  <- raw$result_num[raw$test == test_name]
  lo <- as.numeric(quantile(v, 0.001))
  hi <- as.numeric(quantile(v, 0.999))   # 99.9th percentile x-axis cap
  bw <- (hi - lo) / 60

  sk <- moments::skewness(v)
  ku <- moments::kurtosis(v)

  p <- ggplot(data.frame(x = v), aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = bw, boundary = lo,
                   fill = PAL$base2, alpha = 0.55,
                   colour = "grey55", linewidth = 0.15) +
    geom_vline(xintercept = ri_lines, colour = PAL$refiner,
               linetype = "dashed", linewidth = 0.4)

  if (draw_curve) {
    curve_df <- data.frame(x = seq(lo, hi, length.out = 1500))
    curve_df$y <- np_density(obj, curve_df$x)
    p <- p + geom_line(data = curve_df, aes(x = x, y = y),
                       colour = PAL$refiner, linewidth = 0.8)
  }

  if (!is.null(aptt_ann)) {
    p <- p + geom_text(data = aptt_ann,
                       aes(x = x, y = Inf, label = lab, vjust = vjs),
                       angle = 90, hjust = 1.05, size = 2.2,
                       colour = PAL$refiner, inherit.aes = FALSE)
  }

  p +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4,
             size = 2.7, colour = PAL$text, lineheight = 1.15,
             label = sprintf("skewness = %.2f\nkurtosis = %.2f", sk, ku)) +
    coord_cartesian(xlim = c(lo, hi)) +
    labs(x = xlab_txt, y = "Density", tag = tag) +
    theme_tufte_academic() +
    theme(plot.tag = element_text(size = 12, face = "bold"))
}

# aPTT partition-limit annotations (subtle, along the dashed lines):
# labels flank outward from each close pair of lines so they stay legible.
# Note: for angle = 90 text, vjust > 1 offsets the label to the RIGHT of the
# line and vjust < 0 to the LEFT.
aptt_ann <- data.frame(
  x   = c(24.30, 23.30, 35.30, 33.90),
  lab = c("24.3 (1-12 y)", "23.3 (12-18 y)", "35.3 (1-12 y)", "33.9 (12-18 y)"),
  vjs = c(1.35, -0.35, 1.35, -0.35)
)

p_pt <- make_dist_panel("PT", fit_pt, "PT, seconds", "A",
                        ri_lines = c(8.56, 10.90),
                        draw_curve = gate[["PT"]])
p_aptt <- make_dist_panel("aPTT", fit_aptt, "aPTT, seconds", "B",
                          ri_lines = c(24.30, 35.30, 23.30, 33.90),
                          draw_curve = gate[["aPTT"]],
                          aptt_ann = aptt_ann)
p_fib <- make_dist_panel("Fibrinogen", fit_fib, "Fibrinogen, g/L", "C",
                         ri_lines = c(1.89, 3.79),
                         draw_curve = gate[["Fibrinogen"]])

fig <- p_pt / p_aptt / p_fib

save_fig(fig, "figure_supplement_distribution_REV", width = 180, height = 235)

cat("[REV-DIST] Done.\n")
