# ==============================================================================
# PUBLICATION FIGURES
# Pediatric Coagulation Reference Intervals
# ==============================================================================
# Generates all publication figures using Tufte/Annesley design.
# Requires: source("src/common.R") to be run first (via run_all.R or manually).
# Output: figures/TIFF_600DPI/ and figures/PNG_300DPI/
# Captions: docs/figure_captions.md (not embedded in figures)
# ==============================================================================

if (!exists("PAL")) source("src/common.R")

library(refineR)
library(ggrepel)

# pub.R needs ordered factor for boxplot x-axis
all_data_3sd <- all_data_3sd |> mutate(age_int = ordered(age_int))

# Additional data for pub.R only
pt_ref_sex  <- tryCatch(readRDS("data/processed/ptz_ref_by_age_sex_raw.RDS"),  error = function(e) NULL)
aptt_ref_sex <- tryCatch(readRDS("data/processed/aptt_ref_by_age_sex_raw.RDS"), error = function(e) NULL)

# Alias for backward compatibility
manufacturer_ri <- ri_manufacturer

ri_long <- ri_comparison |>
  mutate(panel_label = paste0(test, "  ", age_group)) |>
  mutate(panel_label = fct_inorder(panel_label))


# ==============================================================================
# FIGURE 1 — Distribution (Boxplot + Jitter)
# ==============================================================================

message("[PLOT] Figure 1: Distribution ...")

make_distribution_panel <- function(df, test_name, y_lab, tag_label) {
  df |>
    filter(test == test_name) |>
    ggplot(aes(x = age_int, y = result_num, colour = sex)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                               dodge.width = 0.6),
               alpha = 0.25, size = 0.6, shape = 16) +
    geom_boxplot(aes(fill = sex), position = position_dodge(width = 0.6),
                 width = 0.45, outlier.shape = NA, alpha = 0.15,
                 linewidth = 0.4, colour = PAL$text) +
    scale_colour_manual(values = c("K" = PAL$female, "E" = PAL$male),
                        labels = c("K" = "Female", "E" = "Male")) +
    scale_fill_manual(values = c("K" = PAL$female, "E" = PAL$male),
                      labels = c("K" = "Female", "E" = "Male")) +
    labs(x = "Age, years", y = y_lab) +
    add_panel_tag(tag_label) +
    theme_tufte_academic() +
    theme(legend.position = "none")
}

p_dist_pt   <- make_distribution_panel(all_data_3sd, "PT",         "PT, seconds",     "A")
p_dist_aptt <- make_distribution_panel(all_data_3sd, "aPTT",       "aPTT, seconds",   "B")
p_dist_fib  <- make_distribution_panel(all_data_3sd, "Fibrinogen", "Fibrinogen, g/L", "C") +
  theme(legend.position = "bottom")

fig_1 <- p_dist_pt / p_dist_aptt / p_dist_fib
save_fig(fig_1, "figure_1_distribution", width = 180, height = 260)


# ==============================================================================
# FIGURE 2 — Sex-specific partitioning CI (PT and aPTT)
# ==============================================================================

message("[PLOT] Figure 2: Partitioning ...")

make_partition_plot <- function(rds_obj, test_label, y_lab, tag_label,
                               overall_ri_list) {
  if (is.null(rds_obj)) {
    message("  [WARN] RDS not found: ", test_label)
    return(ggplot() + annotate("text", x = 1, y = 1, label = "Data not available") +
             theme_void())
  }

  part_data <- rds_obj |>
    mutate(ric = map(ri, function(x)
      getRI(x, RIperc = c(0.025, 0.975), CIprop = 0.90))) |>
    unnest(ric) |>
    select(age_tam_sayi, sex, n, Percentile, PointEst, CILow, CIHigh)

  # Reverse facet order: Upper RL on top, Lower RL on bottom
  part_data <- part_data |>
    mutate(Percentile = factor(Percentile, levels = c("0.975", "0.025")))

  # Build overall RI background bands per facet
  # overall_ri_list: list of tibbles with columns: age_start, age_end, percentile, point_est, ci_lo, ci_hi
  band_data <- overall_ri_list

  # Ensure age is numeric for positioning
  part_data <- part_data |>
    mutate(age_tam_sayi = as.numeric(as.character(age_tam_sayi)))

  # "All" column: overall RI point + CI at x = max_age + 1
  max_age <- max(part_data$age_tam_sayi, na.rm = TRUE)
  all_x <- max_age + 1

  all_points <- band_data |>
    group_by(Percentile) |>
    summarise(
      PointEst = mean(point_est),
      CILow = mean(ci_lo),
      CIHigh = mean(ci_hi),
      .groups = "drop"
    ) |>
    mutate(age_tam_sayi = all_x, sex = "All")

  p <- ggplot(part_data, aes(x = age_tam_sayi)) +
    # Background: overall RI bands (behind everything)
    geom_rect(data = band_data,
              aes(xmin = age_start - 0.5, xmax = age_end + 0.5,
                  ymin = ci_lo, ymax = ci_hi),
              fill = PAL$refiner, alpha = 0.15, inherit.aes = FALSE) +
    geom_hline(data = band_data,
               aes(yintercept = point_est),
               colour = PAL$refiner, linewidth = 0.5, alpha = 0.7, linetype = "solid") +
    # Age-specific points and CIs
    geom_errorbar(aes(ymin = CILow, ymax = CIHigh, colour = sex),
                  linewidth = 0.6, width = 0.3,
                  position = position_dodge(width = 0.5)) +
    geom_point(aes(y = PointEst, colour = sex), size = 1.8,
               position = position_dodge(width = 0.5)) +
    # "All" point at right edge
    geom_errorbar(data = all_points,
                  aes(x = age_tam_sayi, ymin = CILow, ymax = CIHigh),
                  linewidth = 0.7, width = 0.3, colour = PAL$refiner) +
    geom_point(data = all_points,
               aes(x = age_tam_sayi, y = PointEst),
               size = 2.2, colour = PAL$refiner, shape = 18) +
    facet_wrap(~Percentile, ncol = 1, scales = "free_y",
               labeller = labeller(Percentile = c("0.975" = "Upper RL (97.5th)",
                                                  "0.025" = "Lower RL (2.5th)"))) +
    scale_colour_manual(values = c("K" = PAL$female, "E" = PAL$male),
                        labels = c("K" = "Female", "E" = "Male"),
                        guide = guide_legend(order = 1)) +
    scale_x_continuous(breaks = c(seq(1, max_age), all_x),
                       labels = c(as.character(seq(1, max_age)), "All")) +
    labs(x = "Age, years", y = paste0(test_label, ", ", y_lab)) +
    add_panel_tag(tag_label) +
    theme_tufte_academic() +
    theme(legend.position = "bottom",
          panel.spacing = unit(15, "pt"))

  return(p)
}

# Overall RI bands for PT (single group 1-18)
pt_bands <- tribble(
  ~Percentile, ~age_start, ~age_end, ~point_est, ~ci_lo, ~ci_hi,
  "0.025",     1,          17,       8.56,       8.54,   8.58,
  "0.975",     1,          17,       10.90,      10.70,  10.90
) |> mutate(Percentile = factor(Percentile, levels = c("0.975", "0.025")))

# Overall RI bands for aPTT (two groups: 1-12 and 12-18)
aptt_bands <- tribble(
  ~Percentile, ~age_start, ~age_end, ~point_est, ~ci_lo, ~ci_hi,
  "0.025",     1,          11,       24.30,      23.30,  24.50,
  "0.025",     12,         17,       23.30,      22.70,  23.70,
  "0.975",     1,          11,       35.30,      34.40,  35.40,
  "0.975",     12,         17,       33.90,      33.10,  34.10
) |> mutate(Percentile = factor(Percentile, levels = c("0.975", "0.025")))

fig_2 <- make_partition_plot(pt_ref_sex, "PT", "seconds", "A", pt_bands) |
  make_partition_plot(aptt_ref_sex, "aPTT", "seconds", "B", aptt_bands)
save_fig(fig_2, "figure_2_partitioning", width = 260, height = 180)


# ==============================================================================
# FIGURE 3 — Indirect vs Direct RI comparison
# ==============================================================================

message("[PLOT] Figure 3: Indirect vs Direct ...")

ri_with_mfr <- ri_long |>
  left_join(manufacturer_ri, by = "test")

fig_3 <- ri_with_mfr |>
  ggplot(aes(y = fct_rev(method), fill = method)) +
  geom_crossbar(aes(x = (lower.x + upper.x)/2, xmin = lower.x, xmax = upper.x),
                width = 0.4, linewidth = 0.4, alpha = 0.75, colour = NA) +
  geom_errorbarh(aes(xmin = lower_ci_lo, xmax = lower_ci_hi),
                 height = 0.2, linewidth = 0.35, colour = "#444444") +
  geom_errorbarh(aes(xmin = upper_ci_lo, xmax = upper_ci_hi),
                 height = 0.2, linewidth = 0.35, colour = "#444444") +
  geom_vline(aes(xintercept = lower.y), linetype = "dashed",
             colour = PAL$base2, linewidth = 0.5) +
  geom_vline(aes(xintercept = upper.y), linetype = "dashed",
             colour = PAL$base2, linewidth = 0.5) +
  facet_wrap(~panel_label, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = c("refineR" = PAL$refiner, "Direct" = PAL$direct),
                    labels = c("refineR" = "Indirect (refineR)",
                               "Direct"  = "Direct (a posteriori)")) +
  labs(x = "", y = "") +
  theme_tufte_academic() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(10, "pt"))

save_fig(fig_3, "figure_3_indirect_vs_direct", width = 140, height = 200)


# ==============================================================================
# FIGURE 4 — Continuous RI curves (GAMLSS)
# ==============================================================================

message("[PLOT] Figure 4: Continuous RI ...")

plot_continuous_ri <- function(point_est_file, bs_files, input_data,
                               test_label, y_lab, tag_label,
                               covar_name = "Age", value_name = "Value") {

  pe_exists <- file.exists(point_est_file)

  if (!pe_exists) {
    message("  [WARN] Pipeline output not found: ", point_est_file)
    p <- input_data |>
      ggplot(aes(x = as.numeric(as.character(age_int)), y = result_num)) +
      geom_point(alpha = 0.05, size = 0.3, colour = PAL$base2) +
      geom_smooth(method = "loess", span = 0.5, se = TRUE,
                  colour = PAL$contin, fill = PAL$contin, alpha = 0.2) +
      labs(x = "Age, years", y = paste0(test_label, ", ", y_lab)) +
      add_panel_tag(tag_label) +
      theme_tufte_academic()
    return(p)
  }

  load(point_est_file)

  pp <- pointEstGamlss$Params$transformP
  assign("pp", pp, envir = .GlobalEnv)
  assign("ptrans", function(x, p) if (p == 0) log(x) else I(x^p), envir = .GlobalEnv)
  assign("families", pointEstGamlss$Models$family, envir = .GlobalEnv)
  assign("fam", 1, envir = .GlobalEnv)

  bs_results <- lapply(bs_files, function(f) {
    if (file.exists(f)) { load(f); return(gamlssModel) }
    else return(NULL)
  })
  bs_results <- Filter(Negate(is.null), bs_results)

  source("src/pipeline/utils.R")
  source("src/pipeline/plotRICurve.R")
  assign("covarName", covar_name, envir = .GlobalEnv)

  covar_vals <- seq(1, 18, by = 1)

  # Use all 7 default percentiles so CILow/CIHigh indices match RIMat rows
  full_perc <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
  allRes <- tryCatch(
    estimateCIs(pointEst = pointEstGamlss, estBS = bs_results,
                covarValue = covar_vals, withData = TRUE,
                RIperc = full_perc, CIprop = 0.90),
    error = function(e) { message("  [ERROR] estimateCIs: ", e$message); NULL }
  )

  if (is.null(allRes)) {
    p <- input_data |>
      ggplot(aes(x = as.numeric(as.character(age_int)), y = result_num)) +
      geom_point(alpha = 0.05, size = 0.3, colour = PAL$base2) +
      geom_smooth(method = "loess", span = 0.5, se = TRUE,
                  colour = PAL$contin, fill = PAL$contin, alpha = 0.2) +
      labs(x = "Age, years", y = paste0(test_label, ", ", y_lab)) +
      add_panel_tag(tag_label) +
      theme_tufte_academic()
    return(p)
  }

  ri_mat <- allRes$RICurve$RIMat
  # CILow/CIHigh have 7 elements matching full_perc: 1=0.025, 4=0.5, 7=0.975
  ri_df <- tibble(
    age     = covar_vals,
    p2.5    = ri_mat["Perc_0.025", ],
    p50     = ri_mat["Perc_0.5",   ],
    p97.5   = ri_mat["Perc_0.975", ],
    ci_lo_lower  = allRes$CILow[[1]],
    ci_hi_lower  = allRes$CIHigh[[1]],
    ci_lo_median = allRes$CILow[[4]],
    ci_hi_median = allRes$CIHigh[[4]],
    ci_lo_upper  = allRes$CILow[[7]],
    ci_hi_upper  = allRes$CIHigh[[7]]
  )

  # Colors: all continuous RI elements use PAL$contin (green)
  p <- ggplot(ri_df, aes(x = age)) +
    # CI bands (light green)
    geom_ribbon(aes(ymin = ci_lo_upper,  ymax = ci_hi_upper),
                fill = PAL$contin, alpha = 0.12) +
    geom_ribbon(aes(ymin = ci_lo_lower,  ymax = ci_hi_lower),
                fill = PAL$contin, alpha = 0.12) +
    geom_ribbon(aes(ymin = ci_lo_median, ymax = ci_hi_median),
                fill = PAL$contin, alpha = 0.10) +
    # Percentile curves (green, median thicker)
    geom_line(aes(y = p97.5),  colour = PAL$contin, linewidth = 0.7) +
    geom_line(aes(y = p2.5),   colour = PAL$contin, linewidth = 0.7) +
    geom_line(aes(y = p50),    colour = PAL$contin, linewidth = 1.2) +
    # Direct labels (Annesley recommendation)
    annotate("text", x = 17.2, y = last(ri_df$p97.5), label = "97.5th",
             colour = PAL$contin, size = 3, hjust = 0, fontface = "bold") +
    annotate("text", x = 17.2, y = last(ri_df$p50),   label = "50th",
             colour = PAL$contin, size = 3, hjust = 0, fontface = "bold") +
    annotate("text", x = 17.2, y = last(ri_df$p2.5),  label = "2.5th",
             colour = PAL$contin, size = 3, hjust = 0, fontface = "bold") +
    scale_x_continuous(breaks = 1:18) +
    coord_cartesian(xlim = c(1, 19), clip = "off") +
    labs(x = "Age, years", y = paste0(test_label, ", ", y_lab)) +
    add_panel_tag(tag_label) +
    theme_tufte_academic() +
    theme(legend.position = "none")

  return(p)
}

fig_4_pt <- plot_continuous_ri(
  point_est_file = "cont_out/PT_gamlss_PointEst.RData",
  bs_files       = paste0("cont_out/PT_gamlssModel_Est_", 1:5, ".RData"),
  input_data     = filter(all_data_3sd, test == "PT"),
  test_label     = "PT", y_lab = "seconds", tag_label = "A"
)

fig_4_aptt <- plot_continuous_ri(
  point_est_file = "cont_out/aPTT_gamlss_PointEst.RData",
  bs_files       = paste0("cont_out/aPTT_gamlssModel_Est_", 1:5, ".RData"),
  input_data     = filter(all_data_3sd, test == "aPTT"),
  test_label     = "aPTT", y_lab = "seconds", tag_label = "B"
)

fig_4_fib <- plot_continuous_ri(
  point_est_file = "cont_out/Fib_gamlss_PointEst.RData",
  bs_files       = paste0("cont_out/Fib_gamlssModel_Est_", 1:5, ".RData"),
  input_data     = filter(all_data_3sd, test == "Fibrinogen"),
  test_label     = "Fibrinogen", y_lab = "g/L", tag_label = "C"
)

fig_4 <- fig_4_pt / fig_4_aptt / fig_4_fib
save_fig(fig_4, "figure_4_continuous_RI", width = 180, height = 310)


# ==============================================================================
# FIGURE 5 — Three methods combined
# ==============================================================================

message("[PLOT] Figure 5: Combined three methods ...")

make_combined_plot <- function(continuous_plot, ri_comp_subset,
                               mfr_lower, mfr_upper,
                               test_label, y_lab, tag_label) {

  ri_bars <- ri_comp_subset |>
    mutate(
      age_start = case_when(
        age_group == "1-18 y"  ~ 1,
        age_group == "1-12 y"  ~ 1,
        age_group == "12-18 y" ~ 12,
        TRUE ~ 1
      ),
      age_end = case_when(
        age_group == "1-18 y"  ~ 18,
        age_group == "1-12 y"  ~ 12,
        age_group == "12-18 y" ~ 18,
        TRUE ~ 18
      )
    )

  p <- continuous_plot +
    # Manufacturer RI (gray dotted horizontal lines)
    geom_hline(yintercept = mfr_upper, linetype = "dotted",
               colour = PAL$base2, linewidth = 0.6) +
    geom_hline(yintercept = mfr_lower, linetype = "dotted",
               colour = PAL$base2, linewidth = 0.6) +
    # refineR segments (orange solid)
    geom_segment(data = filter(ri_bars, method == "refineR"),
                 aes(x = age_start, xend = age_end,
                     y = upper, yend = upper),
                 colour = PAL$refiner, linewidth = 0.8, linetype = "solid") +
    geom_segment(data = filter(ri_bars, method == "refineR"),
                 aes(x = age_start, xend = age_end,
                     y = lower, yend = lower),
                 colour = PAL$refiner, linewidth = 0.8, linetype = "solid") +
    # Direct segments (blue longdash)
    geom_segment(data = filter(ri_bars, method == "Direct"),
                 aes(x = age_start, xend = age_end,
                     y = upper, yend = upper),
                 colour = PAL$direct, linewidth = 0.8, linetype = "longdash") +
    geom_segment(data = filter(ri_bars, method == "Direct"),
                 aes(x = age_start, xend = age_end,
                     y = lower, yend = lower),
                 colour = PAL$direct, linewidth = 0.8, linetype = "longdash") +
    labs(y = paste0(test_label, ", ", y_lab)) +
    add_panel_tag(tag_label)

  return(p)
}

fig_5_pt <- make_combined_plot(
  continuous_plot = fig_4_pt + theme(plot.tag = element_blank()),
  ri_comp_subset  = filter(ri_comparison, test == "PT"),
  mfr_lower = 8.4, mfr_upper = 10.6,
  test_label = "PT", y_lab = "seconds", tag_label = "A"
)

fig_5_aptt <- make_combined_plot(
  continuous_plot = fig_4_aptt + theme(plot.tag = element_blank()),
  ri_comp_subset  = filter(ri_comparison, test == "aPTT"),
  mfr_lower = 23.6, mfr_upper = 30.6,
  test_label = "aPTT", y_lab = "seconds", tag_label = "B"
)

fig_5_fib <- make_combined_plot(
  continuous_plot = fig_4_fib + theme(plot.tag = element_blank()),
  ri_comp_subset  = filter(ri_comparison, test == "Fibrinogen"),
  mfr_lower = 1.93, mfr_upper = 4.12,
  test_label = "Fibrinogen", y_lab = "g/L", tag_label = "C"
)

# Build shared legend for figure 5
legend_data <- tibble(
  x = 1:4, y = 1:4,
  method = factor(c("Continuous (GAMLSS)", "Indirect (refineR)",
                     "Direct (a posteriori)", "Manufacturer"),
                  levels = c("Continuous (GAMLSS)", "Indirect (refineR)",
                             "Direct (a posteriori)", "Manufacturer"))
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, colour = method, linetype = method)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(
    "Continuous (GAMLSS)"    = PAL$contin,
    "Indirect (refineR)"     = PAL$refiner,
    "Direct (a posteriori)"  = PAL$direct,
    "Manufacturer"           = PAL$base2
  )) +
  scale_linetype_manual(values = c(
    "Continuous (GAMLSS)"    = "solid",
    "Indirect (refineR)"     = "solid",
    "Direct (a posteriori)"  = "longdash",
    "Manufacturer"           = "dotted"
  )) +
  theme_tufte_academic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 9))

fig_5_legend <- cowplot::get_plot_component(legend_plot, "guide-box-bottom", return_all = TRUE)

fig_5 <- (fig_5_pt / fig_5_aptt / fig_5_fib) /
  wrap_elements(fig_5_legend) +
  plot_layout(heights = c(1, 1, 1, 0.08))
save_fig(fig_5, "figure_5_combined_three_methods", width = 180, height = 300)


# ==============================================================================
# DONE
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("[DONE] All figures saved to 'figures/'.")
message("   TIFF (600 DPI): figures/TIFF_600DPI/")
message("   PNG  (300 DPI): figures/PNG_300DPI/")
message(paste(rep("=", 60), collapse = ""))
