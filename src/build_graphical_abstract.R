# ==============================================================================
# GRAPHICAL ABSTRACT
# ==============================================================================
# One 13 x 5 cm panel for journal submission. Layout: the study's own worked
# example on the left (a single aPTT result against each reference interval
# source), the out-of-RI flagging rates on the right.
#
# Elsevier asks for at least 531 x 1328 px (h x w), readable at 5 x 13 cm.
# The TIFF at 600 dpi is 1181 x 3071 px and the PNG at 300 dpi is 591 x 1535,
# both above that.
#
# Every number is read from the analysis outputs; none is written out here.
#
# Usage: setwd(project root) then source("src/build_graphical_abstract.R")
# ==============================================================================

if (!exists("PAL")) source("src/common.R")

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(readr); library(patchwork); library(tidyr)
})

message("[BUILD] Graphical abstract ...")

# --- The worked example -------------------------------------------------------
# aPTT 32 s at age 8: the case the Discussion uses.
GA_TEST  <- "aPTT"
GA_AGE   <- 8
GA_VALUE <- 32

curves_ga <- read_csv("data/processed/tables/PUB_table_S3_continuous_RI_values.csv",
                      show_col_types = FALSE)

src_limits <- bind_rows(
  ri_manufacturer |> filter(test == GA_TEST) |>
    transmute(source = "Manufacturer", lower, upper),
  ri_refiner |> filter(test == GA_TEST, age_group == "1-12") |>
    transmute(source = "Indirect (refineR)", lower, upper),
  ri_direct |> filter(test == GA_TEST, age_group == "1-12") |>
    transmute(source = "Direct (a posteriori)", lower, upper),
  curves_ga |> filter(test == GA_TEST, age == GA_AGE) |>
    transmute(source = "Continuous (GAMLSS)", lower = `p2.5`, upper = `p97.5`),
  guven_ri |> filter(test == GA_TEST, age == GA_AGE) |>
    transmute(source = "Guven et al.", lower, upper)
) |>
  mutate(flagged = GA_VALUE < lower | GA_VALUE > upper)

stopifnot(nrow(src_limits) == 5, !any(is.na(src_limits$lower)))

# --- Flagging rates -----------------------------------------------------------
flag_ga <- read_csv("data/processed/analysis6_fp_bootstrap_ci_RAW.csv",
                    show_col_types = FALSE) |>
  filter(test == GA_TEST) |>
  mutate(source = recode(ri_source,
                         refineR = "Indirect (refineR)",
                         Direct  = "Direct (a posteriori)",
                         Continuous = "Continuous (GAMLSS)",
                         Guven   = "Guven et al.",
                         Manufacturer = "Manufacturer"))

SRC_ORDER <- c("Manufacturer", "Guven et al.", "Indirect (refineR)",
               "Continuous (GAMLSS)", "Direct (a posteriori)")
SRC_COL <- c("Manufacturer"          = PAL$base2,
             "Indirect (refineR)"    = PAL$refiner,
             "Direct (a posteriori)" = PAL$direct,
             "Continuous (GAMLSS)"   = PAL$contin,
             "Guven et al."          = PAL$accent3)

src_limits <- src_limits |>
  mutate(source = factor(source, levels = rev(SRC_ORDER)))
flag_ga <- flag_ga |>
  mutate(source = factor(source, levels = rev(SRC_ORDER))) |>
  filter(!is.na(source))

n_aptt <- sum(all_data$test == GA_TEST)

base_theme <- theme_minimal(base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(colour = "grey15", size = 7.2),
    axis.text.x  = element_text(colour = "grey45", size = 6.5),
    axis.title.x = element_text(colour = "grey45", size = 6.8,
                                margin = margin(t = 1)),
    axis.ticks   = element_blank(),
    plot.title   = element_text(size = 7.6, colour = "grey30",
                                margin = margin(b = 4)),
    plot.margin  = margin(2, 4, 2, 2)
  )

# --- Left panel: one result against every source ------------------------------
xlo <- min(src_limits$lower); xhi <- max(src_limits$upper)
span <- xhi - xlo
lab_x <- xhi + span * 0.07

p_left <- ggplot(src_limits, aes(y = source)) +
  geom_segment(aes(x = lower, xend = upper, yend = source, colour = source),
               linewidth = 3.2, lineend = "round") +
  geom_vline(xintercept = GA_VALUE, colour = "grey25",
             linewidth = 0.4, linetype = "22") +
  geom_point(aes(x = GA_VALUE), shape = 21, size = 1.8, stroke = 0.5,
             fill = "white", colour = "grey20") +
  geom_text(aes(x = lab_x, label = ifelse(flagged, "flagged", "within"),
                fontface = ifelse(flagged, "bold", "plain"),
                colour = ifelse(flagged, "flag", "ok")),
            hjust = 0, size = 2.2, show.legend = FALSE) +
  scale_colour_manual(values = c(SRC_COL, flag = "#C0392B", ok = "grey45"),
                      guide = "none") +
  scale_x_continuous(limits = c(xlo - span * 0.04, xhi + span * 0.34),
                     expand = c(0, 0)) +
  labs(title = sprintf("An aPTT of %g s at age %d", GA_VALUE, GA_AGE),
       x = "aPTT, seconds") +
  base_theme

# --- Right panel: flagging rates ----------------------------------------------
p_right <- ggplot(flag_ga, aes(y = source, x = fp_mean, fill = source)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(xmin = fp_ci_lo, xmax = fp_ci_hi), orientation = "y",
                width = 0.18, colour = "grey25", linewidth = 0.3) +
  geom_vline(xintercept = 5, colour = "grey25", linewidth = 0.35, linetype = "22") +
  geom_text(aes(label = sprintf("%.1f%%", fp_mean)),
            hjust = -0.25, size = 2.3, colour = "grey15") +
  scale_fill_manual(values = SRC_COL, guide = "none") +
  scale_x_continuous(limits = c(0, max(flag_ga$fp_ci_hi) * 1.3),
                     expand = c(0, 0),
                     breaks = c(0, 5, 10, 20, 30),
                     labels = c("0", "5%\nexpected", "10", "20", "30")) +
  labs(title = "Flagged outside the limits",
       x = sprintf("Percent of all %s aPTT results",
                   format(n_aptt, big.mark = ","))) +
  base_theme +
  theme(axis.text.y = element_blank())

ga <- p_left + p_right +
  plot_layout(widths = c(1, 0.85)) +
  plot_annotation(
    title = "Reference interval source changes how a pediatric result is classified",
    theme = theme(plot.title = element_text(size = 8.8, face = "bold",
                                            colour = "grey10",
                                            margin = margin(b = 5)),
                  plot.margin = margin(5, 5, 3, 5))
  ) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

# --- Export -------------------------------------------------------------------
# 13 x 5 cm; 600 dpi TIFF and 300 dpi PNG, as elsewhere in the project.
W_CM <- 13; H_CM <- 5
for (d in c("figures/TIFF_600DPI", "figures/PNG_300DPI"))
  if (!dir.exists(d)) stop("missing output directory: ", d)

ggsave("figures/TIFF_600DPI/graphical_abstract.tiff", ga,
       width = W_CM, height = H_CM, units = "cm", dpi = 600,
       bg = "white", compression = "lzw")
ggsave("figures/PNG_300DPI/graphical_abstract.png", ga,
       width = W_CM, height = H_CM, units = "cm", dpi = 300, bg = "white")

for (f in c("figures/TIFF_600DPI/graphical_abstract.tiff",
            "figures/PNG_300DPI/graphical_abstract.png")) {
  message(sprintf("[OK] %s (%.0f KB)", f, file.size(f) / 1024))
}
message(sprintf("     %d x %d px at 600 dpi; Elsevier asks for at least 1328 x 531",
                round(W_CM / 2.54 * 600), round(H_CM / 2.54 * 600)))
