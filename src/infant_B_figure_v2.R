# ==============================================================================
# INFANT SUPPLEMENTARY FIGURE -- vertical layout, infant in blue
# Updates vs previous:
#  - Vertical bars (categorical x = age group, numeric y = value)
#  - Two-colour fill: infant (30-365 d) blue (PAL$direct), pediatric orange (PAL$refiner)
#  - Errorbars vertical (90% CI)
#  - Tufte theme + bold strip titles + bottom legend (Figure 3 conventions)
# Panel layout:
#   A. PT          - 30-365 d  |  1-18 y
#   B. aPTT        - 30-365 d  |  1-12 y  |  12-18 y
#   C. Fibrinogen  - 30-365 d  |  1-18 y
# Overwrites figures/{TIFF_600DPI,PNG_300DPI}/figure_supplement_infant_RI.{tiff,png}
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(patchwork); library(forcats)
})

if (!exists("PAL")) source("src/common.R")

# --- Pediatric refineR RIs (with 90% CI) -----------------------------------
peds <- read_csv("data/processed/tables/PUB_table_2_reference_intervals.csv",
                 show_col_types = FALSE) |>
  filter(method == "refineR") |>
  mutate(grp = "Pediatric") |>
  select(test, age_group, grp, lower, lower_ci_lo, lower_ci_hi,
         upper, upper_ci_lo, upper_ci_hi)

# --- Infant refineR RIs (with 90% CI) --------------------------------------
infant <- read_csv("data/processed/infant_B_refineR_RI.csv",
                   show_col_types = FALSE) |>
  mutate(age_group = "30-365 d",
         grp = "Infant") |>
  select(test, age_group, grp, lower, lower_ci_lo, lower_ci_hi,
         upper, upper_ci_lo, upper_ci_hi)

# --- Per-panel data --------------------------------------------------------
pt_data <- bind_rows(
  infant |> filter(test == "PT"),
  peds   |> filter(test == "PT")
) |>
  mutate(age_group = factor(age_group, levels = c("30-365 d", "1-18 y")))

aptt_data <- bind_rows(
  infant |> filter(test == "aPTT"),
  peds   |> filter(test == "aPTT")
) |>
  mutate(age_group = factor(age_group, levels = c("30-365 d", "1-12 y", "12-18 y")))

fib_data <- bind_rows(
  infant |> filter(test == "Fibrinogen"),
  peds   |> filter(test == "Fibrinogen")
) |>
  mutate(age_group = factor(age_group, levels = c("30-365 d", "1-18 y")))

# --- Panel builder ---------------------------------------------------------
make_panel <- function(df, panel_title, y_label) {
  ggplot(df, aes(x = age_group, fill = grp)) +
    geom_crossbar(aes(y = (lower + upper)/2, ymin = lower, ymax = upper),
                  width = 0.45, linewidth = 0.4, alpha = 0.78, colour = NA) +
    geom_errorbar(aes(ymin = lower_ci_lo, ymax = lower_ci_hi),
                  width = 0.18, linewidth = 0.4, colour = "#444444") +
    geom_errorbar(aes(ymin = upper_ci_lo, ymax = upper_ci_hi),
                  width = 0.18, linewidth = 0.4, colour = "#444444") +
    scale_fill_manual(values = c("Infant"    = PAL$direct,
                                 "Pediatric" = PAL$refiner),
                      labels = c("Infant"    = "Infant (30-365 d)",
                                 "Pediatric" = "Pediatric (1-18 y)"),
                      breaks = c("Infant", "Pediatric")) +
    labs(title = panel_title, x = NULL, y = y_label, fill = NULL) +
    theme_tufte_academic() +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0),
          axis.text.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          panel.spacing = unit(10, "pt"),
          plot.margin = margin(4, 6, 4, 4),
          legend.position = "none")
}

p_pt   <- make_panel(pt_data,   "PT",         "PT, seconds")
p_aptt <- make_panel(aptt_data, "aPTT",       "aPTT, seconds")
p_fib  <- make_panel(fib_data,  "Fibrinogen", "Fibrinogen, g/L")

fig <- (p_pt | p_aptt | p_fib) +
  plot_layout(widths = c(1, 1.3, 1)) +
  plot_annotation(tag_levels = "A",
                  theme = theme_tufte_academic()) &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_fig(fig, "figure_supplement_infant_RI", width = 220, height = 90)
cat("[OK] Figure saved (vertical, no legend; overwrites previous).\n")
