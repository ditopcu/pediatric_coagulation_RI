# ==============================================================================
# COMMON: Shared libraries, design system, data loading
# Pediatric Coagulation Reference Intervals
# ==============================================================================
# Source this file from pub.R and analysis.R to avoid duplication.
# Provides: PAL, theme_tufte_academic, add_panel_tag, save_fig,
#           all_data, all_data_3sd, filter_3sd, ri_comparison,
#           ri_manufacturer, cont_ri_all, guven_ri, ri_refiner, ri_direct

# --- Libraries ---
library(tidyverse)
library(readxl)
library(gamlss)
library(patchwork)

# ==============================================================================
# DESIGN SYSTEM
# ==============================================================================

PAL <- list(
  highlight = "#C0392B",
  base1     = "#5D6D7E",
  base2     = "#BDC3C7",
  accent1   = "#27AE60",
  accent2   = "#E67E22",
  accent3   = "#8E44AD",
  female    = "#C0392B",
  male      = "#2980B9",
  refiner   = "#E67E22",
  direct    = "#2980B9",
  contin    = "#27AE60",
  ci_band   = "#BDC3C7",
  text      = "#333333",
  grid      = "#E8E8E8"
)

theme_tufte_academic <- function(base_size = 11, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line          = element_line(colour = "#555555", linewidth = 0.5),
      axis.line.x        = element_line(colour = "#555555", linewidth = 0.5),
      axis.line.y        = element_line(colour = "#555555", linewidth = 0.5),
      axis.ticks         = element_line(colour = "#555555", linewidth = 0.3),
      axis.ticks.length  = unit(3, "pt"),
      axis.text          = element_text(size = rel(0.9), colour = PAL$text),
      axis.title         = element_text(size = rel(1.0), colour = PAL$text,
                                        margin = margin(t = 4, r = 4)),
      axis.title.y       = element_text(angle = 90, margin = margin(r = 8)),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.background   = element_rect(fill = "white", colour = NA),
      panel.border       = element_blank(),
      legend.background  = element_blank(),
      legend.key         = element_blank(),
      legend.title       = element_blank(),
      legend.text        = element_text(size = rel(0.85)),
      legend.position    = "bottom",
      strip.background   = element_blank(),
      strip.text         = element_text(size = rel(0.95), face = "bold",
                                        colour = PAL$text, hjust = 0),
      plot.title         = element_text(size = rel(1.1), face = "bold",
                                        colour = PAL$text, hjust = 0,
                                        margin = margin(b = 8)),
      plot.subtitle      = element_text(size = rel(0.85), colour = "#666666",
                                        hjust = 0, margin = margin(b = 10)),
      plot.margin        = margin(12, 12, 12, 12),
      plot.background    = element_rect(fill = "white", colour = NA)
    )
}

add_panel_tag <- function(label) {
  labs(tag = label) +
    theme(plot.tag = element_text(size = 14, face = "bold", hjust = 0))
}

save_fig <- function(plot, filename, width = 180, height = 120,
                     units = "mm", dpi_tiff = 600, dpi_png = 300) {
  dir.create("figures/TIFF_600DPI", recursive = TRUE, showWarnings = FALSE)
  dir.create("figures/PNG_300DPI",  recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path("figures/TIFF_600DPI", paste0(filename, ".tiff")),
         plot = plot, width = width, height = height, units = units,
         dpi = dpi_tiff, bg = "white", compression = "lzw")
  ggsave(file.path("figures/PNG_300DPI", paste0(filename, ".png")),
         plot = plot, width = width, height = height, units = units,
         dpi = dpi_png, bg = "white")
  message("[OK] Saved: ", filename)
}

# ==============================================================================
# DATA LOADING
# ==============================================================================

message("[INFO] Reading data...")

ptz_data <- read_excel("data/coa_results/tce 2026 1-18 pt.xlsx") |>
  janitor::clean_names() |>
  select(sample_id, sex, age_year, result_num) |>
  mutate(age_int = floor(age_year), test = "PT")

aptt_data <- read_excel("data/coa_results/tce 2026 1-18 aptt.xlsx") |>
  janitor::clean_names() |>
  select(sample_id, sex, age_year, result_num) |>
  mutate(age_int = floor(age_year), test = "aPTT")

fib_data <- read_excel("data/coa_results/tce 2026 1-18 fib.xlsx") |>
  janitor::clean_names() |>
  select(sample_id, sex, age_year, result_num) |>
  mutate(age_int = floor(age_year), test = "Fibrinogen")

all_data <- bind_rows(ptz_data, aptt_data, fib_data)

filter_3sd <- function(df) {
  df |>
    group_by(test) |>
    mutate(m = mean(result_num, na.rm = TRUE),
           s = sd(result_num, na.rm = TRUE)) |>
    filter(result_num >= m - 3*s, result_num <= m + 3*s) |>
    select(-m, -s) |>
    ungroup()
}

all_data_3sd <- filter_3sd(all_data)

message(sprintf("[OK] Data loaded: %d observations (%d after 3SD filter).",
                nrow(all_data), nrow(all_data_3sd)))

# ==============================================================================
# REFERENCE INTERVAL DEFINITIONS
# ==============================================================================

ri_comparison <- tribble(
  ~test,         ~age_group,  ~method,    ~n,     ~lower, ~lower_ci_lo, ~lower_ci_hi, ~upper, ~upper_ci_lo, ~upper_ci_hi,
  "PT",          "1-18 y",    "refineR",  17000,  8.56,   8.54,         8.58,         10.90,  10.70,        10.90,
  "PT",          "1-18 y",    "Direct",   305,    8.69,   8.55,         8.84,         10.80,  10.60,        11.20,
  "aPTT",        "1-12 y",   "refineR",  10091,  24.30,  23.30,        24.50,        35.30,  34.40,        35.40,
  "aPTT",        "1-12 y",   "Direct",   170,    22.96,  22.40,        23.40,        34.74,  34.20,        35.60,
  "aPTT",        "12-18 y",  "refineR",  6910,   23.30,  22.70,        23.70,        33.90,  33.10,        34.10,
  "aPTT",        "12-18 y",  "Direct",   150,    24.10,  23.50,        24.30,        32.69,  31.60,        33.40,
  "Fibrinogen",  "1-18 y",   "refineR",  3166,   1.89,   1.84,         1.92,         3.79,   3.55,         4.06,
  "Fibrinogen",  "1-18 y",   "Direct",   198,    2.10,   2.07,         2.16,         3.56,   3.37,         3.89
)

ri_manufacturer <- tribble(
  ~test,         ~lower, ~upper,
  "PT",          8.4,    10.6,
  "aPTT",        23.6,   30.6,
  "Fibrinogen",  1.93,   4.12
)

ri_refiner <- tribble(
  ~test,         ~age_group,  ~lower, ~upper,
  "PT",          "1-18",       8.56,  10.90,
  "aPTT",        "1-12",      24.30,  35.30,
  "aPTT",        "12-18",     23.30,  33.90,
  "Fibrinogen",  "1-18",       1.89,   3.79
)

ri_direct <- tribble(
  ~test,         ~age_group,  ~lower, ~upper,
  "PT",          "1-18",       8.69,  10.80,
  "aPTT",        "1-12",      22.96,  34.74,
  "aPTT",        "12-18",     24.10,  32.69,
  "Fibrinogen",  "1-18",       2.10,   3.56
)

guven_ri <- tribble(
  ~test,  ~age, ~n,   ~lower, ~upper,
  "PT",    1,    68,   7.8,   10.3,
  "PT",    2,   116,   7.8,   10.1,
  "PT",    3,    57,   7.8,   10.0,
  "PT",    4,    26,   7.9,   10.0,
  "PT",    5,    57,   8.0,   10.0,
  "PT",    6,    97,   8.0,   10.1,
  "PT",    7,   137,   8.0,   10.2,
  "PT",    8,   140,   8.1,   10.3,
  "PT",    9,    97,   7.9,   10.1,
  "PT",   10,    51,   7.9,   10.1,
  "PT",   11,    26,   7.9,   10.0,
  "PT",   12,    38,   8.0,   10.0,
  "PT",   13,    27,   8.0,   10.0,
  "PT",   14,    22,   8.1,   10.1,
  "PT",   15,    19,   8.1,   10.2,
  "PT",   16,    15,   8.1,   10.3,
  "PT",   17,    25,   8.2,   10.5,
  "PT",   18,    14,   8.2,   10.7,
  "aPTT",  1,    66,  25.2,   35.6,
  "aPTT",  2,   115,  24.6,   34.9,
  "aPTT",  3,    59,  24.3,   34.5,
  "aPTT",  4,    26,  24.3,   34.4,
  "aPTT",  5,    53,  24.5,   34.6,
  "aPTT",  6,    96,  24.8,   34.8,
  "aPTT",  7,   136,  25.1,   34.9,
  "aPTT",  8,   133,  25.2,   35.1,
  "aPTT",  9,    97,  25.2,   36.2,
  "aPTT", 10,    50,  24.9,   35.2,
  "aPTT", 11,    26,  24.6,   34.4,
  "aPTT", 12,    40,  24.4,   33.9,
  "aPTT", 13,    29,  24.3,   33.6,
  "aPTT", 14,    23,  24.4,   33.5,
  "aPTT", 15,    19,  24.6,   33.7,
  "aPTT", 16,    18,  24.9,   34.3,
  "aPTT", 17,    24,  25.4,   35.1,
  "aPTT", 18,    12,  24.7,   32.4
)

# ==============================================================================
# CONTINUOUS RI FROM GAMLSS PIPELINE
# ==============================================================================

source("src/pipeline/utils.R")

get_continuous_ri <- function(point_est_file, bs_files) {
  load(point_est_file)
  pp <- pointEstGamlss$Params$transformP
  assign("pp", pp, envir = .GlobalEnv)
  assign("ptrans", function(x, p) if (p == 0) log(x) else I(x^p), envir = .GlobalEnv)
  assign("families", pointEstGamlss$Models$family, envir = .GlobalEnv)
  assign("fam", 1, envir = .GlobalEnv)
  assign("covarName", "Age", envir = .GlobalEnv)

  bs_results <- lapply(bs_files, function(f) {
    if (file.exists(f)) { load(f); return(gamlssModel) } else return(NULL)
  })
  bs_results <- Filter(Negate(is.null), bs_results)

  # IMPORTANT: must use all 7 default percentiles so CILow/CIHigh indices
  # match RIMat rows (both are indexed 1..7 inside estimateCIs).
  # Using fewer percentiles causes index mismatch for CIs.
  full_perc <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
  allRes <- estimateCIs(pointEst = pointEstGamlss, estBS = bs_results,
                        covarValue = 1:18, withData = TRUE,
                        RIperc = full_perc, CIprop = 0.90)

  ri_mat <- allRes$RICurve$RIMat
  # Indices for our target percentiles within the 7-element CILow/CIHigh lists
  idx_025 <- 1  # 0.025 is 1st
  idx_50  <- 4  # 0.5   is 4th
  idx_975 <- 7  # 0.975 is 7th

  list(
    ri_df = tibble(age = 1:18,
                   lower = ri_mat["Perc_0.025", ],
                   p50   = ri_mat["Perc_0.5", ],
                   upper = ri_mat["Perc_0.975", ]),
    allRes = allRes,
    ci_idx = list(lower = idx_025, median = idx_50, upper = idx_975)
  )
}

message("[INFO] Loading continuous RI curves...")

cont_pt   <- get_continuous_ri("cont_out/PT_gamlss_PointEst.RData",
                               paste0("cont_out/PT_gamlssModel_Est_", 1:5, ".RData"))
cont_aptt <- get_continuous_ri("cont_out/aPTT_gamlss_PointEst.RData",
                               paste0("cont_out/aPTT_gamlssModel_Est_", 1:5, ".RData"))
cont_fib  <- get_continuous_ri("cont_out/Fib_gamlss_PointEst.RData",
                               paste0("cont_out/Fib_gamlssModel_Est_", 1:5, ".RData"))

cont_ri_all <- bind_rows(
  cont_pt$ri_df   |> mutate(test = "PT"),
  cont_aptt$ri_df |> mutate(test = "aPTT"),
  cont_fib$ri_df  |> mutate(test = "Fibrinogen")
)

message("[OK] Common setup complete.\n")
