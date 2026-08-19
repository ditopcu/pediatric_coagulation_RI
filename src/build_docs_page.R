# ==============================================================================
# BUILD docs/index.html -- interactive reference interval page
# ==============================================================================
# Assembles a single self-contained HTML file from the aggregated analysis
# outputs. The page carries no patient-level data: only reference limits,
# per-age percentile curves, flagging rates and counts.
#
# The numbers are read from the same files the manuscript is traced to in
# docs/PROVENANCE.md, so the page cannot carry a second copy that drifts.
#
# Usage: setwd(project root) then source("src/build_docs_page.R")
# ==============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(jsonlite)
})

# common.R holds ri_manufacturer, ri_refiner, ri_direct and guven_ri.
if (!exists("guven_ri")) source("src/common.R")

message("[BUILD] Assembling docs/index.html ...")

TBL <- "data/processed/tables"
PROC <- "data/processed"

AGES <- 1:17   # the study covers 1 to <18 years

# --- Sources of every panel, shown in the page footer -------------------------
provenance <- tribble(
  ~panel,                    ~file,
  "Result classifier",       "PUB_table_2_reference_intervals.csv, PUB_table_S3_continuous_RI_values.csv, guven_ri in src/common.R",
  "Reference intervals",     "PUB_table_2_reference_intervals.csv",
  "Continuous curves",       "PUB_table_S3_continuous_RI_values.csv",
  "Out-of-RI flagging",      "analysis3_fp_curves_data_RAW.csv, analysis6_fp_bootstrap_ci_RAW.csv, analysis1_fp_summary_RAW.csv",
  "Partitioning",            "PUB_table_S2_refineR_all_by_age_sex.csv, PUB_table_6_harris_boyd.csv",
  "Transferability",         "analysis2_guven_transferability_RAW.csv",
  "Sample size by age",      "analysis1_fp_by_age_RAW.csv",
  "Published comparison",    "PUB_table_5_international_comparison.csv"
)

# --- Reference intervals ------------------------------------------------------
intervals <- read_csv(file.path(TBL, "PUB_table_2_reference_intervals.csv"),
                      show_col_types = FALSE) |>
  select(test, age_group, method, source, n,
         lower, lower_ci_lo, lower_ci_hi, upper, upper_ci_lo, upper_ci_hi)

# --- Continuous per-age curves ------------------------------------------------
curves <- read_csv(file.path(TBL, "PUB_table_S3_continuous_RI_values.csv"),
                   show_col_types = FALSE) |>
  filter(age %in% AGES) |>
  rename(p2_5 = `p2.5`, p2_5_lo = `p2.5_ci_lo`, p2_5_hi = `p2.5_ci_hi`,
         p50_lo = p50_ci_lo, p50_hi = p50_ci_hi,
         p97_5 = `p97.5`, p97_5_lo = `p97.5_ci_lo`, p97_5_hi = `p97.5_ci_hi`)

# --- Flagging rates -----------------------------------------------------------
# Fibrinogen carries Guven et al. rows with no rate, since that study did not
# report fibrinogen. They are dropped rather than plotted as gaps.
flag_age <- read_csv(file.path(PROC, "analysis3_fp_curves_data_RAW.csv"),
                     show_col_types = FALSE) |>
  filter(age_int %in% AGES, !is.na(fp_rate)) |>
  select(test, age = age_int, n, ri_source, fp_rate)

flag_overall <- read_csv(file.path(PROC, "analysis6_fp_bootstrap_ci_RAW.csv"),
                         show_col_types = FALSE) |>
  select(test, ri_source, fp_mean, fp_ci_lo, fp_ci_hi, fp_label)

# Point-estimate rates, from which the absolute number of flagged results is
# recovered exactly. The bootstrap mean in flag_overall is not used for this.
flag_counts <- read_csv(file.path(PROC, "analysis1_fp_summary_RAW.csv"),
                        show_col_types = FALSE) |>
  pivot_longer(starts_with("FP_"), names_to = "ri_source", values_to = "rate",
               names_prefix = "FP_") |>
  filter(!is.na(rate)) |>
  mutate(ri_source = recode(ri_source, refineR = "refineR", Continuous = "Continuous",
                            Direct = "Direct", Guven = "Guven", Manufacturer = "Manufacturer"),
         count = round(rate / 100 * N)) |>
  select(test, ri_source, n_total = N, rate, count)

# --- Partitioning: sex-specific limits by age, and the Harris-Boyd statistic ---
sex_limits <- read_csv(file.path(TBL, "PUB_table_S2_refineR_all_by_age_sex.csv"),
                       show_col_types = FALSE) |>
  filter(age_int %in% AGES) |>
  transmute(test, age = age_int,
            sex = recode(sex, K = "Female", E = "Male", Combined = "Combined"),
            n, lower = `lower_2.5`, lower_lo = lower_ci_lo, lower_hi = lower_ci_hi,
            upper = `upper_97.5`, upper_lo = upper_ci_lo, upper_hi = upper_ci_hi)

harris_boyd <- read_csv(file.path(TBL, "PUB_table_6_harris_boyd.csv"),
                        show_col_types = FALSE) |>
  select(test, cutpoint, n_young, n_old, z_score, justified)

# The sex-combined reference limits and their 90% confidence intervals, drawn
# behind the sex-specific points exactly as in manuscript Figure 2. One band per
# reported partition: a single 1-18 y band for PT and fibrinogen, and two for
# aPTT, split at 12 years. Age spans follow the partition labels in common.R,
# where "1-12 y" is age_year < 12; the data end at age 17.
span_of <- function(group) {
  switch(group, "1-18 y" = c(1L, 17L), "1-12 y" = c(1L, 11L), "12-18 y" = c(12L, 17L),
         stop("unknown age_group: ", group))
}

bands_src <- ri_comparison |> filter(method == "refineR")
partition_bands <- bind_rows(lapply(seq_len(nrow(bands_src)), function(i) {
  r <- bands_src[i, ]; s <- span_of(r$age_group)
  tibble(test = r$test, age_group = r$age_group,
         age_start = s[1], age_end = s[2],
         lower = r$lower, lower_lo = r$lower_ci_lo, lower_hi = r$lower_ci_hi,
         upper = r$upper, upper_lo = r$upper_ci_lo, upper_hi = r$upper_ci_hi)
}))

stopifnot(nrow(partition_bands) == nrow(bands_src))

# --- Transferability of the published Turkish intervals -----------------------
transfer <- read_csv(file.path(PROC, "analysis2_guven_transferability_RAW.csv"),
                     show_col_types = FALSE) |>
  filter(age %in% AGES, !is.na(pass_rate)) |>
  select(test, age, n_total, pass_rate, prop_outside, transferable)

# --- Sample size by single year of age ----------------------------------------
# Taken from the flagging output rather than the descriptive table, because its
# per-age counts sum to the canonical totals the page reports in its header.
sample_size <- read_csv(file.path(PROC, "analysis1_fp_by_age_RAW.csv"),
                        show_col_types = FALSE) |>
  filter(age_int %in% AGES) |>
  select(test, age = age_int, n)

local({
  got <- sample_size |> group_by(test) |> summarise(n = sum(n), .groups = "drop")
  want <- c(PT = 17000L, aPTT = 17001L, Fibrinogen = 3166L)
  for (t in names(want)) {
    have <- got$n[got$test == t]
    if (length(have) != 1 || have != want[[t]])
      stop(sprintf("sample_size for %s sums to %s, expected %d",
                   t, paste(have, collapse = "/"), want[[t]]))
  }
})

# --- Published intervals from other studies -----------------------------------
# This study's own rows are dropped: they already have their own panel. The
# Weidhofer PT row carries no comparable limits and is excluded by the NA filter.
published <- read_csv(file.path(TBL, "PUB_table_5_international_comparison.csv"),
                      show_col_types = FALSE) |>
  filter(!grepl("^This study", study), !is.na(lower), !is.na(upper)) |>
  select(study, year, country, analyzer, method, n, age_group, test,
         lower, upper, unit)

# --- Per-age limits for every source, used by the classifier ------------------
# The partition labels follow common.R: "1-12" is age_year < 12.
expand_by_age <- function(ri, method) {
  ri |>
    rowwise() |>
    mutate(age = list(switch(age_group,
                             "1-18"  = AGES,
                             "1-12"  = AGES[AGES <  12],
                             "12-18" = AGES[AGES >= 12]))) |>
    ungroup() |>
    unnest(age) |>
    transmute(test, age, method = method, lower, upper)
}

limits <- bind_rows(
  ri_manufacturer |> crossing(age = AGES) |>
    transmute(test, age, method = "Manufacturer", lower, upper),
  expand_by_age(ri_refiner, "Indirect (refineR)"),
  expand_by_age(ri_direct,  "Direct (a posteriori)"),
  curves |> transmute(test, age, method = "Continuous (GAMLSS)",
                      lower = p2_5, upper = p97_5),
  guven_ri |> filter(age %in% AGES) |>
    transmute(test, age, method = "Guven et al.", lower, upper)
) |>
  arrange(test, age, method)

stopifnot(!any(is.na(limits$lower)), !any(is.na(limits$upper)))

# --- Study population ---------------------------------------------------------
# Counts come from the canonical files as loaded by common.R.
population <- all_data |>
  count(test, name = "n") |>
  mutate(test = factor(test, levels = c("PT", "aPTT", "Fibrinogen"))) |>
  arrange(test) |>
  mutate(test = as.character(test),
         unit = recode(test, PT = "seconds", aPTT = "seconds",
                       Fibrinogen = "g/L"))

payload <- list(
  built       = format(Sys.Date()),
  ages        = AGES,
  population  = population,
  intervals   = intervals,
  curves      = curves,
  flagAge     = flag_age,
  flagOverall = flag_overall,
  flagCounts  = flag_counts,
  limits      = limits,
  sexLimits   = sex_limits,
  bands       = partition_bands,
  harrisBoyd  = harris_boyd,
  transfer    = transfer,
  sampleSize  = sample_size,
  published   = published,
  provenance  = provenance
)

json <- toJSON(payload, dataframe = "rows", auto_unbox = TRUE,
               na = "null", digits = 6)

message(sprintf("  intervals %d | curves %d | flagAge %d | limits %d",
                nrow(intervals), nrow(curves), nrow(flag_age), nrow(limits)))
message(sprintf("  sexLimits %d | transfer %d | sampleSize %d | published %d",
                nrow(sex_limits), nrow(transfer), nrow(sample_size), nrow(published)))

# ==============================================================================
# HTML TEMPLATE
# ==============================================================================

tmpl <- r"---(<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pediatric Coagulation Reference Intervals</title>
<meta name="description" content="Interactive comparison of reference interval sources for PT, aPTT and fibrinogen in children aged 1 to under 18 years.">
<style>
:root {
  --bg: #ffffff; --panel: #fafaf9; --ink: #1a1a1a; --muted: #6b6b6b;
  --line: #e0ddd8; --rule: #cfcbc4;
  --refiner: #E67E22; --direct: #2980B9; --contin: #27AE60;
  --mfr: #9aa0a6; --guven: #8E44AD;
  --female: #C0392B; --male: #2980B9; --combined: #6b6b6b;
  --in: #27AE60; --out: #C0392B;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #16171a; --panel: #1e2024; --ink: #e8e6e3; --muted: #9b9894;
    --line: #2e3136; --rule: #3a3e44; --mfr: #7d8288;
  }
}
:root[data-theme="dark"] {
  --bg: #16171a; --panel: #1e2024; --ink: #e8e6e3; --muted: #9b9894;
  --line: #2e3136; --rule: #3a3e44; --mfr: #7d8288;
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--bg); color: var(--ink);
  font: 16px/1.6 Georgia, "Iowan Old Style", "Times New Roman", serif;
  -webkit-text-size-adjust: 100%;
}
.wrap { max-width: 62rem; margin: 0 auto; padding: 2.5rem 1.25rem 4rem; }
h1 { font-size: 1.9rem; line-height: 1.25; margin: 0 0 .4rem; font-weight: 600; letter-spacing: -.01em; }
h2 { font-size: 1.3rem; margin: 3.2rem 0 .3rem; font-weight: 600; letter-spacing: -.005em; }
h3 { font-size: .95rem; margin: 1.6rem 0 .5rem; font-weight: 600; }
p { margin: .5rem 0; }
.lede { color: var(--muted); font-size: 1.02rem; max-width: 46rem; }
.note {
  font-size: .82rem; color: var(--muted); border-left: 2px solid var(--rule);
  padding: .1rem 0 .1rem .8rem; margin: 1.2rem 0; max-width: 46rem;
}
.sub { font-size: .85rem; color: var(--muted); margin: .2rem 0 1rem; max-width: 46rem; }
hr { border: 0; border-top: 1px solid var(--line); margin: 0; }
nav { margin: 1.6rem 0 0; font-size: .85rem; }
nav a { color: var(--muted); text-decoration: none; margin-right: 1.1rem; border-bottom: 1px solid var(--line); }
nav a:hover { color: var(--ink); border-color: var(--ink); }
.stats { display: flex; flex-wrap: wrap; gap: 2rem; margin: 1.4rem 0 0; }
.stat b { display: block; font-size: 1.35rem; font-weight: 600; font-variant-numeric: tabular-nums; }
.stat span { font-size: .78rem; color: var(--muted); text-transform: uppercase; letter-spacing: .06em; }

.panel { background: var(--panel); border: 1px solid var(--line); border-radius: 6px; padding: 1.1rem 1.2rem; margin: 1rem 0; }
.controls { display: flex; flex-wrap: wrap; gap: 1rem 1.4rem; align-items: flex-end; }
.field { display: flex; flex-direction: column; gap: .25rem; }
.field label { font-size: .74rem; text-transform: uppercase; letter-spacing: .06em; color: var(--muted); }
select, input[type=number] {
  font: inherit; font-size: .95rem; padding: .4rem .55rem; color: var(--ink);
  background: var(--bg); border: 1px solid var(--rule); border-radius: 4px; min-width: 9rem;
}
input[type=range] { width: 15rem; accent-color: var(--refiner); }
.rangewrap { display: flex; align-items: center; gap: .7rem; }
.agebadge { font-variant-numeric: tabular-nums; font-weight: 600; min-width: 4.5rem; }

table { border-collapse: collapse; width: 100%; font-size: .9rem; }
th, td { text-align: left; padding: .5rem .6rem; border-bottom: 1px solid var(--line); }
th { font-size: .72rem; text-transform: uppercase; letter-spacing: .06em; color: var(--muted); font-weight: 600; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
.scroll { overflow-x: auto; -webkit-overflow-scrolling: touch; }
.swatch { display: inline-block; width: .62rem; height: .62rem; border-radius: 2px; margin-right: .5rem; vertical-align: baseline; }
.verdict { font-weight: 600; font-size: .84rem; letter-spacing: .02em; }
.verdict.in { color: var(--in); }
.verdict.out { color: var(--out); }
.tally { margin-top: .9rem; font-size: .92rem; }
.tally strong { font-variant-numeric: tabular-nums; }

figure { margin: 1rem 0 0; }
figcaption { font-size: .8rem; color: var(--muted); margin-top: .5rem; max-width: 46rem; }
svg { display: block; width: 100%; height: auto; overflow: visible; }
svg text { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; fill: var(--muted); }
svg .axis line, svg .axis path { stroke: var(--rule); }
svg .grid line { stroke: var(--line); }
.legend { display: flex; flex-wrap: wrap; gap: .35rem 1.2rem; font-size: .8rem; margin: .8rem 0 0; color: var(--muted); }

/* sticky section nav */
nav.sticky {
  position: sticky; top: 0; z-index: 20; background: var(--bg);
  border-bottom: 1px solid var(--line); padding: .6rem 0; margin: 1.6rem 0 0;
}
nav a.here { color: var(--ink); border-color: var(--ink); }

/* analyte tabs */
.tabs { display: flex; gap: .2rem; margin: 1rem 0 .2rem; border-bottom: 1px solid var(--line); flex-wrap: wrap; }
.tabs button {
  font: inherit; font-size: .88rem; color: var(--muted); background: none;
  border: 0; border-bottom: 2px solid transparent; padding: .45rem .85rem;
  cursor: pointer; margin-bottom: -1px;
}
.tabs button:hover { color: var(--ink); }
.tabs button[aria-selected="true"] { color: var(--ink); border-bottom-color: var(--refiner); font-weight: 600; }
.tabpanel[hidden] { display: none; }

button.action {
  font: inherit; font-size: .82rem; color: var(--ink); background: var(--panel);
  border: 1px solid var(--rule); border-radius: 4px; padding: .42rem .8rem; cursor: pointer;
}
button.action:hover { border-color: var(--ink); }
.actions { display: flex; gap: .6rem; align-items: center; flex-wrap: wrap; margin: .9rem 0 0; }
.actions .hint { font-size: .78rem; color: var(--muted); }

@media print {
  nav, .actions, input[type=range] { display: none !important; }
  nav.sticky { position: static; }
  .tabs { display: none; }
  .tabpanel[hidden] { display: block !important; }
  body { background: #fff; color: #000; font-size: 10pt; }
  .wrap { max-width: none; padding: 0; }
  h2 { break-after: avoid; }
  figure, .panel { break-inside: avoid; }
  .panel { background: none; border: 1px solid #ccc; }
  a[href^="http"]::after { content: " (" attr(href) ")"; font-size: 8pt; }
}

footer { margin-top: 4rem; padding-top: 1.4rem; border-top: 1px solid var(--line); font-size: .84rem; color: var(--muted); }
footer a { color: inherit; }
code { font: .85em ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; background: var(--panel); padding: .08em .32em; border-radius: 3px; }
@media (max-width: 40rem) {
  .wrap { padding: 1.6rem 1rem 3rem; }
  h1 { font-size: 1.5rem; }
  input[type=range] { width: 11rem; }
}
</style>
</head>
<body>
<div class="wrap">

<h1>Pediatric coagulation reference intervals</h1>
<p class="lede">How the choice of reference interval changes the way a routine
PT, aPTT or fibrinogen result is classified in children aged 1 to under 18
years. Roche Cobas t511/t711, a single tertiary centre.</p>

<div class="stats" id="stats"></div>

<nav class="sticky" id="nav">
  <a href="#classify">Classify a result</a>
  <a href="#intervals">Reference intervals</a>
  <a href="#curves">Continuous curves</a>
  <a href="#flagging">Out-of-RI flagging</a>
  <a href="#partition">Partitioning</a>
  <a href="#transfer">Transferability</a>
  <a href="#published">Published intervals</a>
</nav>

<p class="note">This page contains aggregated results only &mdash; reference
limits, percentile curves, flagging rates and counts. No patient-level record
is included and no individual result can be recovered from it.</p>

<h2 id="classify">Classify a result</h2>
<p class="sub">Enter an age and a measured value. Each reference interval source
is applied to it side by side. This is the comparison the study is about: the
same number can be inside one interval and outside another.</p>

<div class="panel">
  <div class="controls">
    <div class="field">
      <label for="c-test">Analyte</label>
      <select id="c-test"></select>
    </div>
    <div class="field">
      <label for="c-age">Age</label>
      <div class="rangewrap">
        <input type="range" id="c-age" min="1" max="17" step="1" value="8">
        <span class="agebadge" id="c-agelab"></span>
      </div>
    </div>
    <div class="field">
      <label for="c-val">Result <span id="c-unit"></span></label>
      <input type="number" id="c-val" step="0.1">
    </div>
  </div>
  <div class="scroll"><table id="c-table">
    <thead><tr><th>Reference interval source</th><th class="num">Lower</th><th class="num">Upper</th><th>Classification</th></tr></thead>
    <tbody></tbody>
  </table></div>
  <p class="tally" id="c-tally"></p>
  <div class="actions">
    <button class="action" id="c-copy">Copy this comparison</button>
    <button class="action" id="c-link">Copy a link to it</button>
    <span class="hint" id="c-hint"></span>
  </div>
</div>

<h2 id="intervals">Reference intervals</h2>
<p class="sub">The intervals derived in this study alongside the
manufacturer-supplied adult interval. Bars span the 2.5th to 97.5th percentile;
whiskers are the 90% confidence interval on each limit where one was estimated.</p>
<div id="intervals-charts"></div>
<div class="actions">
  <button class="action" id="dl-limits">Download all limits as CSV</button>
  <span class="hint">One row per analyte, single year of age and reference
  interval source. Provided for research and reproducibility; any interval
  requires local verification before clinical use.</span>
</div>

<h2 id="curves">Continuous curves</h2>
<p class="sub">Age-resolved 2.5th, 50th and 97.5th percentiles from the GAMLSS
model, with pointwise 90% confidence bands. Move across a panel to read the
values at a single year of age.</p>
<div id="curve-charts"></div>

<h2 id="flagging">Out-of-RI flagging</h2>
<p class="sub">The proportion of routine results falling outside the reference
limits, by single year of age. Under a correctly specified 95% interval about 5%
of results are expected to fall outside, marked by the dotted line.</p>
<div id="flag-charts"></div>

<h2 id="partition">Partitioning by sex and age</h2>
<p class="sub">Why the intervals are reported the way they are. Male and female
reference limits are shown at each single year of age with their 90% confidence
intervals; overlapping intervals across sexes is what led to a combined-sex
interval. Age partitioning was tested separately with the Harris&ndash;Boyd
statistic.</p>
<div class="panel" id="hb-box"></div>
<div id="partition-charts"></div>

<h2 id="transfer">Transferability of the published Turkish intervals</h2>
<p class="sub">The intervals of Guven et al., derived in the same country on the
same analyzer platform, applied to this population under a bootstrap analogue of
the CLSI EP28-A3c 20-sample rule. A group passes when at most 2 of 20 sampled
results fall outside the published interval; the criterion for transferability is
a pass rate of 90% or more, marked by the dashed line.</p>
<div id="transfer-charts"></div>

<h2 id="published">Published intervals from other sources</h2>
<p class="sub">Reference intervals reported by other pediatric studies, together
with the manufacturer's package insert. Read them with care: the analyzer,
reagent and source population differ in every case, and those differences may
contribute to intervals not transferring between settings rather than being
evidence that any one of them is wrong.</p>
<div class="scroll"><table id="pub-table">
  <thead><tr><th>Study</th><th>Analyzer</th><th>Method</th><th>Age group</th>
  <th class="num">n</th><th class="num">Lower</th><th class="num">Upper</th><th>Unit</th></tr></thead>
  <tbody></tbody>
</table></div>
<p class="note">Weidhofer et al. reported PT in percent activity rather than
seconds, so no PT limits are listed for that study, and reported fibrinogen in
mg/dL. Values were read from published figures where the study reported curves
rather than tabulated limits.</p>

<footer>
  <p><strong>Where the numbers come from.</strong></p>
  <div class="scroll"><table id="prov">
    <thead><tr><th>Panel</th><th>Source file</th></tr></thead><tbody></tbody>
  </table></div>
  <p style="margin-top:1.2rem">Analysis code and intermediate outputs:
  <a href="https://github.com/ditopcu/pediatric_coagulation_RI">github.com/ditopcu/pediatric_coagulation_RI</a>.
  Individual-level data are not published. Reference limits use the 2.5th and
  97.5th percentiles throughout. Built <span id="built"></span>.</p>
</footer>

</div>
<script>
const DATA = __PAYLOAD__;

const COLORS = {
  "Indirect (refineR)":    getVar("--refiner"),
  "Direct (a posteriori)": getVar("--direct"),
  "Continuous (GAMLSS)":   getVar("--contin"),
  "Manufacturer":          getVar("--mfr"),
  "Guven et al.":          getVar("--guven"),
  "Female":                getVar("--female"),
  "Male":                  getVar("--male"),
  "Combined":              getVar("--combined")
};
// The stored analysis files label the sources more tersely than the page does.
const ALIAS = {
  "refineR": "Indirect (refineR)",
  "Direct (a posteriori)": "Direct (a posteriori)",
  "Direct": "Direct (a posteriori)",
  "Continuous (GAMLSS)": "Continuous (GAMLSS)",
  "Continuous": "Continuous (GAMLSS)",
  "Manufacturer": "Manufacturer",
  "Guven et al.": "Guven et al.",
  "Guven": "Guven et al."
};
const ORDER = ["Manufacturer", "Indirect (refineR)", "Direct (a posteriori)",
               "Continuous (GAMLSS)", "Guven et al."];
const TESTS = ["PT", "aPTT", "Fibrinogen"];
const SVGNS = "http://www.w3.org/2000/svg";

function getVar(n) {
  return getComputedStyle(document.documentElement).getPropertyValue(n).trim();
}
function label(s) { return ALIAS[s] || s; }
function unitOf(t) {
  const row = DATA.population.find(p => p.test === t);
  return row ? row.unit : "";
}
function fmt(v) { return (v === null || v === undefined) ? "—" : v.toFixed(2); }
function num(v) { return v.toLocaleString("en-US"); }
function el(tag, attrs, text) {
  const n = document.createElementNS(SVGNS, tag);
  for (const k in attrs) n.setAttribute(k, attrs[k]);
  if (text !== undefined) n.textContent = text;
  return n;
}
function tag(name, cls, text) {
  const n = document.createElement(name);
  if (cls) n.className = cls;
  if (text !== undefined) n.textContent = text;
  return n;
}
function nice(lo, hi) {
  const span = hi - lo, step = Math.pow(10, Math.floor(Math.log10(span / 4)));
  const mults = [1, 2, 2.5, 5, 10];
  let s = step;
  for (const m of mults) { if (span / (step * m) <= 6) { s = step * m; break; } }
  const ticks = [];
  for (let v = Math.ceil(lo / s) * s; v <= hi + 1e-9; v += s) ticks.push(+v.toFixed(6));
  return ticks;
}
// Analyte tabs: one panel per analyte, only the selected one shown.
function tabbed(host, tests, render) {
  const bar = tag("div", "tabs");
  bar.setAttribute("role", "tablist");
  const buttons = [], panels = [];
  tests.forEach((t, i) => {
    const b = tag("button"); b.type = "button"; b.textContent = t;
    b.setAttribute("role", "tab");
    b.setAttribute("aria-selected", i === 0 ? "true" : "false");
    bar.appendChild(b); buttons.push(b);
    const p = tag("div", "tabpanel");
    p.setAttribute("role", "tabpanel");
    if (i > 0) p.hidden = true;
    render(p, t);
    panels.push(p);
  });
  buttons.forEach((b, i) => b.addEventListener("click", () => {
    buttons.forEach((x, j) => x.setAttribute("aria-selected", j === i ? "true" : "false"));
    panels.forEach((p, j) => { p.hidden = j !== i; });
  }));
  host.appendChild(bar);
  panels.forEach(p => host.appendChild(p));
}
function figure(panel, svg, captionText, extra) {
  const fig = tag("figure");
  fig.appendChild(svg);
  if (extra) fig.appendChild(extra);
  fig.appendChild(tag("figcaption", null, captionText));
  panel.appendChild(fig);
}
function axisLabels(svg, ages, x, yBase, W, padL, padR) {
  for (const a of ages) {
    if (a % 2 === 1 || a === ages[ages.length - 1])
      svg.appendChild(el("text", { x: x(a), y: yBase + 16, "text-anchor": "middle",
                                   "font-size": 11 }, a));
  }
  svg.appendChild(el("text", { x: (padL + W - padR) / 2, y: yBase + 34,
                               "text-anchor": "middle", "font-size": 11 }, "Age, years"));
}

/* ---------- header, provenance ---------- */
(function head() {
  const box = document.getElementById("stats");
  const total = DATA.population.reduce((a, p) => a + p.n, 0);
  const items = DATA.population.map(p => [num(p.n), p.test]);
  items.push([num(total), "results in total"]);
  items.push(["1 to <18", "years of age"]);
  for (const [v, k] of items) {
    const d = tag("div", "stat");
    d.appendChild(tag("b", null, v));
    d.appendChild(tag("span", null, k));
    box.appendChild(d);
  }
  document.getElementById("built").textContent = DATA.built;
  const tb = document.querySelector("#prov tbody");
  for (const r of DATA.provenance) {
    const tr = tag("tr");
    tr.appendChild(tag("td", null, r.panel));
    const td = tag("td"); td.appendChild(tag("code", null, r.file));
    tr.appendChild(td); tb.appendChild(tr);
  }
})();

/* ---------- 1. classifier ---------- */
(function classifier() {
  const selTest = document.getElementById("c-test");
  const age = document.getElementById("c-age");
  const ageLab = document.getElementById("c-agelab");
  const val = document.getElementById("c-val");
  const unit = document.getElementById("c-unit");
  const body = document.querySelector("#c-table tbody");
  const tally = document.getElementById("c-tally");
  const hint = document.getElementById("c-hint");

  for (const t of TESTS) {
    const o = tag("option", null, t); o.value = t; selTest.appendChild(o);
  }
  // A starting value near the middle of the observed range for each analyte.
  const SEED = { PT: 10.5, aPTT: 32, Fibrinogen: 3.2 };
  let lastRows = [];

  function seedValue() {
    const t = selTest.value;
    unit.textContent = "(" + unitOf(t) + ")";
    val.step = t === "Fibrinogen" ? 0.05 : 0.1;
    val.value = SEED[t];
  }

  function render(pushState) {
    const t = selTest.value, a = +age.value, v = parseFloat(val.value);
    ageLab.textContent = a + (a === 1 ? " year" : " years");
    const rows = DATA.limits.filter(r => r.test === t && r.age === a)
      .sort((x, y) => ORDER.indexOf(label(x.method)) - ORDER.indexOf(label(y.method)));
    lastRows = rows;
    body.textContent = "";
    let outside = 0, total = 0;
    for (const r of rows) {
      const nm = label(r.method);
      const tr = tag("tr");
      const c1 = tag("td");
      const sw = tag("span", "swatch");
      sw.style.background = COLORS[nm] || "var(--muted)";
      c1.appendChild(sw); c1.appendChild(document.createTextNode(nm));
      const c2 = tag("td", "num", fmt(r.lower));
      const c3 = tag("td", "num", fmt(r.upper));
      const c4 = tag("td");
      if (Number.isFinite(v)) {
        const out = v < r.lower || v > r.upper;
        const s = tag("span", "verdict " + (out ? "out" : "in"),
          out ? (v < r.lower ? "Below reference limit" : "Above reference limit")
              : "Within reference limits");
        c4.appendChild(s);
        outside += out ? 1 : 0; total++;
      } else {
        c4.textContent = "—";
      }
      tr.appendChild(c1); tr.appendChild(c2); tr.appendChild(c3); tr.appendChild(c4);
      body.appendChild(tr);
    }
    if (!Number.isFinite(v)) {
      tally.textContent = "Enter a value to compare the sources.";
    } else {
      const vs = fmt(v) + " " + unitOf(t) + " at " + a + (a === 1 ? " year" : " years");
      if (outside === 0)
        tally.textContent = vs + " is within all " + total + " reference intervals.";
      else if (outside === total)
        tally.textContent = vs + " is outside all " + total + " reference intervals.";
      else
        tally.textContent = vs + " is flagged by " + outside + " of " + total +
          " reference interval sources and not by the other " + (total - outside) + ".";
    }
    if (pushState !== false && window.history && history.replaceState) {
      const p = new URLSearchParams({ test: t, age: String(a), value: val.value });
      history.replaceState(null, "", "?" + p.toString() + location.hash);
    }
  }

  function asText() {
    const t = selTest.value, a = +age.value, v = parseFloat(val.value);
    const lines = [`${t} ${fmt(v)} ${unitOf(t)} at age ${a}`, ""];
    for (const r of lastRows) {
      const out = v < r.lower || v > r.upper;
      lines.push(`${label(r.method).padEnd(24)} ${fmt(r.lower)} - ${fmt(r.upper)}   ` +
                 (out ? "OUTSIDE" : "within"));
    }
    lines.push("", tally.textContent);
    return lines.join("\n");
  }
  function flash(msg) {
    hint.textContent = msg;
    setTimeout(() => { hint.textContent = ""; }, 2500);
  }
  function copy(text, msg) {
    if (navigator.clipboard && navigator.clipboard.writeText)
      navigator.clipboard.writeText(text).then(() => flash(msg),
                                               () => flash("Copy blocked by the browser."));
    else flash("Copying is not available in this browser.");
  }

  // restore state from the address bar
  const q = new URLSearchParams(location.search);
  selTest.value = TESTS.includes(q.get("test")) ? q.get("test") : "PT";
  seedValue();
  const qa = parseInt(q.get("age"), 10);
  if (Number.isFinite(qa) && qa >= 1 && qa <= 17) age.value = String(qa);
  const qv = parseFloat(q.get("value"));
  if (Number.isFinite(qv)) val.value = String(qv);

  selTest.addEventListener("change", () => { seedValue(); render(); });
  age.addEventListener("input", () => render());
  val.addEventListener("input", () => render());
  document.getElementById("c-copy").addEventListener("click",
    () => copy(asText(), "Comparison copied."));
  document.getElementById("c-link").addEventListener("click",
    () => { render(); copy(location.href, "Link copied."); });
  render(false);
})();

/* ---------- 2. reference intervals ---------- */
(function intervals() {
  tabbed(document.getElementById("intervals-charts"), TESTS, (panel, t) => {
    const rows = DATA.intervals.filter(r => r.test === t)
      .sort((x, y) => ORDER.indexOf(label(x.method)) - ORDER.indexOf(label(y.method)));
    const W = 760, rowH = 34, padT = 26, padB = 34, padL = 168, padR = 22;
    const H = padT + rows.length * rowH + padB;
    let lo = Infinity, hi = -Infinity;
    for (const r of rows) {
      lo = Math.min(lo, r.lower_ci_lo === null || r.lower_ci_lo === undefined ? r.lower : r.lower_ci_lo);
      hi = Math.max(hi, r.upper_ci_hi === null || r.upper_ci_hi === undefined ? r.upper : r.upper_ci_hi);
    }
    const pad = (hi - lo) * 0.06; lo -= pad; hi += pad;
    const x = v => padL + (v - lo) / (hi - lo) * (W - padL - padR);

    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
                            "aria-label": t + " reference intervals by source" });
    const g = el("g", { class: "grid" });
    for (const v of nice(lo, hi)) {
      g.appendChild(el("line", { x1: x(v), x2: x(v), y1: padT - 8, y2: H - padB }));
      svg.appendChild(el("text", { x: x(v), y: H - padB + 16, "text-anchor": "middle",
                                   "font-size": 11 }, v));
    }
    svg.appendChild(g);
    svg.appendChild(el("text", { x: (padL + W - padR) / 2, y: H - 4,
                                 "text-anchor": "middle", "font-size": 11 },
                       t + ", " + unitOf(t)));
    rows.forEach((r, i) => {
      const y = padT + i * rowH + rowH / 2;
      const nm = label(r.method);
      const col = COLORS[nm] || "var(--muted)";
      const name = nm + (r.age_group && r.age_group !== "All ages" ? "  " + r.age_group : "");
      svg.appendChild(el("text", { x: padL - 12, y: y + 4, "text-anchor": "end",
                                   "font-size": 12 }, name));
      svg.appendChild(el("rect", { x: x(r.lower), y: y - 7,
                                   width: Math.max(1, x(r.upper) - x(r.lower)),
                                   height: 14, rx: 2, fill: col, opacity: .82 }));
      for (const [a, b] of [[r.lower_ci_lo, r.lower_ci_hi], [r.upper_ci_lo, r.upper_ci_hi]]) {
        if (a === null || a === undefined || b === null || b === undefined) continue;
        svg.appendChild(el("line", { x1: x(a), x2: x(b), y1: y, y2: y,
                                     stroke: "var(--ink)", "stroke-width": 1.2 }));
        for (const e of [a, b])
          svg.appendChild(el("line", { x1: x(e), x2: x(e), y1: y - 4, y2: y + 4,
                                       stroke: "var(--ink)", "stroke-width": 1.2 }));
      }
      if (r.n) svg.appendChild(el("text", { x: W - padR, y: y - 10, "text-anchor": "end",
                                            "font-size": 10 }, "n = " + num(r.n)));
    });
    figure(panel, svg,
      "Whiskers show the 90% confidence interval on each reference limit. The " +
      "manufacturer interval is an adult interval supplied in the package insert and " +
      "carries no confidence interval; the continuous row is the envelope of the " +
      "age-resolved curves over the ages the partition covers.");
  });

  document.getElementById("dl-limits").addEventListener("click", () => {
    const head = ["test", "age_years", "source", "lower", "upper", "unit"];
    const lines = [head.join(",")];
    for (const r of DATA.limits)
      lines.push([r.test, r.age, '"' + label(r.method) + '"', r.lower, r.upper,
                  unitOf(r.test)].join(","));
    const blob = new Blob([lines.join("\n") + "\n"], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = tag("a");
    a.href = url; a.download = "pediatric_coagulation_reference_limits.csv";
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });
})();

/* ---------- 3. continuous curves, with the sample size beneath ---------- */
(function curves() {
  tabbed(document.getElementById("curve-charts"), TESTS, (panel, t) => {
    const rows = DATA.curves.filter(r => r.test === t).sort((a, b) => a.age - b.age);
    const W = 760, H = 300, padT = 16, padB = 42, padL = 54, padR = 74;
    let lo = Infinity, hi = -Infinity;
    for (const r of rows) { lo = Math.min(lo, r.p2_5_lo); hi = Math.max(hi, r.p97_5_hi); }
    const pad = (hi - lo) * 0.08; lo -= pad; hi += pad;
    const ages = DATA.ages;
    const x = a => padL + (a - ages[0]) / (ages[ages.length - 1] - ages[0]) * (W - padL - padR);
    const y = v => H - padB - (v - lo) / (hi - lo) * (H - padT - padB);

    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
                            "aria-label": t + " continuous reference interval curves" });
    const grid = el("g", { class: "grid" });
    for (const v of nice(lo, hi)) {
      grid.appendChild(el("line", { x1: padL, x2: W - padR, y1: y(v), y2: y(v) }));
      svg.appendChild(el("text", { x: padL - 8, y: y(v) + 4, "text-anchor": "end", "font-size": 11 }, v));
    }
    svg.appendChild(grid);
    axisLabels(svg, ages, x, H - padB, W, padL, padR);
    svg.appendChild(el("text", { x: 14, y: (padT + H - padB) / 2, "font-size": 11,
                                 "text-anchor": "middle",
                                 transform: `rotate(-90 14 ${(padT + H - padB) / 2})` },
                       t + ", " + unitOf(t)));

    const col = COLORS["Continuous (GAMLSS)"];
    const series = [["p2_5", "p2_5_lo", "p2_5_hi", "2.5th", 1.6],
                    ["p50",  "p50_lo",  "p50_hi",  "50th",  2.4],
                    ["p97_5","p97_5_lo","p97_5_hi","97.5th",1.6]];
    for (const [k, klo, khi, name, w] of series) {
      const up = rows.map(r => `${x(r.age)},${y(r[khi])}`).join(" ");
      const dn = rows.slice().reverse().map(r => `${x(r.age)},${y(r[klo])}`).join(" ");
      svg.appendChild(el("polygon", { points: up + " " + dn, fill: col, opacity: .13 }));
      svg.appendChild(el("polyline", { points: rows.map(r => `${x(r.age)},${y(r[k])}`).join(" "),
                                       fill: "none", stroke: col, "stroke-width": w,
                                       "stroke-linejoin": "round" }));
      const last = rows[rows.length - 1];
      svg.appendChild(el("text", { x: x(last.age) + 8, y: y(last[k]) + 4, "font-size": 11,
                                   fill: col }, name));
    }
    const hair = el("line", { y1: padT, y2: H - padB, stroke: "var(--rule)",
                              "stroke-width": 1, opacity: 0 });
    svg.appendChild(hair);
    const read = el("text", { x: padL, y: padT + 2, "font-size": 11, opacity: 0 });
    svg.appendChild(read);
    const hit = el("rect", { x: padL, y: padT, width: W - padL - padR,
                             height: H - padT - padB, fill: "transparent" });
    svg.appendChild(hit);
    function move(ev) {
      const box = svg.getBoundingClientRect();
      const px = (ev.touches ? ev.touches[0].clientX : ev.clientX) - box.left;
      const frac = (px / box.width * W - padL) / (W - padL - padR);
      const a = Math.min(ages[ages.length - 1], Math.max(ages[0],
                Math.round(ages[0] + frac * (ages[ages.length - 1] - ages[0]))));
      const r = rows.find(q => q.age === a);
      if (!r) return;
      hair.setAttribute("x1", x(a)); hair.setAttribute("x2", x(a));
      hair.setAttribute("opacity", 1);
      read.setAttribute("opacity", 1);
      read.textContent = `age ${a}:  ${fmt(r.p2_5)} – ${fmt(r.p97_5)}  (median ${fmt(r.p50)})`;
    }
    hit.addEventListener("mousemove", move);
    hit.addEventListener("touchmove", move, { passive: true });
    hit.addEventListener("mouseleave", () => {
      hair.setAttribute("opacity", 0); read.setAttribute("opacity", 0);
    });
    figure(panel, svg,
      "Shaded bands are pointwise 90% confidence intervals from the bootstrap " +
      "implemented in the published GAMLSS pipeline.");

    // sample size strip, which is what the band width tracks
    const ss = DATA.sampleSize.filter(r => r.test === t).sort((a, b) => a.age - b.age);
    const H2 = 96, p2T = 14, p2B = 40;
    const maxN = Math.max(...ss.map(r => r.n));
    const bw = (W - padL - padR) / ages.length * 0.62;
    const svg2 = el("svg", { viewBox: `0 0 ${W} ${H2}`, role: "img",
                             "aria-label": t + " number of results by age" });
    for (const r of ss) {
      const h = (r.n / maxN) * (H2 - p2T - p2B);
      svg2.appendChild(el("rect", { x: x(r.age) - bw / 2, y: H2 - p2B - h,
                                    width: bw, height: Math.max(1, h), rx: 1,
                                    fill: "var(--muted)", opacity: .5 }));
    }
    svg2.appendChild(el("text", { x: padL - 8, y: p2T + 8, "text-anchor": "end",
                                  "font-size": 11 }, num(maxN)));
    svg2.appendChild(el("line", { x1: padL, x2: W - padR, y1: H2 - p2B, y2: H2 - p2B,
                                  stroke: "var(--rule)", "stroke-width": 1 }));
    axisLabels(svg2, ages, x, H2 - p2B, W, padL, padR);
    const total = ss.reduce((a, r) => a + r.n, 0);
    figure(panel, svg2,
      "Number of results at each single year of age, " + num(total) + " in total. " +
      "The confidence bands above widen where the data thin out.");
  });
})();

/* ---------- 4. flagging ---------- */
(function flagging() {
  tabbed(document.getElementById("flag-charts"), TESTS, (panel, t) => {
    const rows = DATA.flagAge.filter(r => r.test === t);
    const names = ORDER.filter(n => rows.some(r => label(r.ri_source) === n));
    const W = 760, H = 290, padT = 16, padB = 42, padL = 48, padR = 96;
    let hi = 0;
    for (const r of rows) hi = Math.max(hi, r.fp_rate);
    hi = Math.ceil(hi / 5) * 5;
    const ages = DATA.ages;
    const x = a => padL + (a - ages[0]) / (ages[ages.length - 1] - ages[0]) * (W - padL - padR);
    const y = v => H - padB - v / hi * (H - padT - padB);

    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
                            "aria-label": t + " out-of-RI flagging rate by age" });
    const grid = el("g", { class: "grid" });
    for (const v of nice(0, hi)) {
      grid.appendChild(el("line", { x1: padL, x2: W - padR, y1: y(v), y2: y(v) }));
      svg.appendChild(el("text", { x: padL - 8, y: y(v) + 4, "text-anchor": "end",
                                   "font-size": 11 }, v + "%"));
    }
    svg.appendChild(grid);
    axisLabels(svg, ages, x, H - padB, W, padL, padR);
    svg.appendChild(el("line", { x1: padL, x2: W - padR, y1: y(5), y2: y(5),
                                 stroke: "var(--ink)", "stroke-width": 1,
                                 "stroke-dasharray": "2 4", opacity: .55 }));
    svg.appendChild(el("text", { x: W - padR + 6, y: y(5) + 4, "font-size": 10 }, "5% nominal"));
    for (const nm of names) {
      const s = rows.filter(r => label(r.ri_source) === nm).sort((a, b) => a.age - b.age);
      svg.appendChild(el("polyline", { points: s.map(r => `${x(r.age)},${y(r.fp_rate)}`).join(" "),
                                       fill: "none", stroke: COLORS[nm], "stroke-width": 2,
                                       "stroke-linejoin": "round" }));
      const last = s[s.length - 1];
      svg.appendChild(el("circle", { cx: x(last.age), cy: y(last.fp_rate), r: 2.6,
                                     fill: COLORS[nm] }));
    }

    const tbl = tag("table");
    const thead = tag("thead");
    const hr = tag("tr");
    for (const [h, c] of [["Reference interval source", ""], ["Flagged", "num"],
                          ["of total", "num"], ["Rate (95% CI)", "num"]]) {
      const th = tag("th", c, h); hr.appendChild(th);
    }
    thead.appendChild(hr); tbl.appendChild(thead);
    const tb = tag("tbody");
    for (const nm of names) {
      const o = DATA.flagOverall.find(r => r.test === t && label(r.ri_source) === nm);
      const c = DATA.flagCounts.find(r => r.test === t && label(r.ri_source) === nm);
      const tr = tag("tr");
      const c1 = tag("td");
      const sw = tag("span", "swatch"); sw.style.background = COLORS[nm];
      c1.appendChild(sw); c1.appendChild(document.createTextNode(nm));
      tr.appendChild(c1);
      tr.appendChild(tag("td", "num", c ? num(c.count) : "—"));
      tr.appendChild(tag("td", "num", c ? num(c.n_total) : "—"));
      tr.appendChild(tag("td", "num", o ? o.fp_label : "—"));
      tb.appendChild(tr);
    }
    tbl.appendChild(tb);
    const wrap = tag("div", "scroll"); wrap.appendChild(tbl);
    figure(panel, svg,
      "Counts are the number of results outside the limits over the whole age range. " +
      "Rates carry 95% bootstrap confidence intervals from 1,000 resamples." +
      (t === "Fibrinogen"
        ? " Guven et al. did not report fibrinogen, so it is compared across the remaining sources."
        : ""), wrap);
  });
})();

/* ---------- 5. partitioning ---------- */
(function partitioning() {
  const box = document.getElementById("hb-box");
  const p = tag("p", "sub",
    "Harris–Boyd statistic at the candidate cut-point of 12 years, computed on " +
    "the full data set. Partitioning is taken as justified above a threshold of 3.");
  box.appendChild(p);
  const tbl = tag("table");
  const thead = tag("thead"); const hr = tag("tr");
  for (const [h, c] of [["Analyte", ""], ["Cut-point", "num"], ["n below", "num"],
                        ["n above", "num"], ["z", "num"], ["Partitioned", ""]])
    hr.appendChild(tag("th", c, h));
  thead.appendChild(hr); tbl.appendChild(thead);
  const tb = tag("tbody");
  for (const r of DATA.harrisBoyd) {
    const tr = tag("tr");
    tr.appendChild(tag("td", null, r.test));
    tr.appendChild(tag("td", "num", r.cutpoint + " y"));
    tr.appendChild(tag("td", "num", num(r.n_young)));
    tr.appendChild(tag("td", "num", num(r.n_old)));
    tr.appendChild(tag("td", "num", r.z_score.toFixed(2)));
    tr.appendChild(tag("td", null, r.justified ? "Yes" : "No"));
    tb.appendChild(tr);
  }
  tbl.appendChild(tb);
  const wrap = tag("div", "scroll"); wrap.appendChild(tbl);
  box.appendChild(wrap);

  tabbed(document.getElementById("partition-charts"), TESTS, (panel, t) => {
    const rows = DATA.sexLimits.filter(r => r.test === t);
    const bands = DATA.bands.filter(r => r.test === t);
    const groups = ["Female", "Male", "Combined"].filter(g => rows.some(r => r.sex === g));
    const W = 760, H = 330, padT = 16, padB = 42, padL = 54, padR = 96;
    let lo = Infinity, hi = -Infinity;
    for (const r of rows) {
      lo = Math.min(lo, r.lower_lo); hi = Math.max(hi, r.upper_hi);
    }
    for (const b of bands) {
      lo = Math.min(lo, b.lower_lo); hi = Math.max(hi, b.upper_hi);
    }
    const pad = (hi - lo) * 0.06; lo -= pad; hi += pad;
    const ages = DATA.ages;
    const x = a => padL + (a - ages[0]) / (ages[ages.length - 1] - ages[0]) * (W - padL - padR);
    const y = v => H - padB - (v - lo) / (hi - lo) * (H - padT - padB);
    const off = g => groups.length === 1 ? 0 : (g === "Female" ? -4 : 4);
    const half = (x(ages[1]) - x(ages[0])) / 2;

    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
                            "aria-label": t + " sex-specific reference limits by age" });
    const grid = el("g", { class: "grid" });
    for (const v of nice(lo, hi)) {
      grid.appendChild(el("line", { x1: padL, x2: W - padR, y1: y(v), y2: y(v) }));
      svg.appendChild(el("text", { x: padL - 8, y: y(v) + 4, "text-anchor": "end",
                                   "font-size": 11 }, v));
    }
    svg.appendChild(grid);

    // Sex-combined reference limits with their 90% confidence bands, one span
    // per reported partition, drawn behind the sex-specific points.
    const bcol = COLORS["Indirect (refineR)"];
    for (const b of bands) {
      const x0 = x(b.age_start) - half, x1 = x(b.age_end) + half;
      for (const [pt, clo, chi] of [[b.lower, b.lower_lo, b.lower_hi],
                                    [b.upper, b.upper_lo, b.upper_hi]]) {
        svg.appendChild(el("rect", { x: x0, y: y(chi), width: Math.max(1, x1 - x0),
                                     height: Math.max(1, y(clo) - y(chi)),
                                     fill: bcol, opacity: .16 }));
        svg.appendChild(el("line", { x1: x0, x2: x1, y1: y(pt), y2: y(pt),
                                     stroke: bcol, "stroke-width": 1.1, opacity: .85 }));
      }
      // label each partition's limits at the right edge of its own span
      const atEnd = b.age_end === ages[ages.length - 1];
      if (atEnd) {
        for (const pt of [b.upper, b.lower])
          svg.appendChild(el("text", { x: x1 + 6, y: y(pt) + 4, "font-size": 10,
                                       fill: bcol }, fmt(pt)));
      }
      svg.appendChild(el("text", { x: (x0 + x1) / 2, y: padT + 2,
                                   "text-anchor": "middle", "font-size": 10,
                                   fill: bcol }, b.age_group));
    }

    axisLabels(svg, ages, x, H - padB, W, padL, padR);
    svg.appendChild(el("text", { x: 14, y: (padT + H - padB) / 2, "font-size": 11,
                                 "text-anchor": "middle",
                                 transform: `rotate(-90 14 ${(padT + H - padB) / 2})` },
                       t + ", " + unitOf(t)));
    for (const g of groups) {
      const s = rows.filter(r => r.sex === g).sort((a, b) => a.age - b.age);
      const col = COLORS[g];
      for (const [k, klo, khi] of [["lower", "lower_lo", "lower_hi"],
                                   ["upper", "upper_lo", "upper_hi"]]) {
        for (const r of s) {
          const cx = x(r.age) + off(g);
          svg.appendChild(el("line", { x1: cx, x2: cx, y1: y(r[klo]), y2: y(r[khi]),
                                       stroke: col, "stroke-width": 1.1, opacity: .65 }));
          svg.appendChild(el("circle", { cx: cx, cy: y(r[k]), r: 2.3, fill: col }));
        }
      }
    }
    const leg = tag("div", "legend");
    for (const g of groups) {
      const s = tag("span");
      const sw = tag("span", "swatch"); sw.style.background = COLORS[g];
      s.appendChild(sw);
      s.appendChild(document.createTextNode(
        g === "Combined" ? "Both sexes combined" : g));
      leg.appendChild(s);
    }
    const sb = tag("span");
    const ssw = tag("span", "swatch"); ssw.style.background = bcol; ssw.style.opacity = ".55";
    sb.appendChild(ssw);
    sb.appendChild(document.createTextNode(
      "Reported interval and its 90% CI" +
      (bands.length > 1 ? ", one per age partition" : "")));
    leg.appendChild(sb);

    figure(panel, svg,
      "Points are the reference limits at each single year of age with their 90% " +
      "confidence intervals. The shaded spans are the intervals actually reported, " +
      (bands.length > 1
        ? "one for each age partition, split at 12 years. "
        : "estimated across the whole 1 to <18 year range. ") +
      (groups.length === 1
        ? "Fibrinogen was reported for both sexes combined, because single-year strata " +
          "held too few results to support stable sex-specific limits."
        : "Confidence intervals overlap between the sexes at the great majority of " +
          "ages, which is why a combined-sex interval is reported."), leg);
  });
})();

/* ---------- 6. transferability ---------- */
(function transferability() {
  const tests = TESTS.filter(t => DATA.transfer.some(r => r.test === t));
  tabbed(document.getElementById("transfer-charts"), tests, (panel, t) => {
    const rows = DATA.transfer.filter(r => r.test === t).sort((a, b) => a.age - b.age);
    const W = 760, H = 280, padT = 16, padB = 42, padL = 48, padR = 90;
    const ages = DATA.ages;
    const x = a => padL + (a - ages[0]) / (ages[ages.length - 1] - ages[0]) * (W - padL - padR);
    const y = v => H - padB - v / 100 * (H - padT - padB);
    const bw = (W - padL - padR) / ages.length * 0.6;

    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
                            "aria-label": t + " bootstrap pass rate by age" });
    const grid = el("g", { class: "grid" });
    for (const v of [0, 25, 50, 75, 100]) {
      grid.appendChild(el("line", { x1: padL, x2: W - padR, y1: y(v), y2: y(v) }));
      svg.appendChild(el("text", { x: padL - 8, y: y(v) + 4, "text-anchor": "end",
                                   "font-size": 11 }, v + "%"));
    }
    svg.appendChild(grid);
    axisLabels(svg, ages, x, H - padB, W, padL, padR);
    for (const r of rows) {
      const h = r.pass_rate / 100 * (H - padT - padB);
      svg.appendChild(el("rect", { x: x(r.age) - bw / 2, y: y(r.pass_rate),
                                   width: bw, height: Math.max(1, h), rx: 1,
                                   fill: COLORS["Guven et al."], opacity: .75 }));
    }
    svg.appendChild(el("line", { x1: padL, x2: W - padR, y1: y(90), y2: y(90),
                                 stroke: "var(--ink)", "stroke-width": 1,
                                 "stroke-dasharray": "3 3", opacity: .7 }));
    svg.appendChild(el("text", { x: W - padR + 6, y: y(90) + 4, "font-size": 10 },
                       "90% required"));
    const mean = rows.reduce((a, r) => a + r.pass_rate, 0) / rows.length;
    const passed = rows.filter(r => r.pass_rate >= 90).length;
    figure(panel, svg,
      `No age group reached the threshold: ${passed} of ${rows.length} passed, ` +
      `at a mean pass rate of ${mean.toFixed(1)}%. Each bar is the proportion of ` +
      `200 bootstrap samples of 20 results in which at most 2 fell outside the ` +
      `published interval for that age.`);
  });
})();

/* ---------- 7. published intervals ---------- */
(function published() {
  const tb = document.querySelector("#pub-table tbody");
  const rows = DATA.published.slice().sort((a, b) =>
    (a.study + a.test).localeCompare(b.study + b.test));
  for (const r of rows) {
    const tr = tag("tr");
    tr.appendChild(tag("td", null, r.study + (r.year ? " " + r.year : "") +
                       (r.country ? ", " + r.country : "")));
    tr.appendChild(tag("td", null, r.analyzer));
    tr.appendChild(tag("td", null, r.method));
    tr.appendChild(tag("td", null, r.test + "  " + r.age_group));
    tr.appendChild(tag("td", "num", r.n ? num(r.n) : "—"));
    tr.appendChild(tag("td", "num", String(r.lower)));
    tr.appendChild(tag("td", "num", String(r.upper)));
    tr.appendChild(tag("td", null, r.unit));
    tb.appendChild(tr);
  }
})();

/* ---------- sticky nav highlighting ---------- */
(function navHighlight() {
  if (!window.IntersectionObserver) return;
  const links = {};
  for (const a of document.querySelectorAll("#nav a"))
    links[a.getAttribute("href").slice(1)] = a;
  const seen = new Set();
  const io = new IntersectionObserver(entries => {
    for (const e of entries) {
      if (e.isIntersecting) seen.add(e.target.id); else seen.delete(e.target.id);
    }
    for (const id in links) links[id].classList.remove("here");
    for (const id in links) { if (seen.has(id)) { links[id].classList.add("here"); break; } }
  }, { rootMargin: "-10% 0px -80% 0px" });
  for (const id in links) {
    const s = document.getElementById(id);
    if (s) io.observe(s);
  }
})();
</script>
</body>
</html>
)---"

html <- sub("__PAYLOAD__", json, tmpl, fixed = TRUE)
if (!dir.exists("docs")) stop("docs/ not found; run from the project root.")
writeLines(html, "docs/index.html", useBytes = TRUE)

message(sprintf("[OK] docs/index.html written (%.0f KB).",
                file.size("docs/index.html") / 1024))
