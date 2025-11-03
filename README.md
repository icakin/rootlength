
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rootlength

**Last updated:** 2025-11-03 15:42 GMT

<!-- badges: start -->

[![R-CMD-check](https://github.com/icakin/rootlength/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/icakin/rootlength/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/rootlength)](https://CRAN.R-project.org/package=rootlength)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/rootlength)](https://cranlogs.r-pkg.org/badges/rootlength)
<!-- badges: end -->

The goal of **rootlength** is to analyze seed/root-length bioassay
experiments from a simple CSV file and produce tidy tables and
publication-ready plots. It converts root-length measurements to
centimeters, groups by treatment × dose, and runs either Welch ANOVA +
Games–Howell (all-pairs) or Dunnett tests versus a user-defined control
(with HC3-robust fallback). It also computes relative germination,
growth, and Germination Index (GI) % change versus the chosen control.

## Installation

Once on CRAN, you will be able to install the released version with:

``` r
install.packages("rootlength")
```

You can install the development version of **rootlength** from GitHub
with:

``` r
# install.packages("remotes")
remotes::install_github("icakin/rootlength")
```

Then load the package:

``` r
library(rootlength)
```

## Input format

`run_pipeline()` expects a CSV file with **columns**, not a special file
format.

### Required columns

- `treatment`  
  Character; treatment / group name, e.g. `"Control"`, `"RB"`, `"MMB"`.

- `root_length`  
  Numeric; root length for each seed/root. Zero or `NA` can represent
  non-germination.

- `length_unit`  
  Character; length unit. Accepted values include:  
  `"cm"`, `"centimeter"`, `"centimeters"`,  
  `"mm"`, `"millimeter"`, `"millimeters"`,  
  `"m"`, `"meter"`, `"meters"`,  
  `"in"`, `"inch"`, `"inches"`,  
  `"ft"`, `"foot"`, `"feet"`.

### Optional columns

If present, they are used; if omitted, the analysis still works:

- `concentration`  
  Numeric dose (e.g. `0`, `1`, `4`).

- `dose_unit`  
  Character; dose unit, e.g. `"g/L"`.

A typical header might look like:

``` text
treatment,concentration,dose_unit,root_length,length_unit
```

## Example: from CSV to results

Assume you have a CSV file named `rootlength_data.csv` in your working
directory.

``` r
library(rootlength)

file    <- "rootlength_data.csv"   # path to your CSV
out_dir <- "rootlength_outputs"    # folder where results will be written
```

### 1. Explore all groups (Welch + Games–Howell)

First, explore all groups and see what labels exist:

``` r
res_explore <- run_pipeline(
  file    = file,
  out_dir = out_dir,
  mode    = "explore",   # ignore control, all-pairs comparisons
  boot    = FALSE,
  quiet   = FALSE        # show control-selection banner + messages
)

# Group-level summaries
head(res_explore$results)

# Games–Howell multiple comparisons
head(res_explore$games_howell)
```

`run_pipeline()` prints a banner listing detected groups, for example:

    ================ CONTROL SELECTION ================
    Steps:
      1) Review the detected group labels below.
      2) Choose mode & control_name via function arguments.
         - Exact label (recommended): e.g., RB (4 g/L)
         - Pooled treatment: e.g., RB
      Tip: mode = 'explore' ignores control_name.
    Detected options:
      - BMC (1 g/L)
      - BMC (4 g/L)
      - Control (0 g/L)
      - MMB (1 g/L)
      - MMB (4 g/L)
      - RB (1 g/L)
      - RB (4 g/L)
      - Spent Lees A
    ==================================================

### 2. Compare vs a chosen control (Dunnett + GI)

Pick one of these labels as the control, e.g. `"Control (0 g/L)"`, and
run:

``` r
res_control <- run_pipeline(
  file         = file,
  out_dir      = out_dir,
  mode         = "control",
  control_name = "Control (0 g/L)",  # must match one of the printed labels
  boot         = TRUE,               # TRUE = bootstrap GI % change CIs
  boot_B       = 2000,
  quiet        = FALSE
)

# Main table: germination, growth, GI, % change vs control
head(res_control$results)
```

You can also inspect the test tables:

``` r
# Dunnett vs control (standard and robust HC3)
head(res_control$dunnett_standard)
head(res_control$dunnett_hc3)

# Games–Howell results
head(res_control$games_howell)
```

All tables and plots are written to the folder you passed as `out_dir`
(here `"rootlength_outputs"`), including:

- CSV summaries (e.g. `rootlength_results_vs_control.csv`)
- Dunnett and Games–Howell tables
- Boxplots of root length by group
- GI bar plots and Dunnett-style forest plots

## Minimal toy example (no dose columns)

You can also run `rootlength` on a tiny dataset without doses:

``` r
library(rootlength)

dat <- data.frame(
  treatment   = c("Control", "Treatment"),
  root_length = c(2.1, 1.7),
  length_unit = c("cm", "cm")
)
tmp <- tempfile(fileext = ".csv")
write.csv(dat, tmp, row.names = FALSE)

run_pipeline(
  file    = tmp,
  mode    = "explore",
  boot    = FALSE,
  quiet   = TRUE
)
```
