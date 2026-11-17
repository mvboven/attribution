# Respiratory virus attribution

R and Stan analysis of weekly respiratory-tract-infection (RTI) hospital-admission counts, 2013–2019, for eight age groups. Separate analyses use all diagnoses and the main diagnosis. The four pathogen series are influenza (A and B combined), respiratory syncytial virus (RSV), human metapneumovirus (HMPV) and parainfluenza (PIV, all types combined).

## Install

Open `attribution.Rproj` in RStudio. Install the R packages and CmdStan once:

```r
install.packages(c("readr", "dplyr", "tidyr", "tibble", "ggplot2",
                   "posterior", "loo", "cowplot", "withr"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev",
                                      "https://cloud.r-project.org"))
cmdstanr::install_cmdstan()
```

CmdStan requires a working C++ toolchain; see the [CmdStanR installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html). For an existing installation, use `cmdstanr::set_cmdstan_path()` or the `CMDSTAN` environment variable.

## Run

From the project root, run both diagnosis scopes:

```r
for (scope in c("all_diagnoses", "main_diagnosis")) {
  options(attribution.diagnosis_scope = scope)
  source("scripts/0_source.R")
  run_analysis()
}
run_manuscript_figures()
```

Settings are in `scripts/0_source.R`; the model is in `scripts/attribution.stan`. Sourcing the settings alone does not fit a model. `run_analysis()` runs new MAP/HMC fits, diagnostics, tables, figures, LOO and validation. Outputs are saved separately by scope in `output/` and `figures/`.

Stages can also be selected:

```r
run_analysis(1:2)   # Read data and prepare Stan inputs
run_analysis(3:4)   # Fit MAP and HMC
run_analysis(5:15)  # Diagnostics, tables, figures, LOO and validation
```

To reuse saved fits, configure the same scope and settings, then run stages 1–2 and 5–15. Keep the generated chain CSVs and manifests together in the project. Script 16 exports the eight manuscript figures from the saved fits for both scopes.

## Data

All input files are directly in `data/`.

| File | Contents |
|---|---|
| `attribution_all_diagnoses.csv` | Weekly admissions and virology, all diagnoses |
| `attribution_main_diagnosis.csv` | Weekly admissions and virology, main diagnosis |
| `validation_confirmed_all_diagnoses.csv` | Weekly confirmed pathogen counts, all diagnoses |
| `validation_confirmed_main_diagnosis.csv` | Weekly confirmed pathogen counts, main diagnosis |

Admission files have 365 rows. Columns are ISO `year` and `week`; eight age-specific count columns (`0-2m`, `3-5m`, `6-11m`, `1-2y`, `3-4y`, `65-74y`, `75-84y`, `85+y`); and four virological counts: `VIRINF` (influenza), `VIRRS`, `VIRHMPV`, `VIRPIV`. The age suffixes `m` and `y` mean months and years. Virology is shared between diagnosis scopes. Hospital counts sum RTI diagnosis categories, including the supplied handling of coinfections.

Validation files have 1,460 rows and three columns: `week` (ISO Monday), `component` (pathogen), and `confirmed` (count across the eight age groups). Fractional confirmed counts result from replacing suppressed cells ("1-4" instead of "1", "2","3", or "4") by their midpoint. Validation compares these counts with the sum of age-specific posterior median fitted pathogen intensities.

## Licence

Code: BSD 3-Clause; see [LICENSE](LICENSE).
