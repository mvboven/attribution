# Respiratory-virus attribution of hospital admissions

R and Stan code for attributing weekly respiratory-tract-infection (RTI) hospital-admission counts to influenza, RSV, human metapneumovirus (HMPV), and parainfluenza (PIV).

The analysis covers ISO years 2013–2019: 365 weeks, including 2015-W53, and eight age groups. Separate analyses use all diagnoses and the main diagnosis. The public workflow reads the four aggregate CSVs in `data/compact/`.

## Installation

Open `attribution.Rproj` in RStudio, or start R in the repository root. Install the R packages once:

```r
source("scripts/install_dependencies.R")
```

Model fitting also requires [CmdStan and its C++ toolchain](https://mc-stan.org/cmdstanr/articles/cmdstanr.html):

```r
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
```

For an existing CmdStan installation, set the `CMDSTAN` environment variable to its installation directory, or use `cmdstanr::set_cmdstan_path()`.

## Check the supplied data

```r
source("scripts/check_compact_data.R")
```

This checks both diagnosis scopes, constructs the Stan input lists, and writes the prepared inputs under `output/`. It does not fit a model or require a CmdStan installation; the `cmdstanr` R package must be installed.

## Run the analyses

Configure and run the all-diagnosis analysis:

```r
options(attribution.diagnosis_scope = "all_diagnoses")
source("scripts/0_source.R")
run_analysis()
```

Configure and run the main-diagnosis analysis:

```r
options(attribution.diagnosis_scope = "main_diagnosis")
source("scripts/0_source.R")
run_analysis()
```

Without the option, the default is `main_diagnosis`. Re-source `0_source.R` after changing the scope: this updates the file paths and resets the fit state. Other settings are defined in that file.

`0_source.R` only configures the analysis. `run_analysis()` runs the selected stages, including new MAP and HMC fits. The supplied settings use 200 knots, spline degree 3, lags −2 to 2 weeks, five chains, 1,000 warm-up and 1,000 sampling iterations per chain, and thinning by 2. Validation uses 2,000 circular block-bootstrap replicates with block length 13 weeks. The model, including its Gamma(2,2) prior for the scaling coefficients, is in `scripts/attribution_20260821.stan`.

To run stages separately:

```r
run_analysis(c(1, 2))  # Read compact data and prepare Stan inputs
run_analysis(c(3, 4))  # New MAP and HMC fits
run_analysis(5:15)    # Diagnostics, figures, tables, LOO and validation
```

To reuse saved chains in a later session, configure the same scope and model settings, run stages 1–2, then run stages 5–15. The helpers reload the saved fit and verify its data/model signature. Keep the saved chain CSVs, prepared-data bundle, and manifest together within the project.

## Manuscript figures

After fitting both diagnosis scopes:

```r
run_manuscript_figures()
```

Script 16 reads the saved fits and the compact validation counts and exports eight figures as PDF and PNG. Its default is 200 knots; edit its `num_knots` setting if different fits are being used. No model is refitted by script 16.

## Repository contents

| Location | Contents |
|---|---|
| `data/compact/` | Four public input CSVs, their data dictionary and SHA-256 checksums |
| `scripts/0_source.R` | Configuration and ordered analysis entry point |
| `scripts/1_read_data_compact.R` | Public admission/virology reader, used by stage 1 |
| `scripts/2_prepare_stan_data.R` through `12_compare_loo_models.R` | Model preparation, fitting, diagnostics, attribution and model comparisons |
| `scripts/13_prepare_validation_data_compact.R` | Public validation reader and bootstrap statistics, used by stage 13 |
| `scripts/14_validation_figure_correlation.R`, `15_validation_figure_time_series.R` | Validation figures |
| `scripts/16_figures_manuscript.R` | Combined manuscript figure export from saved fits |
| `scripts/analysis_helpers.R` | Shared input checks, posterior summaries and restart helpers |
| `scripts/attribution_20260821.stan` | Supplied Stan model |

The original readers `1_read_data.R` and `13_prepare_validation_data.R` are retained to document the compact-data exports. They require the original source files and are not selected by the public entry point. The default analysis does not need those source files.

The additional standalone scripts `17_figure_PIV_types_monthly.R` and `18_figure_linked_time.R` require the original virological workbook. Script 17 uses PIV subtypes; script 18 uses the workbook's full available period. Those inputs are not recoverable from the four compact CSVs and are not included. These two optional scripts are not called by `run_analysis()` or `run_manuscript_figures()`. The earlier `6_figure_attribution_old.R` is retained as an alternative plot and is not selected by the entry point.

## Outputs and provenance

Analysis outputs are written to `output/<scope>/`, figures to `figures/<scope>/`, and temporary fit files to `tmp/<scope>/`. These generated directories are excluded from version control. Posterior samples and fitted results are generated by running the analyses; they are not included in this repository.

The compact CSVs preserve the supplied preprocessing of suppressed counts, age-group aggregation, and absent hospital week/age cells. The admission counts are sums of RTI diagnosis categories; coinfections may contribute to more than one category. The validation counts aggregate linked pathogen categories across the same eight age groups. See [the data dictionary](data/compact/README.md) for definitions.

The public packaging keeps the supplied Stan model unchanged. It selects the compact readers, uses project-relative compact paths, and updates script 16 to read compact validation counts. Script 10's prior plot uses Gamma(2,2) to match the supplied Stan model. See [VALIDATION.md](VALIDATION.md) for the checks performed.

## Licence

Code is distributed under the BSD 3-Clause licence; see [LICENSE](LICENSE).
