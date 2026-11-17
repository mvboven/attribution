# Settings and entry point. From the project root: source("scripts/0_source.R").
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
# This file does not run MAP/HMC automatically.

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(project_dir, "attribution.Rproj")) &&
    file.exists(file.path(project_dir, "..", "attribution.Rproj"))) {
  project_dir <- normalizePath(file.path(project_dir, ".."), winslash = "/", mustWork = TRUE)
}
scripts_dir <- file.path(project_dir, "scripts")
if (!file.exists(file.path(project_dir, "attribution.Rproj")) ||
    !file.exists(file.path(scripts_dir, "analysis_helpers.R"))) {
  stop("Open attribution.Rproj, then source scripts/0_source.R from the project root.")
}

packages <- c("readr", "dplyr", "tidyr", "tibble", "ggplot2",
              "posterior", "cmdstanr", "loo", "cowplot", "withr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install the missing R packages: ", paste(missing_packages, collapse = ", "))
invisible(lapply(packages, function(package) suppressPackageStartupMessages(
  library(package, character.only = TRUE))))
source(file.path(scripts_dir, "analysis_helpers.R"), local = TRUE)

# Change analysis settings here, or set options(attribution.diagnosis_scope = ...)
# before sourcing this file. The public workflow reads the supplied compact CSVs.
diagnosis_scope <- getOption("attribution.diagnosis_scope", "main_diagnosis")
first_year <- 2013L
last_year <- 2019L
num_knots <- 200L
spline_degree <- 3L
min_lag <- -2L
max_lag <- 2L

# Notice that some prior distributions are still defined in the Stan code
prior_sigma <- rep(0.5, 8L)
prior_season_log_effect <- rep(0.38, 8L)
# Plotting parameters for the Gamma(2, 2) beta prior in the supplied Stan file.
# Editing these values changes the prior plot, not the prior hardcoded in Stan.
beta_prior_shape <- 2
beta_prior_rate <- 2

winter_start_week <- 44L
winter_end_week <- 13L
make_hmc_figure <- TRUE
make_parameter_figures <- TRUE
show_figures <- interactive()
run_loo <- TRUE
run_validation <- TRUE
bootstrap_replicates <- 2000L
bootstrap_block_length <- 13L

map_seed <- 1234L
hmc_seed <- 97L
hmc_init_seed <- 20260720L
n_chains <- 5L
iter_warmup <- 1000L
iter_sampling <- 1000L
thin <- 2L
adapt_delta <- 0.85
max_treedepth <- 12L
n_cores <- parallel::detectCores(logical = FALSE)
if (is.na(n_cores)) n_cores <- 1L
n_parallel_chains <- min(n_chains, n_cores)
loo_cores <- n_parallel_chains
options(cmdstanr_force_recompile = TRUE) # can be disabled; just to make sure

compact_data_dir <- file.path(project_dir, "data")
stan_file <- file.path(scripts_dir, "attribution.stan")
stan_build_dir <- Sys.getenv("ATTRIBUTION_STAN_BUILD_DIR", unset = "")
# Use the configured CmdStan installation, or set CMDSTAN to a local path.
cmdstan_directory <- Sys.getenv("CMDSTAN", unset = "")
if (nzchar(cmdstan_directory)) cmdstanr::set_cmdstan_path(cmdstan_directory)

stopifnot(diagnosis_scope %in% c("all_diagnoses", "main_diagnosis"),
          first_year == 2013L, last_year == 2019L,
          num_knots >= 2L, num_knots == as.integer(num_knots),
          min_lag < max_lag, length(prior_sigma) == 8L,
          length(prior_season_log_effect) == 8L, all(prior_sigma > 0),
          all(prior_season_log_effect > 0), beta_prior_shape > 0, beta_prior_rate > 0)
age_groups <- c("0-2m", "3-5m", "6-11m", "1-2y", "3-4y", "65-74y", "75-84y", "85+y")
pathogen_lookup <- tibble::tibble(
  pathogen_index = 1:4, pathogen_variable = c("VIRINF", "VIRRS", "VIRHMPV", "VIRPIV"),
  pathogen = c("Influenza", "RSV", "HMPV", "Parainfluenza"))
pathogen_vars <- pathogen_lookup$pathogen_variable
pathogen_table_order <- c("Influenza", "RSV", "HMPV", "Parainfluenza")
pathogen_colors <- c(Influenza = "#0072B2", RSV = "#D55E00", HMPV = "#009E73", Parainfluenza = "#CC79A7")
cbPalette <- c(Overall = "grey70", pathogen_colors)
num_types <- length(age_groups)
num_path <- nrow(pathogen_lookup)
season_labels <- paste0(2012:2019, "/", 2013:2020)
analysis_tag <- diagnosis_scope
model_id <- sprintf("hmc_%s_background_knots_%02d", analysis_tag, num_knots)
analysis_output_dir <- file.path(project_dir, "output", analysis_tag, model_id)
analysis_tmp_dir <- file.path(project_dir, "tmp", analysis_tag, model_id)
analysis_hmc_csv_dir <- file.path(analysis_output_dir, "hmc_csv")
analysis_loo_dir <- file.path(project_dir, "output", analysis_tag, "loo")
analysis_figure_dir <- file.path(project_dir, "figures", analysis_tag, model_id)
hmc_manifest_file <- file.path(analysis_output_dir, "hmc_manifest.rds")
map_init_file <- file.path(analysis_output_dir, "map_initial_values.rds")
for (directory in c(analysis_output_dir, analysis_tmp_dir, analysis_hmc_csv_dir,
                    analysis_loo_dir, analysis_figure_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
# Reset only analysis state, so a fit from another scope cannot be used silently.
fit_hmc <- NULL
fit_run_signature <- NULL
run_signature <- NULL
map_initial_values <- NULL
validation_run_signature <- NULL

# Numeric ordering avoids running 10_... before 2_.... Any stages can be selected.
analysis_scripts <- c("1_read_data_compact.R", "2_prepare_stan_data.R", "3_run_map.R", "4_run_hmc.R",
  "5_hmc_diagnostics.R", "6_figure_attribution.R", "7_figure_seasonal_multipliers.R",
  "8_figure_delay_distributions.R", "9_table_attributions.R", "10_figure_beta_prior_posterior.R",
  "11_calculate_loo.R", "12_compare_loo_models.R", "13_prepare_validation_data_compact.R",
  "14_validation_figure_correlation.R", "15_validation_figure_time_series.R")
run_analysis <- function(stages = seq_along(analysis_scripts), envir = .GlobalEnv) {
  stopifnot(length(stages) > 0L, all(stages %in% seq_along(analysis_scripts)),
            !anyDuplicated(stages))
  for (stage in stages) sys.source(file.path(scripts_dir, analysis_scripts[stage]), envir = envir)
  invisible(NULL)
}
message("Configured ", diagnosis_scope, ", ", first_year, "–", last_year,
        ", ", num_knots, " knots. Run run_analysis() or source selected stages.")

# Export all eight manuscript figures after both diagnosis groups are fitted.
run_manuscript_figures <- function() {
  source(file.path(scripts_dir, "16_figures_manuscript.R"))
}
