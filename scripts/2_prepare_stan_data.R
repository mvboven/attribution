# Prepare matrices, season indices, priors, and the Stan data list.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
require_objects(c("data", "age_groups", "pathogen_vars"),
                instruction = "Source 0_source.R and 1_read_data.R first.")
if (!identical(data_diagnosis_scope, diagnosis_scope)) stop("Diagnosis scope changed; reread the data.")

y_data <- as.matrix(data[, age_groups])  # week x age group
if (any(!is.finite(y_data)) || any(y_data < 0) || any(y_data != round(y_data)) ||
    any(y_data > .Machine$integer.max)) stop("Invalid hospital admissions in Stan input.")
storage.mode(y_data) <- "integer"
n_data <- nrow(y_data)
x_data <- seq_len(n_data)
num_basis <- num_knots + spline_degree - 1L
knots <- seq(1, n_data, length.out = num_knots)
VIR <- t(as.matrix(data[, pathogen_vars]))  # pathogen x week
rownames(VIR) <- pathogen_vars
if (any(!is.finite(VIR)) || any(VIR < 0)) stop("Virological counts are missing, non-finite, or negative.")

# Seasons run from week 27 through week 26 of the following ISO year.
data <- data %>% mutate(
  season_start_year = if_else(week >= 27L, year, year - 1L),
  season = paste0(season_start_year, "/", season_start_year + 1L),
  season_index = match(season_start_year, (first_year - 1L):last_year))
season_index <- data$season_index
week_seq <- as.Date(data$week_start)
stopifnot(n_data == 365L, "2015-W53" %in% data$ISOweek,
          all(diff(week_seq) == 7), !anyNA(season_index), all(season_index %in% 1:8),
          identical(pathogen_vars, pathogen_lookup$pathogen_variable))
data_list <- list(
  num_weeks = n_data, num_knots = num_knots, knots = knots,
  num_types = num_types, num_basis = num_basis, num_path = num_path,
  spline_degree = spline_degree, x_data = x_data, y_data = t(y_data), VIR = VIR,
  prior_sigma = prior_sigma, prior_season_log_effect = prior_season_log_effect,
  min_lag = min_lag, max_lag = max_lag, season_index = season_index)

# Fingerprints protect restarts and comparisons from mixing different observations.
observation_signature <- object_md5(list(diagnosis_scope, week_seq, age_groups, y_data))
stan_md5 <- if (file.exists(stan_file)) unname(tools::md5sum(stan_file)) else NA_character_
run_signature <- object_md5(list(data_list, stan_md5, age_groups, pathogen_lookup))
saveRDS(list(data = data, data_list = data_list, observation_signature = observation_signature,
             run_signature = run_signature, stan_md5 = stan_md5,
             age_groups = age_groups, pathogen_lookup = pathogen_lookup),
        file.path(analysis_output_dir, "prepared_stan_data.rds"))
message("Prepared Stan data: ", n_data, " weeks, ", num_types, " age groups, ", num_path,
        " pathogens, ", num_basis, " spline coefficients per age group.")
