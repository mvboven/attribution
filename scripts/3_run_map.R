# MAP fit and named initial estimates for HMC.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
require_objects(c("data_list", "run_signature"),
                instruction = "Run run_analysis(c(1, 2)) after sourcing scripts/0_source.R first.")
mod <- compile_model()
init_fun <- function() withr::with_seed(map_seed, list(
  coefficients = matrix(0, num_types, num_basis),
  background_intercept = log(pmax(colMeans(y_data), 1)),
  beta = matrix(0.002, num_types, num_path),
  mu = matrix((min_lag + max_lag) / 2, num_types, num_path),
  alpha = matrix(runif(num_types * num_path, 0.1, 2), num_types, num_path),
  season_log_effect = matrix(0, num_path, 8L)))
message("Finding the MAP estimate for HMC initialisation.")
fit_map <- mod$optimize(
  data = data_list, init = init_fun, seed = map_seed, iter = 10000L,
  algorithm = "lbfgs", jacobian = TRUE, refresh = 1000L,
  output_dir = analysis_tmp_dir, save_cmdstan_config = TRUE)
map_initial_values <- list(
  mu = indexed_point_estimate(fit_map$mle("mu"), "mu", c(num_types, num_path)),
  alpha = indexed_point_estimate(fit_map$mle("alpha"), "alpha", c(num_types, num_path)),
  beta = indexed_point_estimate(fit_map$mle("beta"), "beta", c(num_types, num_path)),
  coefficients = indexed_point_estimate(fit_map$mle("coefficients"), "coefficients", c(num_types, num_basis)),
  background_intercept = indexed_point_estimate(fit_map$mle("background_intercept"), "background_intercept", num_types),
  season_log_effect = indexed_point_estimate(fit_map$mle("season_log_effect"), "season_log_effect", c(num_path, 8L)))
saveRDS(list(run_signature = run_signature, initial_values = map_initial_values), map_init_file)
rm(fit_map)
message("Saved MAP initial values for all six parameter blocks.")
