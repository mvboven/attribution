# HMC sampling and restart manifest. Sourcing this stage runs a new fit.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
require_objects(c("data_list", "run_signature"),
                instruction = "Run run_analysis(c(1, 2)) after sourcing scripts/0_source.R first.")
if (!file.exists(map_init_file)) stop("No saved MAP values; run stage 3 first.")
map_bundle <- readRDS(map_init_file)
if (!identical(map_bundle$run_signature, run_signature)) stop("MAP values do not match the current model/data.")
map_initial_values <- map_bundle$initial_values
mod <- compile_model()
init_fun_hmc <- function(chain_id) withr::with_seed(hmc_init_seed + chain_id, {
  initial <- map_initial_values
  initial$beta <- pmax(initial$beta * exp(matrix(rnorm(num_types * num_path, 0, 0.05), num_types)), 1e-8)
  initial$alpha <- pmax(initial$alpha * exp(matrix(rnorm(num_types * num_path, 0, 0.05), num_types)), 1e-8)
  initial$mu <- pmin(pmax(initial$mu + matrix(rnorm(num_types * num_path, 0, 0.05), num_types),
                         min_lag + 1e-6), max_lag - 1e-6)
  initial$coefficients <- initial$coefficients + matrix(rnorm(num_types * num_basis, 0, 0.02), num_types)
  initial$background_intercept <- initial$background_intercept + rnorm(num_types, 0, 0.02)
  initial$season_log_effect <- initial$season_log_effect + matrix(rnorm(num_path * 8L, 0, 0.03), num_path)
  initial$season_log_effect <- sweep(initial$season_log_effect, 1,
                                     rowMeans(initial$season_log_effect), "-")
  initial
})
message("Sampling ", n_chains, " HMC chains; ", n_parallel_chains, " run simultaneously.")
fit_hmc <- mod$sample(
  data = data_list, init = init_fun_hmc, seed = hmc_seed,
  chains = n_chains, parallel_chains = n_parallel_chains,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, thin = thin,
  adapt_delta = adapt_delta, max_treedepth = max_treedepth,
  refresh = 200L, show_messages = TRUE, show_exceptions = FALSE,
  output_dir = analysis_hmc_csv_dir, output_basename = model_id,
  save_cmdstan_config = TRUE)
fit_run_signature <- run_signature
csv_paths <- normalizePath(fit_hmc$output_files(), winslash = "/", mustWork = TRUE)
prefix <- paste0(project_dir, "/")
if (!all(startsWith(csv_paths, prefix))) stop("HMC outputs must be inside the project folder.")
saveRDS(list(model_id = model_id, run_signature = run_signature,
             observation_signature = observation_signature, stan_md5 = stan_md5,
             csv_files = substring(csv_paths, nchar(prefix) + 1L),
             sampling = list(chains = n_chains, warmup = iter_warmup, sampling = iter_sampling,
                             thin = thin, adapt_delta = adapt_delta, max_treedepth = max_treedepth)),
        hmc_manifest_file)
writeLines(capture.output(sessionInfo()), file.path(analysis_output_dir, "session_info.txt"))
rm(map_bundle, map_initial_values)
message("HMC finished. The chain CSVs and run manifest are saved for reuse.")
