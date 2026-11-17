# Sampler and parameter diagnostics, with detailed CSVs.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
require_fit()
diagnostic_summary <- fit_hmc$diagnostic_summary()
parameter_summary <- fit_hmc$summary(variables = c(
  "beta", "mu", "alpha", "coefficients", "background_intercept", "season_log_effect"))
finite_extreme <- function(x, fun) if (any(is.finite(x))) fun(x[is.finite(x)]) else NA_real_
diagnostic_overview <- parameter_summary %>% summarise(
  maximum_rhat = finite_extreme(rhat, max),
  minimum_ess_bulk = finite_extreme(ess_bulk, min),
  minimum_ess_tail = finite_extreme(ess_tail, min))
export_table(as_tibble(diagnostic_summary), "hmc_sampler_diagnostics", analysis_output_dir)
export_table(parameter_summary, "hmc_parameter_diagnostics", analysis_output_dir)
export_table(diagnostic_overview, "hmc_diagnostic_overview", analysis_output_dir)
message("HMC diagnostics: max R-hat ", round(diagnostic_overview$maximum_rhat, 3),
        "; minimum bulk/tail ESS ", round(diagnostic_overview$minimum_ess_bulk), "/",
        round(diagnostic_overview$minimum_ess_tail), ". Detailed diagnostics saved as CSV.")
if (is.na(diagnostic_overview$maximum_rhat) || diagnostic_overview$maximum_rhat > 1.01) {
  warning("R-hat is missing or above 1.01; inspect the diagnostic CSVs before interpreting estimates.")
}
if (any(diagnostic_summary$num_divergent > 0) || any(diagnostic_summary$num_max_treedepth > 0)) {
  warning("Divergences or maximum-treedepth events occurred; inspect the sampler diagnostics.")
}
