# One PSIS-LOO calculation, with scope-specific, observation-matched bundles.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
if (run_loo) {
  require_fit()
  message("Calculating PSIS-LOO for ", model_id, ".")
  loo_hmc <- fit_hmc$loo(variables = "log_lik", cores = loo_cores)
  pareto_k <- loo::pareto_k_values(loo_hmc)
  loo_summary <- tibble(
    model_id = model_id, diagnosis_scope = diagnosis_scope, num_knots = num_knots,
    observation_signature = observation_signature, run_signature = run_signature,
    elpd_loo = loo_hmc$estimates["elpd_loo", "Estimate"],
    se_elpd_loo = loo_hmc$estimates["elpd_loo", "SE"],
    p_loo = loo_hmc$estimates["p_loo", "Estimate"],
    looic = loo_hmc$estimates["looic", "Estimate"],
    se_looic = loo_hmc$estimates["looic", "SE"],
    maximum_pareto_k = max(pareto_k),
    n_pareto_k_gt_07 = sum(pareto_k > 0.7), n_pareto_k_gt_10 = sum(pareto_k > 1.0))
  loo_bundle <- list(model_id = model_id, num_knots = num_knots, diagnosis_scope = diagnosis_scope,
                     observation_signature = observation_signature, run_signature = run_signature, stan_md5 = stan_md5,
                     week_seq = week_seq, age_groups = age_groups,
                     loo = loo_hmc, summary = loo_summary)
  saveRDS(loo_bundle, file.path(analysis_loo_dir, paste0(model_id, "_loo.rds")))
  export_table(loo_summary, "loo_model_summary", analysis_output_dir)
  export_table(tibble(observation_index = seq_along(pareto_k), pareto_k = pareto_k),
               "loo_pareto_k", analysis_output_dir)
  message("ELPD-LOO ", round(loo_summary$elpd_loo, 2), " (SE ", round(loo_summary$se_elpd_loo, 2),
          "); LOOIC ", round(loo_summary$looic, 2), ". Pareto k > 1.0: ",
          loo_summary$n_pareto_k_gt_10, ".")
  if (any(!is.finite(pareto_k)) || any(pareto_k > 1.0)) {
    warning("PSIS-LOO has problematic Pareto k values; inspect them before comparing models.")
  }
} else message("PSIS-LOO skipped.")
