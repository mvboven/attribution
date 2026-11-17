# Compare saved LOO objects on the same observations, in the same order.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
if (run_loo) {
  require_objects("observation_signature", instruction = "Run run_analysis(c(1, 2)) after sourcing scripts/0_source.R first.")
  loo_files <- list.files(analysis_loo_dir, pattern = "_loo\\.rds$", full.names = TRUE)
  loo_bundles <- lapply(loo_files, readRDS)
  compatible <- vapply(loo_bundles, function(bundle) {
    is.list(bundle) && identical(bundle$diagnosis_scope, diagnosis_scope) &&
      identical(bundle$observation_signature, observation_signature) &&
      identical(bundle$week_seq, week_seq) && identical(bundle$age_groups, age_groups) &&
      identical(bundle$stan_md5, stan_md5) && !is.null(bundle$loo)
  }, logical(1))
  if (any(!compatible)) message("Skipped ", sum(!compatible), " LOO bundles with different or unverified observations/model source.")
  loo_bundles <- loo_bundles[compatible]
  if (length(loo_bundles)) {
    summaries <- bind_rows(lapply(loo_bundles, function(bundle) bundle$summary)) %>% arrange(num_knots)
    export_table(summaries, "loo_model_summary", dirname(analysis_loo_dir))
  }
  if (length(loo_bundles) >= 2L) {
    loo_list <- lapply(loo_bundles, function(bundle) bundle$loo)
    names(loo_list) <- vapply(loo_bundles, function(bundle) bundle$model_id, character(1))
    if (anyDuplicated(names(loo_list))) stop("Duplicate LOO model IDs; retain only one bundle per model.")
    loo_comparison <- loo::loo_compare(loo_list)
    comparison_table <- tibble::as_tibble(as.data.frame(loo_comparison), rownames = "model_id")
    export_table(comparison_table, "loo_model_comparison", dirname(analysis_loo_dir))
    message("Compared ", length(loo_list), " models. Best ELPD-LOO: ", rownames(loo_comparison)[1],
            ". Paired differences and their SEs are saved as CSV.")
  } else message("LOO comparison needs at least two saved models fitted to identical observations.")
} else message("LOO comparison skipped.")
