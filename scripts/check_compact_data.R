# Check both public data scopes and prepare their Stan inputs without fitting.
# Run from the project root: source("scripts/check_compact_data.R").
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
(function() {
  old_options <- options("attribution.diagnosis_scope")
  on.exit(options(old_options), add = TRUE)
  totals <- c(all_diagnoses = 515069, main_diagnosis = 289355)
  virology <- list()
  for (scope in names(totals)) {
    options(attribution.diagnosis_scope = scope)
    env <- new.env(parent = .GlobalEnv)
    source("scripts/0_source.R", local = env)
    env$run_analysis(stages = c(1, 2), envir = env)
    stopifnot(nrow(env$data) == 365L, nrow(env$y_data) == 365L,
              ncol(env$y_data) == 8L, sum(env$y_data) == totals[[scope]],
              "2015-W53" %in% env$data$ISOweek,
              identical(env$data_diagnosis_scope, scope),
              is.null(env$fit_hmc), is.null(env$fit_run_signature))
    file <- file.path(env$compact_data_dir, paste0("validation_confirmed_", scope, ".csv"))
    confirmed <- env$read_compact_confirmed_weekly(file, env$week_seq, names(env$pathogen_colors))
    stopifnot(nrow(confirmed) == 1460L)
    virology[[scope]] <- env$data[, env$pathogen_vars]
    message("PASS: ", scope, "; 365 weeks, eight age groups, four pathogens; ",
            sum(env$y_data), " RTI diagnosis-category counts; 1,460 validation rows.")
  }
  stopifnot(identical(virology[[1]], virology[[2]]))
  message("PASS: virology is identical across diagnosis scopes. No MAP/HMC/LOO was run.")
})()
