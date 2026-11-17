# Align linked pathogen counts and fitted intensities; block-bootstrap correlations.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
if (run_validation) {
  require_fit()
  message("Preparing linked-pathogen validation and block-bootstrap intervals.")
  # Match linked validation to the selected hospital diagnosis scope.
  confirmed_weekly <- read_confirmed_weekly(
    linked_validation_file, age_group_lookup, week_seq, age_groups)
  
  # Same estimand as supplied code: sum the age-specific posterior median fitted
  # pathogen intensities. It is intentionally distinct from y * contribution/lambda.
  pathogen_summary <- summarise_pathogen_intensities(
    fit_hmc, week_seq, age_groups, pathogen_lookup$pathogen)
  attributed_weekly <- pathogen_summary %>%
    transmute(
      week = as.Date(week),
      age_group = as.character(age_group),
      component = as.character(component),
      attributed = as.numeric(.data$y_data)
    ) %>%
    filter(
      week %in% week_seq,
      age_group %in% age_groups,
      component %in% names(pathogen_colors)
    ) %>%
    group_by(week, component) %>%
    summarise(
      attributed = sum(attributed),
      .groups = "drop"
    )
  
  # Combine complete weekly time series.
  validation_weekly <- tidyr::expand_grid(
    week = week_seq,
    component = names(pathogen_colors)
  ) %>%
    left_join(confirmed_weekly, by = c("week", "component")) %>%
    mutate(confirmed = replace_na(confirmed, 0)) %>%
    left_join(attributed_weekly, by = c("week", "component")) %>%
    mutate(
      component = factor(component, levels = names(pathogen_colors))
    ) %>%
    arrange(component, week)
  
  stopifnot(
    nrow(validation_weekly) == length(week_seq) * length(pathogen_colors),
    !anyNA(validation_weekly$confirmed),
    !anyNA(validation_weekly$attributed),
    all(validation_weekly$confirmed >= 0),
    all(validation_weekly$attributed >= 0)
  )
  
  # Minimal public validation input: confirmed counts by ISO Monday and pathogen.
  # Counts already combine the selected age_groups and include explicit zeros.
  # Export each diagnosis_scope separately, using the same directory as script 1.
  # Attributed counts are recomputed from the current fit by the compact reader.
  if (!exists("compact_data_dir", inherits = FALSE)) {
    compact_data_dir <- file.path("data", "compact")
  }
  validation_compact_file <- file.path(
    compact_data_dir,
    paste0("validation_confirmed_", diagnosis_scope, ".csv")
  )
  validation_confirmed_compact <- validation_weekly %>%
    transmute(
      week = as.Date(week),
      component = as.character(component),
      confirmed = as.numeric(confirmed)
    )
  check_unique(validation_confirmed_compact, c("week", "component"), "Compact validation data")
  if (any(!is.finite(validation_confirmed_compact$confirmed)) ||
      any(validation_confirmed_compact$confirmed < 0)) {
    stop("Compact confirmed counts must be finite and non-negative.")
  }
  dir.create(compact_data_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(validation_confirmed_compact, validation_compact_file)
  message("Wrote ", nrow(validation_confirmed_compact),
          " weekly pathogen counts to ", validation_compact_file, ".")
  
  validation_statistics <- validation_weekly %>%
    group_by(component) %>%
    group_modify(~ bind_cols(
      tibble(n_weeks = nrow(.x),
             scaling_factor = if (sum(.x$attributed) > 0) sum(.x$confirmed) / sum(.x$attributed) else NA_real_),
      block_spearman(.x$confirmed, .x$attributed, n_boot = bootstrap_replicates,
                     block_length = bootstrap_block_length,
                     seed = 20260903L + match(as.character(.y$component), names(pathogen_colors))))) %>%
    ungroup()
  if (anyNA(validation_statistics$spearman_rho) || anyNA(validation_statistics$scaling_factor)) {
    warning("Some validation series are constant or have zero fitted totals; corresponding statistics are undefined.")
  }
  validation_run_signature <- run_signature
  validation_theme <- theme_minimal() + theme(
    panel.grid.minor = element_blank(), strip.text = element_text(size = 11))
  export_table(validation_statistics, "validation_spearman_confirmed", analysis_output_dir)
  export_table(validation_weekly, "validation_weekly_numeric", analysis_output_dir)
  saveRDS(list(weekly = validation_weekly, statistics = validation_statistics,
               run_signature = run_signature), file.path(analysis_output_dir, "validation_data.rds"))
  message("Saved weekly validation data and correlations for four pathogens. Scaling matches totals for a shape comparison.")
} else message("Linked validation skipped.")
