# Compare RSV table percentages with fitted and plotted RSV shares.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
#
# In the configured analysis session, after stages 1 and 2, run:
#   source("check_rsv_attribution.R")
# require_fit() reloads existing saved chains when necessary; it never refits.
# The comparison is stored as rsv_attribution_check and printed for children.
# The supplied table and plotting scripts are not modified.

rsv_attribution_check <- local({
  require_fit(envir = .GlobalEnv)
  require_objects(c("data", "y_data", "n_data", "age_groups",
                    "pathogen_lookup", "winter_start_week", "winter_end_week"),
                  envir = .GlobalEnv)
  
  rsv_index <- which(as.character(pathogen_lookup$pathogen) == "RSV")
  stopifnot(length(rsv_index) == 1L,
            nrow(y_data) == n_data,
            ncol(y_data) == length(age_groups),
            identical(colnames(y_data), as.character(age_groups)))
  periods <- make_periods(data$week, winter_start_week, winter_end_week)
  child_indices <- match(c("0-2m", "3-5m", "6-11m", "1-2y", "3-4y"),
                         as.character(age_groups))
  if (anyNA(child_indices)) stop("Expected child age groups are missing.")
  
  # Request and then explicitly reorder every age/week column by its name.
  lambda_names <- unlist(lapply(child_indices, function(a) {
    sprintf("lambda[%d,%d]", a, seq_len(n_data))
  }), use.names = FALSE)
  rsv_names <- unlist(lapply(child_indices, function(a) {
    sprintf("lambda_contributions[%d,%d,%d]", rsv_index, a, seq_len(n_data))
  }), use.names = FALSE)
  lambda <- get_draw_matrix(fit_hmc, lambda_names)
  rsv <- get_draw_matrix(fit_hmc, rsv_names)
  stopifnot(identical(dim(lambda), dim(rsv)),
            identical(colnames(lambda), lambda_names),
            identical(colnames(rsv), rsv_names),
            all(lambda > 0), all(rsv >= 0),
            all(rsv <= lambda * (1 + 1e-8)))
  
  # Pointwise posterior medians reproduce the lines used by the figure.
  lambda_weekly_median <- apply(lambda, 2, stats::median)
  rsv_weekly_median <- apply(rsv, 2, stats::median)
  summarise_draws <- function(x) {
    c(median = stats::median(x),
      lower = unname(stats::quantile(x, 0.025)),
      upper = unname(stats::quantile(x, 0.975)))
  }
  
  rows <- lapply(seq_along(child_indices), function(j) {
    a <- child_indices[j]
    columns <- (j - 1L) * n_data + seq_len(n_data)
    fitted <- lambda[, columns, drop = FALSE]
    contribution <- rsv[, columns, drop = FALSE]
    observed <- y_data[, a]
    
    # Recompute independently of allocate_observed(), using the table's formula.
    allocated <- sweep(contribution / fitted, 2, observed, "*")
    lapply(names(periods), function(period) {
      weeks <- periods[[period]]
      observed_total <- sum(observed[weeks])
      if (observed_total <= 0) stop("No observed admissions in a selected period.")
      rsv_total <- rowSums(contribution[, weeks, drop = FALSE])
      fitted_total <- rowSums(fitted[, weeks, drop = FALSE])
      allocated_total <- rowSums(allocated[, weeks, drop = FALSE])
      allocated_pct <- summarise_draws(100 * allocated_total / observed_total)
      fitted_pct <- summarise_draws(100 * rsv_total / fitted_total)
      
      data.frame(
        age_group = as.character(age_groups[a]),
        period = period,
        n_weeks = length(weeks),
        observed_total = observed_total,
        allocated_RSV_total_median = stats::median(allocated_total),
        fitted_RSV_total_median = stats::median(rsv_total),
        fitted_total_median = stats::median(fitted_total),
        allocated_pct_median = unname(allocated_pct["median"]),
        allocated_pct_lower = unname(allocated_pct["lower"]),
        allocated_pct_upper = unname(allocated_pct["upper"]),
        fitted_pct_median = unname(fitted_pct["median"]),
        fitted_pct_lower = unname(fitted_pct["lower"]),
        fitted_pct_upper = unname(fitted_pct["upper"]),
        curve_area_pct = 100 * sum(rsv_weekly_median[columns[weeks]]) /
          sum(lambda_weekly_median[columns[weeks]]),
        fitted_total_vs_observed_pct = stats::median(
          100 * fitted_total / observed_total),
        stringsAsFactors = FALSE
      )
    })
  })
  result <- do.call(rbind, unlist(rows, recursive = FALSE))
  rownames(result) <- NULL
  
  # When stage 9 has already run, compare with its actual numeric table object.
  table_summary <- get0("attribution_summary", envir = .GlobalEnv,
                        inherits = FALSE)
  required <- c("age_group", "period", "pathogen", "fraction_observed_median")
  if (!is.null(table_summary) && all(required %in% names(table_summary))) {
    existing <- as.data.frame(table_summary)
    existing <- existing[as.character(existing$pathogen) == "RSV", , drop = FALSE]
    key <- paste(as.character(existing$age_group), as.character(existing$period),
                 sep = "|")
    if (anyDuplicated(key)) stop("Duplicate RSV rows in attribution_summary.")
    index <- match(paste(result$age_group, result$period, sep = "|"), key)
    result$table_pct_median <- existing$fraction_observed_median[index]
    result$table_minus_recomputed_pp <- result$table_pct_median -
      result$allocated_pct_median
    if (anyNA(index)) {
      warning("Some period/age rows are absent from the existing RSV table.")
    }
    discrepancy <- abs(result$table_minus_recomputed_pp)
    if (any(discrepancy > 1e-7, na.rm = TRUE)) {
      warning("Existing RSV table differs from the directly recomputed allocation.")
    } else if (!anyNA(index)) {
      message("Existing RSV table matches the directly recomputed allocation.")
    }
  }
  
  cat("\nRSV percentages by child age group and period\n",
      "allocated_pct: median of 100 * sum(y * RSV/lambda) / sum(y) [table]\n",
      "fitted_pct: median of 100 * sum(RSV) / sum(lambda)\n",
      "curve_area_pct: 100 * sum(weekly RSV medians) / sum(weekly lambda medians)\n",
      "fitted_total_vs_observed_pct: fitted/observed total; 100 means agreement\n",
      "Winter uses the configured weeks; overall includes all study weeks.\n\n",
      sep = "")
  print_columns <- c("age_group", "period", "allocated_pct_median",
                     "fitted_pct_median", "curve_area_pct",
                     "fitted_total_vs_observed_pct")
  print_columns <- c(print_columns,
                     intersect(c("table_pct_median", "table_minus_recomputed_pp"),
                               names(result)))
  display <- result[, print_columns, drop = FALSE]
  numeric_columns <- vapply(display, is.numeric, logical(1))
  display[numeric_columns] <- lapply(display[numeric_columns], round, digits = 2)
  print(display, row.names = FALSE)
  result
})
