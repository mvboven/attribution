# Shared helpers for the pathogen-attribution analysis.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.

require_objects <- function(objects, envir = parent.frame(), instruction = "Run the earlier stages first.") {
  missing <- objects[!vapply(objects, function(name) {
    exists(name, envir = envir, inherits = FALSE) && !is.null(get(name, envir = envir, inherits = FALSE))
  }, logical(1))]
  if (length(missing)) stop("Missing objects: ", paste(missing, collapse = ", "), ". ", instruction)
}

check_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  invisible(x)
}

check_unique <- function(x, keys, label) {
  if (anyDuplicated(x[, keys, drop = FALSE])) {
    stop(label, " contains duplicate ", paste(keys, collapse = "/"), " combinations.")
  }
  invisible(x)
}

convert_hospital_counts <- function(x, seed = 20240610L) {
  x <- trimws(as.character(x))
  suppressed <- x %in% c("1-4", "1–4")
  if (any(suppressed)) {
    x[suppressed] <- withr::with_seed(seed, sample(1:4, sum(suppressed), replace = TRUE))
  }
  counts <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts)) ||
      any(counts > .Machine$integer.max)) {
    stop("Hospital counts must be finite, non-negative integers or suppressed cells ('1-4').")
  }
  as.integer(counts)
}

input_files_for_scope <- function(scope, directory = "data") {
  stopifnot(length(scope) == 1L, scope %in% c("all_diagnoses", "main_diagnosis", "primary_diagnosis"))
  suffix <- switch(scope, all_diagnoses = "", main_diagnosis = "_main", primary_diagnosis = "_primary")
  list(hospital = file.path(directory, paste0("counts_wk_admissions", suffix, ".csv")),
       linked = file.path(directory, paste0("counts_wk_linked", suffix, ".csv")))
}

# Linked validation uses the midpoint of suppressed 1-4 cells. This is an
# expected count for plotting/comparison, not an integer observation for Stan.
read_confirmed_weekly <- function(file, age_group_lookup, week_seq, age_groups) {
  if (!file.exists(file)) stop("Linked validation file not found: ", file)
  lookup <- tibble::tribble(
    ~source_pathogen, ~component,
    "Influenza_single", "Influenza", "Influenza_coinfection", "Influenza",
    "RSV_single", "RSV", "RSV_coinfection", "RSV",
    "HMPV_single", "HMPV", "HMPV_coinfection", "HMPV",
    "PIV_single", "Parainfluenza", "PIV_coinfection", "Parainfluenza")
  raw <- readr::read_csv(file, col_types = readr::cols(.default = readr::col_character()),
                         show_col_types = FALSE, name_repair = "minimal")
  raw <- raw[, nzchar(names(raw)), drop = FALSE]
  check_columns(raw, c("week_start", "age_group", "pathogen", "n_admissions_v"), "Linked validation data")
  check_unique(raw, c("week_start", "age_group", "pathogen"), "Linked validation data")
  dates <- as.Date(raw$week_start)
  if (anyNA(dates) || any(as.integer(format(dates, "%u")) != 1L)) {
    stop("Linked validation week_start values must be valid ISO Mondays.")
  }
  result <- raw %>% dplyr::transmute(
    week = dates, source_age_group = trimws(age_group), source_pathogen = trimws(pathogen),
    count_character = trimws(n_admissions_v)) %>%
    dplyr::inner_join(age_group_lookup, by = "source_age_group") %>%
    dplyr::inner_join(lookup, by = "source_pathogen") %>%
    dplyr::filter(week %in% week_seq, age_group %in% age_groups) %>%
    dplyr::mutate(confirmed = ifelse(count_character %in% c("1-4", "1–4"), 2.5,
                                    suppressWarnings(as.numeric(count_character))))
  if (any(!is.finite(result$confirmed)) || any(result$confirmed < 0)) {
    stop("Invalid confirmed-pathogen counts in linked validation data.")
  }
  result %>% dplyr::group_by(week, component) %>%
    dplyr::summarise(confirmed = sum(confirmed), .groups = "drop")
}

# Read the complete public validation grid without the original linked extract.
# These counts combine the eight age groups selected when the CSV was exported.
read_compact_confirmed_weekly <- function(file, week_seq, pathogens) {
  if (!file.exists(file)) stop("Compact validation file not found: ", file)
  result <- readr::read_csv(
    file,
    col_types = readr::cols(week = readr::col_date(format = "%Y-%m-%d"),
                           component = readr::col_character(),
                           confirmed = readr::col_double(),
                           .default = readr::col_skip()),
    show_col_types = FALSE, name_repair = "check_unique")
  check_columns(result, c("week", "component", "confirmed"), "Compact validation data")
  if (nrow(readr::problems(result))) stop("Could not parse compact validation values in ", file)
  if (anyNA(result$week) || anyNA(result$component)) stop("Missing compact validation dates/pathogens.")
  check_unique(result, c("week", "component"), "Compact validation data")
  if (any(!is.finite(result$confirmed)) || any(result$confirmed < 0)) {
    stop("Compact confirmed counts must be finite and non-negative.")
  }
  expected <- tidyr::expand_grid(week = as.Date(week_seq), component = pathogens)
  missing <- dplyr::anti_join(expected, result, by = c("week", "component"))
  unexpected <- dplyr::anti_join(result, expected, by = c("week", "component"))
  if (nrow(missing) || nrow(unexpected)) {
    stop("Compact validation data do not match the configured weeks/pathogens: ",
         nrow(missing), " missing and ", nrow(unexpected), " unexpected combinations.")
  }
  result %>% dplyr::arrange(factor(component, levels = pathogens), week)
}

# Construct Mondays from ISO years, retaining leap weeks such as 2015-W53.
make_iso_calendar <- function(first_year, last_year) {
  dates <- seq(as.Date(sprintf("%d-01-01", first_year)) - 7,
               as.Date(sprintf("%d-12-31", last_year)) + 7, by = "day")
  dates <- dates[as.integer(format(dates, "%u")) == 1L]
  year <- as.integer(format(dates, "%G"))
  week <- as.integer(format(dates, "%V"))
  keep <- year >= first_year & year <= last_year
  data.frame(week_start = dates[keep], year = year[keep], week = week[keep],
             ISOweek = sprintf("%04d-W%02d", year[keep], week[keep]))
}

# Winter and other_weeks partition the ENTIRE study period, including the
# winter fragments at its boundaries. This differs from complete-winters-only.
make_periods <- function(week, winter_start = 44L, winter_end = 13L) {
  stopifnot(length(week) > 0L, all(is.finite(week)), all(week %in% 1:53),
            winter_start %in% 1:53, winter_end %in% 1:53, winter_start > winter_end)
  winter <- which(week >= winter_start | week <= winter_end)
  other <- setdiff(seq_along(week), winter)
  if (!length(winter) || !length(other)) stop("Winter and other-weeks periods must both contain weeks.")
  list(overall = seq_along(week), winter = winter, other_weeks = other)
}

object_md5 <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(x, path, version = 2)
  unname(tools::md5sum(path))
}

# Select explicit names before reshaping: never assume CmdStan's output order.
indexed_point_estimate <- function(values, variable, dimensions) {
  grid <- do.call(expand.grid, lapply(dimensions, seq_len))
  expected <- apply(grid, 1, function(i) paste0(variable, "[", paste(i, collapse = ","), "]"))
  if (is.null(names(values)) || !all(expected %in% names(values))) {
    stop("MAP output lacks the expected indices for ", variable, ". Check the Stan model.")
  }
  result <- as.numeric(values[expected])
  if (any(!is.finite(result))) stop("Non-finite MAP estimate for ", variable)
  if (length(dimensions) == 1L) result else array(result, dim = dimensions)
}

get_draw_matrix <- function(fit, variables) {
  draws <- as.matrix(fit$draws(variables = variables, format = "draws_matrix"))
  if (!all(variables %in% colnames(draws))) stop("Expected Stan variables are missing from the fit.")
  draws <- draws[, variables, drop = FALSE]
  if (!nrow(draws) || any(!is.finite(draws))) stop("Posterior draws are empty or non-finite.")
  draws
}

# Conditional expected allocations, not a simulation of latent multinomial counts.
# Fractions are based on sums of allocated counts / sums of observed counts.
allocate_observed <- function(contribution, intensity, observed) {
  stopifnot(is.matrix(contribution), is.matrix(intensity),
            identical(dim(contribution), dim(intensity)), ncol(intensity) == length(observed))
  if (any(!is.finite(intensity)) || any(intensity <= 0) ||
      any(!is.finite(contribution)) || any(!is.finite(observed)) || any(observed < 0)) {
    stop("Invalid intensities, contributions, or observations in attribution allocation.")
  }
  probability <- contribution / intensity
  if (any(probability < -1e-10) || any(probability > 1 + 1e-10)) {
    stop("Pathogen allocation probabilities are outside [0, 1].")
  }
  sweep(probability, 2, observed, "*")
}

interval <- function(x) {
  if (!length(x) || any(!is.finite(x))) stop("Cannot summarise empty or non-finite values.")
  c(median = stats::median(x), lower = unname(stats::quantile(x, 0.025)),
    upper = unname(stats::quantile(x, 0.975)))
}

format_interval <- function(median, lower, upper, digits = 1L, percent = FALSE) {
  pattern <- paste0("%.", digits, "f", if (percent) "%%" else "",
                    " (%.", digits, "f–%.", digits, "f)")
  sprintf(pattern, median, lower, upper)
}

# UTF-8 BOM makes intervals readable when opening the CSV directly in Excel.
# The tab-delimited companion can be copied into Word or converted to a table.
export_table <- function(x, stem, directory, word = FALSE) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  csv <- file.path(directory, paste0(stem, ".csv"))
  if (word) {
    readr::write_excel_csv(x, csv, na = "")
    readr::write_tsv(x, file.path(directory, paste0(stem, ".tsv")), na = "")
  } else {
    readr::write_csv(x, csv, na = "")
  }
  invisible(csv)
}

compile_model <- function(model_builder = cmdstanr::cmdstan_model, cmdstan_root = NULL) {
  if (!file.exists(stan_file)) stop("Stan model not found: ", stan_file, ". Add the original model file.")
  build_root <- get0("stan_build_dir", ifnotfound = "", inherits = TRUE)
  if (!nzchar(build_root)) {
    if (is.null(cmdstan_root)) cmdstan_root <- cmdstanr::cmdstan_path()
    build_root <- file.path(cmdstan_root, "attribution_models")
  }
  # Key the cache by source content so separate projects cannot accidentally
  # use a different model with the same filename.
  model_hash <- unname(tools::md5sum(stan_file))
  executable_dir <- file.path(path.expand(build_root), paste0(
    tools::file_path_sans_ext(basename(stan_file)), "_", model_hash))
  dir.create(executable_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(executable_dir) || file.access(executable_dir, 2L) != 0L) {
    stop("Cannot write to Stan build directory: ", executable_dir,
         ". Set stan_build_dir in scripts/0_source.R to an allowed directory.")
  }
  executable_dir <- normalizePath(executable_dir, winslash = "/", mustWork = TRUE)
  model_builder(stan_file, dir = executable_dir)
}

# Recover a saved fit without repeating MAP or HMC. The data/model signature
# prevents accidentally analysing chains fitted to a different data set or prior.
require_fit <- function(envir = parent.frame()) {
  require_objects(c("data_list", "run_signature"), envir,
                  "Run run_analysis(c(1, 2)) after sourcing scripts/0_source.R first.")
  current <- get("run_signature", envir)
  fit <- get0("fit_hmc", envir = envir, inherits = FALSE)
  signature <- get0("fit_run_signature", envir = envir, inherits = FALSE)
  if (!is.null(fit)) {
    if (!identical(signature, current)) stop("Loaded HMC fit does not match the current data/model settings.")
    return(invisible(fit))
  }
  if (!file.exists(hmc_manifest_file)) stop("No saved HMC run. Run run_analysis(c(3, 4)) first.")
  manifest <- readRDS(hmc_manifest_file)
  if (!identical(manifest$run_signature, current)) {
    stop("Saved HMC fit does not match the current data/model settings. Refit or restore those settings.")
  }
  paths <- file.path(project_dir, manifest$csv_files)
  if (!all(file.exists(paths))) stop("Some HMC CSV files are missing. Keep the output/hmc_csv files together.")
  fit <- cmdstanr::as_cmdstan_fit(paths, check_diagnostics = FALSE)
  assign("fit_hmc", fit, envir)
  assign("fit_run_signature", current, envir)
  message("Loaded the saved HMC chains for ", model_id, ".")
  invisible(fit)
}

# Return posterior summaries without depending on any figure script.
summarise_pathogen_intensities <- function(fit, week_seq, age_groups, pathogen_names) {
  result <- fit$summary(variables = "lambda_contributions", "median",
                       quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975)))
  result <- tidyr::extract(result, variable,
                          into = c("pathogen_index", "age_index", "week_index"),
                          regex = "^lambda_contributions\\[(\\d+),(\\d+),(\\d+)\\]$", convert = TRUE)
  if (anyNA(result[, c("pathogen_index", "age_index", "week_index")]) ||
      any(!result$pathogen_index %in% seq_along(pathogen_names)) ||
      any(!result$age_index %in% seq_along(age_groups)) ||
      any(!result$week_index %in% seq_along(week_seq))) stop("Unexpected lambda_contributions indices.")
  result <- dplyr::transmute(result, week = week_seq[week_index], y_data = median,
                            y_lo = q2.5, y_hi = q97.5, age_group = age_groups[age_index],
                            component = pathogen_names[pathogen_index])
  check_unique(result, c("week", "age_group", "component"), "Pathogen posterior summary")
  stopifnot(nrow(result) == length(week_seq) * length(age_groups) * length(pathogen_names),
            all(is.finite(result$y_data)), all(result$y_data >= 0),
            all(is.finite(result$y_lo)), all(is.finite(result$y_hi)))
  result
}

# Preserve the original validation estimand: SUM OF AGE-SPECIFIC MEDIANS of
# fitted pathogen intensities, rather than the conditional observed allocation.
# This is a timing/shape comparison, not an independent calibration of burden.
block_spearman <- function(x, y, n_boot = 2000L, block_length = 13L, seed = 20260903L) {
  stopifnot(length(x) == length(y), length(x) >= 2L, all(is.finite(x)), all(is.finite(y)),
            n_boot >= 2L, block_length >= 1L, block_length <= length(x))
  rho <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  if (!is.finite(rho)) {
    return(data.frame(spearman_rho = NA_real_, rho_ci_low = NA_real_, rho_ci_high = NA_real_))
  }
  rho_boot <- withr::with_seed(seed, replicate(n_boot, {
    n <- length(x)
    n_blocks <- ceiling(n / block_length)
    starts <- sample.int(n, n_blocks, replace = TRUE)
    indices <- (rep(starts, each = block_length) - 1L +
                  rep(seq_len(block_length) - 1L, n_blocks)) %% n + 1L
    indices <- indices[seq_len(n)]
    suppressWarnings(stats::cor(x[indices], y[indices], method = "spearman"))
  }))
  valid <- rho_boot[is.finite(rho_boot)]
  if (length(valid) < 0.9 * n_boot) warning("More than 10% of bootstrap correlations were undefined.")
  ci <- if (length(valid)) unname(stats::quantile(valid, c(0.025, 0.975))) else c(NA_real_, NA_real_)
  data.frame(spearman_rho = rho, rho_ci_low = ci[1], rho_ci_high = ci[2])
}

round_log_down <- function(x) {
  stopifnot(is.finite(x), x > 0)
  exponent <- floor(log10(x))
  values <- c(1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10)
  max(values[values <= x / 10^exponent]) * 10^exponent
}
round_log_up <- function(x) {
  stopifnot(is.finite(x), x > 0)
  exponent <- floor(log10(x))
  values <- c(1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10)
  min(values[values >= x / 10^exponent]) * 10^exponent
}
format_multiplier <- function(x) {
  vapply(x, function(z) format(round(z, if (z >= 10) 1 else if (z >= 1) 2 else 3),
                              trim = TRUE, scientific = FALSE), character(1))
}
format_beta_axis <- function(x) {
  vapply(x, function(z) {
    if (!is.finite(z)) NA_character_ else if (z < 0.01) {
      format(z, scientific = TRUE, digits = 1, trim = TRUE)
    } else format(round(z, 2), scientific = FALSE, trim = TRUE)
  }, character(1))
}
