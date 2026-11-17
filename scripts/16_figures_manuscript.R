# Standalone manuscript figures from saved CmdStan output.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
# Save as scripts/16_figures_manuscript.R and source from the project root.
# Loads shared helper functions; does not compile Stan or run MAP/HMC/LOO.
# Requires the saved manifests, chain CSVs, prepared_stan_data.rds, and compact validation CSVs.
# Current Stan source is NOT used. Saved data are checked against each fit.

(function() {
  # -------------------------- USER SETTINGS -------------------------------
  project_root <- getwd()
  num_knots <- 200L
  output_folder <- file.path(project_root, "figures", "manuscript",
                             sprintf("knots_%03d", num_knots))
  show_figures <- FALSE
  png_dpi <- 600
  pathogen_colours <- c(Influenza = "#0072B2", RSV = "#D55E00",
                        HMPV = "#009E73", Parainfluenza = "#CC79A7")
  # Set these explicitly if your saved fits are elsewhere.
  run_dirs <- c(
    all_diagnoses = file.path(project_root, "output", "all_diagnoses",
                              sprintf("hmc_all_diagnoses_background_knots_%02d", num_knots)),
    main_diagnosis = file.path(project_root, "output", "main_diagnosis",
                               sprintf("hmc_main_diagnosis_background_knots_%02d", num_knots))
  )
  validation_files <- c(
    all_diagnoses = file.path(project_root, "data", "compact", "validation_confirmed_all_diagnoses.csv"),
    main_diagnosis = file.path(project_root, "data", "compact", "validation_confirmed_main_diagnosis.csv")
  )
  # ------------------------------------------------------------------------
  packages <- c("cmdstanr", "posterior", "dplyr", "tidyr", "tibble", "readr", "ggplot2")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  invisible(lapply(packages, function(p) suppressPackageStartupMessages(
    library(p, character.only = TRUE))))
  source(file.path(project_root, "scripts", "analysis_helpers.R"), local = TRUE)
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  output_folder <- normalizePath(output_folder, winslash = "/", mustWork = TRUE)
  message("Manuscript output directory: ", output_folder)
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  
  pathogen_order <- names(pathogen_colours)
  panel_labels <- setNames(c("A  Influenza", "B  RSV", "C  HMPV", "D  PIV"), pathogen_order)
  age_labels <- c("0-2m" = "0-2 months", "3-5m" = "3-5 months",
                  "6-11m" = "6-11 months", "1-2y" = "1-2 years", "3-4y" = "3-4 years",
                  "65-74y" = "65-74 years", "75-84y" = "75-84 years", "85+y" = "85+ years")
  age_alpha <- setNames(c(.25, .35, .45, .55, .65, .75, .87, 1), names(age_labels))
  
  hash_object <- function(x) {
    f <- tempfile(fileext = ".rds")
    on.exit(unlink(f))
    saveRDS(x, f, version = 2)
    unname(tools::md5sum(f))
  }
  theme_ms <- function() {
    theme_minimal(base_size = 9, base_family = "sans") +
      theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(colour = "grey92", linewidth = .2),
            axis.line = element_line(colour = "grey40", linewidth = .25),
            axis.ticks = element_line(colour = "grey40", linewidth = .25),
            axis.ticks.length = grid::unit(1, "mm"),
            axis.text = element_text(size = 8, colour = "grey20"),
            axis.title = element_text(size = 9, colour = "grey20"),
            axis.title.x = element_text(margin = margin(t = 6)),
            axis.title.y = element_text(margin = margin(r = 6)),
            strip.background = element_blank(),
            strip.text = element_text(size = 10, face = "bold", hjust = 0,
                                      margin = margin(t = 3, b = 6)),
            panel.spacing = grid::unit(5, "mm"),
            plot.margin = margin(5, 4, 3, 4),
            legend.position = "bottom", legend.title = element_blank(),
            legend.text = element_text(size = 8),
            legend.key.width = grid::unit(7, "mm"))
  }
  figure_index <- list()
  save_figure <- function(plot, stem, height, scope, signature) {
    for (extension in c("pdf", "png")) {
      args <- list(filename = file.path(output_folder, paste0(stem, ".", extension)),
                   plot = plot, width = 180, height = height, units = "mm", bg = "white")
      if (extension == "pdf") args$device <- pdf_device else args$dpi <- png_dpi
      do.call(ggsave, args)
    }
    if (show_figures) print(plot)
    figure_index[[length(figure_index) + 1L]] <<- data.frame(
      figure = stem, diagnosis_scope = scope, run_signature = signature,
      width_mm = 180, height_mm = height)
    write.csv(bind_rows(figure_index), file.path(output_folder, "manuscript_figure_index.csv"),
              row.names = FALSE)
    message("Saved ", stem, " (PDF and PNG).")
  }
  
  # Load only the named variables, and summarise before reshaping.
  summarise_variable <- function(fit, variable, indices) {
    draws <- as.matrix(fit$draws(variables = variable, format = "draws_matrix"))
    if (!nrow(draws) || any(!is.finite(draws))) stop("Invalid draws for ", variable)
    q <- apply(draws, 2, quantile, probs = c(.025, .5, .975), names = FALSE)
    result <- tibble(variable = colnames(draws), mean = colMeans(draws),
                     lower = q[1, ], median = q[2, ], upper = q[3, ])
    pattern <- paste0("^", variable, "\\[",
                      paste(rep("(\\d+)", length(indices)), collapse = ","), "\\]$")
    result <- tidyr::extract(result, variable, into = indices, regex = pattern, convert = TRUE)
    if (anyNA(result[, indices, drop = FALSE])) stop("Unexpected indices for ", variable)
    result
  }
  read_linked <- function(file, weeks, ages) {
    if (!setequal(ages, names(age_labels))) {
      stop("The compact validation counts were exported for the eight study age groups.")
    }
    read_compact_confirmed_weekly(file, weeks, pathogen_order) %>%
      rename(pathogen = component)
  }

  # Verify saved plotting inputs using the SOURCE HASH RECORDED AT FIT TIME.
  # This allows current source-code changes without attaching unrelated data.
  bundles <- list()
  for (scope in names(run_dirs)) {
    folder <- run_dirs[[scope]]
    paths <- c(file.path(folder, "hmc_manifest.rds"),
               file.path(folder, "prepared_stan_data.rds"), validation_files[[scope]])
    if (!all(file.exists(paths))) stop("Missing required saved input:\n",
                                       paste(paths[!file.exists(paths)], collapse = "\n"))
    m <- readRDS(paths[1]); b <- readRDS(paths[2])
    required <- c("data", "data_list", "age_groups", "pathogen_lookup")
    if (!all(required %in% names(b))) stop("Incomplete prepared data: ", paths[2])
    if (is.null(m$stan_md5) || is.null(m$run_signature) || !length(m$csv_files)) {
      stop("Incomplete saved fit manifest: ", paths[1])
    }
    signature <- hash_object(list(b$data_list, m$stan_md5, b$age_groups, b$pathogen_lookup))
    if (!identical(signature, m$run_signature)) {
      stop("Saved plotting data do not match the fitted inputs for ", scope, ".\n",
           "Restore prepared_stan_data.rds from that fit (not from a later preparation run).\n",
           "A current Stan-file change alone does not cause this check to fail.")
    }
    weeks <- as.Date(b$data$week_start)
    observations <- t(b$data_list$y_data)
    observation_hash <- hash_object(list(scope, weeks, b$age_groups, observations))
    if (!is.null(m$observation_signature) &&
        !identical(observation_hash, m$observation_signature)) {
      stop("Observed admissions do not match the saved fit for ", scope)
    }
    csv <- vapply(m$csv_files, function(p) {
      if (grepl("^(/|[A-Za-z]:[/\\\\])", p)) return(p)
      file.path(project_root, p)
    }, character(1))
    if (!all(file.exists(csv))) stop("Missing chain CSVs for ", scope, ":\n",
                                     paste(csv[!file.exists(csv)], collapse = "\n"))
    if (!setequal(as.character(b$pathogen_lookup$pathogen), pathogen_order) ||
        !all(b$age_groups %in% names(age_labels)) || anyNA(weeks) ||
        nrow(observations) != length(weeks) || any(diff(weeks) != 7)) {
      stop("Unexpected saved data layout for ", scope)
    }
    bundles[[scope]] <- list(manifest = m, prepared = b, csv = csv,
                             weeks = weeks, observations = observations)
    message("Verified saved plotting inputs for ", scope, ".")
  }
  if (!identical(bundles[[1]]$manifest$stan_md5, bundles[[2]]$manifest$stan_md5)) {
    warning("The two fits have different saved Stan source hashes. They are exported as fitted; ",
            "check that their model/prior differences are intentional.")
  }
  
  all_scaling <- list()
  provenance <- list()
  for (scope in names(bundles)) {
    bundle <- bundles[[scope]]; b <- bundle$prepared; m <- bundle$manifest
    weeks <- bundle$weeks; ages <- b$age_groups
    pathogen_names <- as.character(b$pathogen_lookup$pathogen)
    fit <- cmdstanr::as_cmdstan_fit(bundle$csv, check_diagnostics = FALSE)
    message("Reading posterior summaries for ", scope, ".")
    suffix <- if (scope == "all_diagnoses") "all" else "main"
    main_number <- if (suffix == "all") "1" else "2"
    seasonal_number <- if (suffix == "all") "S1" else "S2"
    delay_number <- if (suffix == "all") "S3" else "S4"
    validation_number <- if (suffix == "all") "3" else "S5"
    date_scale <- function() scale_x_date(
      breaks = as.Date(paste0(seq(2013, 2019, 2), "-01-01")), date_labels = "%Y",
      limits = range(weeks), expand = expansion(mult = c(.01, .01)))
    
    total <- summarise_variable(fit, "lambda", c("a", "w")) %>%
      mutate(week = weeks[w], age = factor(ages[a], levels = ages))
    contributions <- summarise_variable(fit, "lambda_contributions", c("p", "a", "w")) %>%
      mutate(week = weeks[w], age = factor(ages[a], levels = ages),
             pathogen = factor(pathogen_names[p], levels = pathogen_order))
    if (anyNA(total$age) || anyNA(contributions$pathogen) || anyNA(contributions$age) ||
        anyNA(total$week) || anyNA(contributions$week)) stop("Invalid posterior indices.")
    observed <- expand_grid(a = seq_along(ages), w = seq_along(weeks)) %>%
      mutate(week = weeks[w], age = factor(ages[a], levels = ages),
             count = bundle$observations[cbind(w, a)])
    p <- ggplot() +
      geom_ribbon(data = contributions, aes(week, ymin = lower, ymax = upper,
                                            fill = pathogen, group = pathogen), alpha = .15, colour = NA) +
      geom_line(data = total, aes(week, median), colour = "grey65", linewidth = .5) +
      geom_line(data = contributions, aes(week, median, colour = pathogen), linewidth = .4) +
      geom_point(data = observed, aes(week, count), colour = "black", size = .3, alpha = .8) +
      facet_wrap(~age, ncol = 3, scales = "free_y", labeller = as_labeller(age_labels)) +
      scale_colour_manual(values = pathogen_colours, breaks = pathogen_order,
                          labels = c("Influenza", "RSV", "HMPV", "PIV")) +
      scale_fill_manual(values = pathogen_colours, guide = "none") + date_scale() +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .06))) +
      labs(x = "Year", y = "Weekly RTI hospital admissions") + theme_ms()
    save_figure(p, paste0("Figure_", main_number, "_attribution_", suffix), 160, scope, m$run_signature)
    
    season_labels <- paste0(sort(unique(b$data$season_start_year)), "/",
                            sort(unique(b$data$season_start_year)) + 1L)
    if (!length(season_labels)) stop("Saved data lack season_start_year.")
    sd <- as.data.frame(fit$draws("seasoneffect", format = "draws_df")) %>%
      select(starts_with("seasoneffect[")) %>%
      pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
      tidyr::extract(variable, c("p", "s"), "^seasoneffect\\[(\\d+),(\\d+)\\]$", convert = TRUE) %>%
      mutate(pathogen = factor(pathogen_names[p], levels = pathogen_order),
             season = factor(season_labels[s], levels = season_labels))
    if (anyNA(sd) || any(sd$value <= 0)) stop("Invalid seasonal multipliers.")
    sm <- sd %>% group_by(pathogen, season) %>% summarise(median = median(value), .groups = "drop")
    p <- ggplot(sd, aes(season, value)) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = .3) +
      geom_violin(aes(fill = pathogen, colour = pathogen), alpha = .25,
                  scale = "width", trim = TRUE, linewidth = .35) +
      geom_point(data = sm, aes(y = median, colour = pathogen), size = 1.4) +
      facet_wrap(~pathogen, ncol = 2, scales = "free_y", labeller = as_labeller(panel_labels)) +
      scale_colour_manual(values = pathogen_colours, guide = "none") +
      scale_fill_manual(values = pathogen_colours, guide = "none") +
      scale_y_log10(labels = function(x) format(x, trim = TRUE, scientific = FALSE)) +
      labs(x = "Season", y = "Seasonal multiplier") + theme_ms() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_figure(p, paste0("Figure_", seasonal_number, "_seasonal_", suffix), 145, scope, m$run_signature)
    
    lag_values <- seq.int(b$data_list$min_lag, b$data_list$max_lag)
    ds <- summarise_variable(fit, "pmf", c("a", "p", "l")) %>%
      mutate(age = factor(ages[a], levels = ages),
             pathogen = factor(pathogen_names[p], levels = pathogen_order),
             lag = factor(lag_values[l], levels = lag_values))
    sums <- ds %>% group_by(age, pathogen) %>% summarise(total = sum(mean), .groups = "drop")
    if (anyNA(ds) || any(abs(sums$total - 1) > 1e-6)) stop("Invalid lag probabilities.")
    dodge <- position_dodge(width = .85)
    p <- ggplot(ds, aes(lag, mean, group = age)) +
      geom_col(aes(fill = pathogen, alpha = age), position = dodge, width = .78) +
      geom_errorbar(aes(ymin = lower, ymax = upper, colour = pathogen),
                    position = dodge, width = .10, linewidth = .25, show.legend = FALSE) +
      facet_wrap(~pathogen, ncol = 2, labeller = as_labeller(panel_labels)) +
      scale_fill_manual(values = pathogen_colours, guide = "none") +
      scale_colour_manual(values = pathogen_colours, guide = "none") +
      scale_alpha_manual(values = age_alpha, breaks = ages, name = "Age group",
                         guide = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(fill = "grey20"))) +
      scale_y_continuous(limits = c(0, 1 + 1e-6), breaks = seq(0, 1, .25),
                         expand = expansion(mult = c(0, .03))) +
      labs(x = "Lag (weeks)", y = "Probability") + theme_ms() +
      theme(legend.title = element_text(size = 9))
    save_figure(p, paste0("Figure_", delay_number, "_delay_", suffix), 150, scope, m$run_signature)
    
    # Preserve the validation estimand: sum of age-specific posterior medians
    # of fitted intensities, not observed-count allocation and not median of sum.
    attributed <- contributions %>% group_by(week, pathogen) %>%
      summarise(attributed = sum(median), .groups = "drop") %>%
      mutate(pathogen = as.character(pathogen))
    confirmed <- read_linked(validation_files[[scope]], weeks, ages)
    vd <- expand_grid(week = weeks, pathogen = pathogen_order) %>%
      left_join(confirmed, by = c("week", "pathogen")) %>%
      mutate(confirmed = replace_na(confirmed, 0)) %>%
      left_join(attributed, by = c("week", "pathogen"))
    scaling <- vd %>% group_by(pathogen) %>%
      summarise(confirmed_total = sum(confirmed), attributed_total = sum(attributed),
                .groups = "drop") %>%
      mutate(confirmed_to_attributed = attributed_total / confirmed_total,
             attributed_to_confirmed = confirmed_total / attributed_total,
             diagnosis_scope = scope,
             pathogen = factor(pathogen, levels = pathogen_order)) %>% arrange(pathogen)
    if (any(!is.finite(scaling$confirmed_to_attributed)) ||
        any(scaling$confirmed_to_attributed <= 0)) stop("Invalid validation scaling.")
    vd <- vd %>% left_join(scaling %>% mutate(pathogen = as.character(pathogen)) %>%
                             select(pathogen, attributed_to_confirmed), by = "pathogen") %>%
      mutate(scaled_attributed = attributed * attributed_to_confirmed,
             pathogen = factor(pathogen, levels = pathogen_order))
    p <- ggplot(vd, aes(week)) +
      geom_line(aes(y = confirmed, linetype = "Confirmed"),
                colour = "grey25", linewidth = .3, alpha = .85) +
      geom_line(aes(y = scaled_attributed, colour = pathogen,
                    linetype = "Scaled fitted attribution"), linewidth = .55) +
      facet_wrap(~pathogen, ncol = 2, scales = "free_y", labeller = as_labeller(panel_labels)) +
      scale_colour_manual(values = pathogen_colours, guide = "none") +
      scale_linetype_manual(values = c("Confirmed" = "solid", "Scaled fitted attribution" = "solid"),
                            breaks = c("Confirmed", "Scaled fitted attribution"),
                            guide = guide_legend(nrow = 1, override.aes = list(
                              colour = c("grey25", "grey40"), linewidth = c(.3, .55), alpha = c(.85, 1)))) +
      date_scale() + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .06))) +
      labs(x = "Year", y = "Pathogen-linked admissions per week") + theme_ms()
    save_figure(p, paste0("Figure_", validation_number, "_validation_", suffix), 140, scope, m$run_signature)
    all_scaling[[scope]] <- scaling
    write.csv(bind_rows(all_scaling), file.path(output_folder, "validation_scaling_factors.csv"), row.names = FALSE)
    cat("\nConfirmed-to-attributed multipliers: ", scope, "\n", sep = "")
    print(as.data.frame(scaling[, c("pathogen", "confirmed_to_attributed")]), row.names = FALSE)
    write.csv(vd, file.path(output_folder, paste0("validation_weekly_", suffix, ".csv")), row.names = FALSE)
    provenance[[scope]] <- list(manifest = m, chain_files = bundle$csv,
                                validation_file = validation_files[[scope]], validation_md5 = unname(tools::md5sum(validation_files[[scope]])))
    saveRDS(provenance, file.path(output_folder, "manuscript_figure_provenance.rds"))
    rm(fit, sd, contributions, total, ds, p)
    invisible(gc())
  }
  message("\nCompleted all eight manuscript figures in:\n", output_folder)
  message("Current Stan source and scripts/0_source.R were not used; no models were refitted.")
})()
