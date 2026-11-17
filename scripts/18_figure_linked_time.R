# Monthly virological surveillance series for four respiratory pathogens.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
# Standalone replacement for the earlier linked-admissions script 18.
# Source from the project root or scripts/; no Stan fit or script 17 required.
# These are laboratory detections, without age or hospital-diagnosis strata.

(function() {
  packages <- c("readxl", "dplyr", "tidyr", "ggplot2")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install these packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
  
  # Change these settings through options(), if necessary.
  root <- getOption("attribution.project_root", getwd())
  if (basename(root) == "scripts" && dir.exists(file.path(root, "..", "data"))) root <- file.path(root, "..")
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  input_file <- getOption("attribution.piv_virology_file",
                          file.path(root, "data", "virologische_weekstaten_20240322.xlsx"))
  input_sheet <- getOption("attribution.piv_virology_sheet", "VW_2013-2022")
  figure_dir <- getOption("attribution.supplement_dir", file.path(root, "figures", "manuscript", "supplement"))
  data_dir <- getOption("attribution.virology_data_dir", file.path(root, "data"))
  show_figures <- isTRUE(getOption("attribution.supplement_show", FALSE))
  # Use the same workbook and sheet options as script 17.
  types <- c("Influenza", "RSV", "HMPV", "Parainfluenza")
  influenza_columns <- c("VIRINFA", "VIRINFB")
  piv_columns <- c("VIRPAR1", "VIRPAR2", "VIRPAR3", "VIRPAR4",
                   "VIRPARANGETYP", "VIRPARANTYP")
  source_columns <- c(influenza_columns, "VIRRS", "VIRHMPV", piv_columns)
  if (!file.exists(input_file)) stop("Virological workbook not found: ", input_file)
  
  # All PIV categories are pooled, including untyped/non-typeable detections.
  raw <- readxl::read_excel(input_file, sheet = input_sheet)
  required <- c("Jaar", "Week", "Aantal_laboratoria", source_columns)
  if (!all(required %in% names(raw))) stop("Missing columns: ", paste(setdiff(required, names(raw)), collapse = ", "))
  for (name in required) {
    z <- raw[[name]]
    value <- suppressWarnings(as.numeric(z))
    if (any(!is.na(z) & is.na(value)) || any(!is.na(value) &
                                             (!is.finite(value) | value < 0 | value != floor(value)))) stop("Invalid counts in ", name)
    raw[[name]] <- value
  }
  if (anyNA(raw$Jaar) || anyNA(raw$Week)) stop("Missing year or week.")
  
  # ISO week 1 contains 4 January. Reject non-existent week-53 entries.
  january4 <- as.Date(sprintf("%04d-01-04", as.integer(raw$Jaar)))
  raw$week_start <- january4 - (as.integer(format(january4, "%u")) - 1L) + 7L * (raw$Week - 1L)
  valid <- raw$Week >= 1 & raw$Week <= 53 & !is.na(raw$week_start) &
    as.integer(format(raw$week_start + 3, "%Y")) == raw$Jaar
  if (anyNA(valid)) stop("Could not construct ISO dates.")
  omitted <- raw[!valid, required, drop = FALSE]
  if (nrow(omitted)) {message("Omitting invalid ISO weeks:"); print(as.data.frame(omitted), row.names = FALSE)}
  raw <- raw[valid, , drop = FALSE]
  if (!nrow(raw) || anyDuplicated(raw$week_start)) stop("No valid weeks, or duplicate ISO weeks.")
  weekly <- tibble(
    week_start = raw$week_start,
    iso_year = raw$Jaar,
    iso_week = raw$Week,
    n_laboratories = raw$Aantal_laboratoria,
    Influenza = rowSums(as.data.frame(raw[, influenza_columns, drop = FALSE])),
    RSV = raw$VIRRS,
    HMPV = raw$VIRHMPV,
    Parainfluenza = rowSums(as.data.frame(raw[, piv_columns, drop = FALSE]))
  )
  weekly <- weekly %>% arrange(week_start)
  # Missing weeks/counts are retained as NA rather than treated as zero.
  weekly <- tibble(week_start = seq(min(weekly$week_start), max(weekly$week_start), by = "week")) %>%
    left_join(mutate(weekly, source_week_present = TRUE), by = "week_start") %>%
    mutate(source_week_present = coalesce(source_week_present, FALSE),
           iso_year = as.integer(format(week_start + 3, "%Y")),
           iso_week = as.integer(format(week_start, "%V")))
  for (path in c(data_dir, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(weekly, file.path(data_dir, "pathogens_virological_weekly.csv"), row.names = FALSE, na = "NA")
  if (nrow(omitted)) utils::write.csv(omitted, file.path(data_dir, "pathogens_virological_invalid_ISO_weeks.csv"), row.names = FALSE)
  
  # Approximate monthly totals: assign each entire week to the month of its
  # Thursday. Exact calendar-month counts require individual detection dates.
  monthly <- weekly %>% mutate(month = as.Date(format(week_start + 3, "%Y-%m-01"))) %>%
    group_by(month) %>% summarise(across(all_of(types), ~sum(.x)),
                                  n_weeks_present = sum(source_week_present), .groups = "drop")
  monthly$expected_weeks <- vapply(monthly$month, function(m) {
    first <- as.Date(m, origin = "1970-01-01")
    last <- seq(first, by = "month", length.out = 2)[2] - 1
    as.integer(sum(format(seq(first, last, by = "day"), "%u") == "4"))
  }, integer(1))
  monthly$complete_month <- monthly$n_weeks_present == monthly$expected_weeks
  utils::write.csv(monthly, file.path(data_dir, "pathogens_virological_monthly.csv"), row.names = FALSE, na = "NA")
  long <- monthly %>% pivot_longer(all_of(types), names_to = "component", values_to = "detections") %>%
    mutate(component = factor(component, levels = types))
  long$detections[!long$complete_month] <- NA_real_
  
  if (!any(is.finite(long$detections))) stop("No complete monthly counts to plot.")
  totals <- data.frame(
    pathogen = types,
    detections = vapply(weekly[types], sum, numeric(1)),
    weeks_with_counts = vapply(weekly[types], function(x) sum(!is.na(x)), integer(1)),
    row.names = NULL
  )
  cat("\nVIROLOGICAL DETECTIONS OVER THE AVAILABLE PERIOD\n")
  cat("Weeks:", format(min(weekly$week_start)), "to", format(max(weekly$week_start)), "\n")
  print(totals, row.names = FALSE)
  stem <- "Figure_S7_virological_monthly_full_series"
  utils::write.csv(totals, file.path(figure_dir, paste0(stem, "_totals.csv")),
                   row.names = FALSE, na = "NA")
  utils::write.csv(long, file.path(figure_dir, paste0(stem, "_plotted_data.csv")),
                   row.names = FALSE, na = "NA")
  
  colours <- c(Influenza = "#0072B2", RSV = "#D55E00", HMPV = "#009E73",
               Parainfluenza = "#CC79A7")
  panels <- c(Influenza = "A  Influenza", RSV = "B  RSV", HMPV = "C  HMPV",
              Parainfluenza = "D  Parainfluenza")
  lancet_theme <- theme_minimal(base_size = 9, base_family = "sans") +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.2),
          axis.line = element_line(colour = "grey50", linewidth = 0.2),
          axis.ticks = element_line(colour = "grey50", linewidth = 0.2),
          axis.ticks.length = grid::unit(1, "mm"),
          axis.text = element_text(size = 8, colour = "grey20"),
          axis.title = element_text(size = 9, colour = "grey20"),
          axis.title.x = element_text(margin = margin(t = 6)),
          axis.title.y = element_text(margin = margin(r = 6)),
          strip.text = element_text(size = 10, face = "bold", hjust = 0, colour = "grey20"),
          panel.spacing = grid::unit(4, "mm"), plot.margin = margin(5, 4, 4, 4),
          legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 9), legend.key.width = grid::unit(8, "mm"))
  
  # Dotted vertical lines delineate half-years, as in script 17.
  period_boundaries <- seq(
    as.Date(paste0(format(min(long$month), "%Y"), "-01-01")),
    as.Date(paste0(format(max(long$month), "%Y"), "-07-01")),
    by = "6 months"
  )
  period_boundaries <- period_boundaries[
    period_boundaries >= min(long$month) & period_boundaries <= max(long$month)
  ]
  plot_time <- ggplot(long, aes(month, detections, colour = component)) +
    geom_vline(xintercept = as.numeric(period_boundaries),
               linetype = "dotted", colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.55, na.rm = TRUE) +
    facet_wrap(~component, ncol = 2, scales = "free_y", labeller = as_labeller(panels)) +
    scale_colour_manual(values = colours, guide = "none") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y",
                 limits = range(long$month), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(x = "Year", y = "Reported virological detections per month") + lancet_theme
  
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot_time, device = pdf_device,
         width = 180, height = 140, units = "mm", bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot_time, dpi = 600,
         width = 180, height = 140, units = "mm", bg = "white")
  if (show_figures) print(plot_time)
  
  notes <- c(
    paste("Source:", normalizePath(input_file, winslash = "/")),
    paste("Source MD5:", unname(tools::md5sum(input_file))), paste("Sheet:", input_sheet),
    paste("Available weeks:", min(weekly$week_start), "to", max(weekly$week_start)),
    "Laboratory surveillance counts, without age or hospital diagnosis stratification.",
    "Influenza: VIRINFA + VIRINFB. RSV: VIRRS. HMPV: VIRHMPV.",
    paste("Parainfluenza:", paste(piv_columns, collapse = " + ")),
    "PIV includes types 1-4 plus untyped/non-typeable detections.",
    "Coinfection status is unavailable; detections cannot be deduplicated by patient.",
    "Whole weeks assigned to their Thursday's month; monthly totals are approximate.",
    "Missing component counts propagate to pooled and monthly counts; no missing counts are set to zero.",
    "Incomplete months remain in the data CSV but are omitted from the figure.",
    "Counts unadjusted for month length, reporting laboratory numbers and testing intensity.",
    "Dotted vertical lines mark 1 January and 1 July each year.",
    "Panels use different vertical scales. No model has been fitted."
  )
  writeLines(notes, file.path(figure_dir, paste0(stem, "_notes.txt")))
  message("\nSaved weekly/monthly data in:\n", normalizePath(data_dir, winslash = "/"),
          "\nSaved figure and summaries in:\n", normalizePath(figure_dir, winslash = "/"))
})()
