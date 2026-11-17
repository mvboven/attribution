# PIV1-PIV4 in the Dutch virological laboratory surveillance.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
# Standalone: source from the project root or scripts/; no Stan fit required.
# These are laboratory detections. Coinfection status is unavailable.

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
  profile_years <- getOption("attribution.piv_profile_years", 2013:2019)
  figure_dir <- getOption("attribution.supplement_dir", file.path(root, "figures", "manuscript", "supplement"))
  data_dir <- getOption("attribution.piv_data_dir", file.path(root, "data"))
  show_figures <- isTRUE(getOption("attribution.supplement_show", FALSE))
  types <- c("PIV1", "PIV2", "PIV3", "PIV4")
  source_columns <- paste0("VIRPAR", 1:4)
  if (!file.exists(input_file)) stop("Virological workbook not found: ", input_file)
  if (!is.numeric(profile_years) || !length(profile_years) || anyNA(profile_years) ||
      any(!is.finite(profile_years)) || any(profile_years != floor(profile_years))) {
    stop("Use integer profile years, e.g. 2013:2019.")
  }
  profile_years <- sort(unique(as.integer(profile_years)))
  
  # Only the four type-specific columns are used; untyped PIV is excluded.
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
  weekly <- raw %>% select(week_start, Jaar, Week, Aantal_laboratoria, all_of(source_columns))
  names(weekly) <- c("week_start", "iso_year", "iso_week", "n_laboratories", types)
  weekly <- weekly %>% arrange(week_start)
  # Missing weeks/counts are retained as NA rather than treated as zero.
  weekly <- tibble(week_start = seq(min(weekly$week_start), max(weekly$week_start), by = "week")) %>%
    left_join(mutate(weekly, source_week_present = TRUE), by = "week_start") %>%
    mutate(source_week_present = coalesce(source_week_present, FALSE),
           iso_year = as.integer(format(week_start + 3, "%Y")),
           iso_week = as.integer(format(week_start, "%V")))
  for (path in c(data_dir, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(weekly, file.path(data_dir, "piv_types_virological_weekly.csv"), row.names = FALSE, na = "NA")
  if (nrow(omitted)) utils::write.csv(omitted, file.path(data_dir, "piv_types_invalid_ISO_weeks.csv"), row.names = FALSE)
  
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
  utils::write.csv(monthly, file.path(data_dir, "piv_types_virological_monthly.csv"), row.names = FALSE, na = "NA")
  long <- monthly %>% pivot_longer(all_of(types), names_to = "type", values_to = "detections") %>%
    mutate(type = factor(type, levels = types), year = as.integer(format(month, "%Y")),
           month_number = as.integer(format(month, "%m")))
  long$detections[!long$complete_month] <- NA_real_
  selected <- long %>% filter(year %in% profile_years)
  if (nrow(selected) != 12L * length(types) * length(profile_years) ||
      !setequal(unique(selected$year), profile_years) || anyNA(selected$detections)) {
    stop("The calendar profile requires complete data for the requested years. CSVs have been saved. ",
         "Set attribution.piv_profile_years to complete years.")
  }
  profile <- selected %>% group_by(type, month_number) %>%
    summarise(mean_monthly_detections = mean(detections), n_years = n_distinct(year), .groups = "drop")
  utils::write.csv(profile, file.path(figure_dir, "Figure_S6_PIV_calendar_month_profile.csv"), row.names = FALSE)
  printed <- selected %>% group_by(type) %>% summarise(detections = sum(detections), .groups = "drop") %>%
    mutate(percent_of_typed_detections = if (sum(detections) > 0) 100 * detections / sum(detections) else NA_real_)
  cat("\nPIV TYPE COUNTS, YEARS ", paste(profile_years, collapse = ", "), "\n", sep = "")
  print(as.data.frame(printed), row.names = FALSE)
  utils::write.csv(printed, file.path(figure_dir, "Figure_S6_PIV_type_totals.csv"), row.names = FALSE)
  
  colours <- c(PIV1 = "#0072B2", PIV2 = "#D55E00", PIV3 = "#009E73", PIV4 = "#CC79A7")
  lines <- c(PIV1 = "solid", PIV2 = "dashed", PIV3 = "solid", PIV4 = "dotdash")
  labels <- c(PIV1 = "PIV-1", PIV2 = "PIV-2", PIV3 = "PIV-3", PIV4 = "PIV-4")
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
  
  # Main figure: common scale so the relative numbers remain visible.
  plot_profile <- ggplot(profile, aes(month_number, mean_monthly_detections,
                                      colour = type, linetype = type, group = type)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.3) +
    scale_colour_manual(values = colours, breaks = types, labels = labels) +
    scale_linetype_manual(values = lines, breaks = types, labels = labels) +
    scale_x_continuous(breaks = 1:12, labels = month.abb, expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(x = "Month", y = "Mean reported PIV detections per month") + lancet_theme
  # Companion figure: full available period; each type has its own y-axis scale.
  # Dotted vertical lines mark 1 January and 1 July of each displayed year.
  period_boundaries <- seq(
    as.Date(paste0(format(min(long$month), "%Y"), "-01-01")),
    as.Date(paste0(format(max(long$month), "%Y"), "-07-01")),
    by = "6 months"
  )
  period_boundaries <- period_boundaries[
    period_boundaries >= min(long$month) & period_boundaries <= max(long$month)
  ]
  panels <- setNames(paste(LETTERS[1:4], unname(labels), sep = "  "), types)
  plot_time <- ggplot(long, aes(month, detections, colour = type)) +
    geom_vline(xintercept = as.numeric(period_boundaries),
               linetype = "dotted", colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.55, na.rm = TRUE) +
    facet_wrap(~type, ncol = 2, scales = "free_y", labeller = as_labeller(panels)) +
    scale_colour_manual(values = colours, guide = "none") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y",
                 limits = range(long$month), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(x = "Year", y = "Reported PIV detections per month") + lancet_theme
  save_figure <- function(p, stem, height) {
    pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggsave(file.path(figure_dir, paste0(stem, ".pdf")), p, device = pdf_device,
           width = 180, height = height, units = "mm", bg = "white")
    ggsave(file.path(figure_dir, paste0(stem, ".png")), p, dpi = 600,
           width = 180, height = height, units = "mm", bg = "white")
    if (show_figures) print(p)
  }
  save_figure(plot_profile, "Figure_S6_PIV_types_by_month", 110)
  save_figure(plot_time, "Figure_S6_PIV_types_full_series", 140)
  notes <- c(
    paste("Source:", normalizePath(input_file, winslash = "/")),
    paste("Source MD5:", unname(tools::md5sum(input_file))), paste("Sheet:", input_sheet),
    paste("Available weeks:", min(weekly$week_start), "to", max(weekly$week_start)),
    paste("Calendar-profile years:", paste(profile_years, collapse = ", ")),
    "Laboratory surveillance counts, without age or hospital diagnosis stratification.",
    "Included: VIRPAR1, VIRPAR2, VIRPAR3, VIRPAR4. Untyped/non-typeable PIV excluded.",
    "Coinfection status is unavailable. Counts cannot be restricted to single infections or deduplicated by patient.",
    "Whole weeks assigned to their Thursday's month; monthly totals are approximate.",
    "Calendar profile: mean monthly counts across selected complete years.",
    "Counts unadjusted for month length, reporting laboratory numbers and testing intensity.",
    "Dotted vertical lines in the full-series figure mark 1 January and 1 July each year.",
    "Full-series figure uses different vertical scales by type. No model has been fitted."
  )
  writeLines(notes, file.path(figure_dir, "Figure_S6_PIV_notes.txt"))
  message("\nSaved weekly/monthly data in:\n", normalizePath(data_dir, winslash = "/"),
          "\nSaved figures and summaries in:\n", normalizePath(figure_dir, winslash = "/"))
})()
