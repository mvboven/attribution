# Read the minimal weekly hospital-admission and virological data.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
# Source scripts/0_source.R first, then use this script instead of 1_read_data.R.
# The matching diagnosis_scope, age_groups and study years must be configured.
# This reader needs neither the original hospital CSV nor the virology workbook.
require_objects(
  c("diagnosis_scope", "age_groups", "first_year", "last_year", "analysis_output_dir"),
  instruction = "Source scripts/0_source.R first."
)

if (!exists("compact_data_dir", inherits = FALSE)) {
  compact_data_dir <- file.path("data", "compact")
}
compact_data_file <- file.path(
  compact_data_dir,
  paste0("attribution_", diagnosis_scope, ".csv")
)
if (!file.exists(compact_data_file)) {
  stop("Compact data file not found: ", compact_data_file,
       ". Generate it with 1_read_data.R using the original inputs.")
}
message("Reading compact hospital and virological data for ", diagnosis_scope, ".")

analysis_calendar <- tibble::as_tibble(make_iso_calendar(first_year, last_year))
compact_virology_columns <- c(
  "VIRINFA", "VIRINFB", "VIRRS", "VIRHMPV", "VIRPIV"
)
compact_count_columns <- c(age_groups, compact_virology_columns)
compact_data <- readr::read_csv(
  compact_data_file,
  col_types = readr::cols(
    .default = readr::col_double(),
    year = readr::col_integer(),
    week = readr::col_integer()
  ),
  show_col_types = FALSE,
  name_repair = "check_unique"
)
check_columns(compact_data, c("year", "week", compact_count_columns), "Compact data")
if (nrow(readr::problems(compact_data))) {
  stop("Some compact-data values could not be read as numbers; see readr::problems(compact_data).")
}
check_unique(compact_data, c("year", "week"), "Compact data")
if (anyNA(compact_data$year) || anyNA(compact_data$week)) {
  stop("Invalid compact-data year/week values.")
}

missing_compact_weeks <- anti_join(
  analysis_calendar, compact_data, by = c("year", "week")
)
unexpected_compact_weeks <- anti_join(
  compact_data, analysis_calendar, by = c("year", "week")
)
if (nrow(missing_compact_weeks) || nrow(unexpected_compact_weeks)) {
  stop("Compact data do not match the configured ISO calendar: ",
       nrow(missing_compact_weeks), " missing and ",
       nrow(unexpected_compact_weeks), " unexpected weeks.")
}
if (any(!is.finite(as.matrix(compact_data[compact_count_columns]))) ||
    any(as.matrix(compact_data[compact_count_columns]) < 0)) {
  stop("Compact-data counts must be finite and non-negative.")
}

# Reconstruct the prepared inputs and totals used by the attribution analysis.
# Separate PIV subtypes and original hospital categories cannot be recovered
# from this aggregate dataset.
data <- analysis_calendar %>%
  left_join(
    compact_data %>%
      select(year, week, all_of(compact_count_columns)),
    by = c("year", "week")
  ) %>%
  mutate(VIRINF = VIRINFA + VIRINFB) %>%
  arrange(week_start)

hospital_data <- data %>%
  select(year, week, all_of(age_groups))
hospital_data_long <- hospital_data %>%
  pivot_longer(
    cols = all_of(age_groups),
    names_to = "age_group",
    values_to = "number"
  ) %>%
  arrange(year, week, factor(age_group, levels = age_groups))
observed_totals_by_age <- hospital_data_long %>%
  group_by(age_group) %>%
  summarise(observed_total = sum(number), .groups = "drop") %>%
  mutate(age_group = factor(age_group, levels = age_groups)) %>%
  arrange(age_group)
vir_data <- data %>%
  select(year, week, all_of(compact_virology_columns), VIRINF)

check_unique(data, c("year", "week"), "Merged study data")
data_diagnosis_scope <- diagnosis_scope
run_signature <- NULL  # Reprepare Stan inputs after every successful data read.
export_table(observed_totals_by_age, "observed_totals_by_age", analysis_output_dir)
message("Aligned ", nrow(data), " ISO weeks; ", sum(observed_totals_by_age$observed_total),
        " summed RTI diagnosis-category counts across the eight age groups.")
