# Read and align hospital admissions and virological data.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
require_objects(c("diagnosis_scope", "hospital_file", "age_groups"),
                instruction = "Source scripts/0_source.R first.")
message("Reading hospital and virological data for ", diagnosis_scope, ".")
analysis_calendar <- tibble::as_tibble(make_iso_calendar(first_year, last_year))
# load hospital data
#
# The first column in the CSV is an exported row-number column with an empty
# header. readr repairs this to a name such as "...1"; it is removed below.
#
# Read all columns as character initially because suppressed cells are stored
# as "1-4".
if (!file.exists(hospital_file)) {
  stop(
    "Hospital file not found: ",
    hospital_file
  )
}

hospital_raw <- readr::read_csv(
  hospital_file,
  col_types = readr::cols(
    .default = readr::col_character()
  ),
  show_col_types = FALSE
) %>%
  select(
    -matches("^\\.\\.\\.[0-9]+$")
  )

# check required columns
required_hospital_columns <- c(
  "week_start",
  "age_group",
  "pathogen",
  "n_admissions_v"
)

missing_hospital_columns <- setdiff(
  required_hospital_columns,
  names(hospital_raw)
)

if (length(missing_hospital_columns) > 0L) {
  stop(
    "Required hospital-data columns are missing: ",
    paste(
      missing_hospital_columns,
      collapse = ", "
    )
  )
}

# standardise the variables used below
hospital_raw <- hospital_raw %>%
  mutate(
    week_start = as.Date(week_start),
    age_group = trimws(
      as.character(age_group)
    ),
    pathogen = trimws(
      as.character(pathogen)
    )
  )

if (anyNA(hospital_raw$week_start)) {
  stop(
    "Some hospital-data week_start values could not be converted to dates."
  )
}

check_unique(hospital_raw, c("week_start", "age_group", "pathogen"), "Hospital extract")
if (any(as.integer(format(hospital_raw$week_start, "%u")) != 1L)) {
  stop("Hospital week_start values must be ISO Mondays.")
}
hospital_raw <- hospital_raw %>%
  mutate(
    number = convert_hospital_counts(
      n_admissions_v
    )
  )

# Retain the supplied code's sum of non-COVID diagnosis categories.
# IMPORTANT: the export explanation states that coinfections are counted under
# each pathogen. This category sum is therefore not a deduplicated admission
# total. See data/README.md; deduplicated totals require additional source data.
#
# COVID19 categories are deliberately excluded. They are not part of the
# respiratory-virus attribution model and the analysis period is pre-COVID.
rti_components <- c(
  "Background",
  "Influenza_single",
  "Influenza_coinfection",
  "RSV_single",
  "RSV_coinfection",
  "HMPV_single",
  "HMPV_coinfection",
  "PIV_single",
  "PIV_coinfection"
)

missing_rti_components <- setdiff(
  rti_components,
  unique(hospital_raw$pathogen)
)

if (length(missing_rti_components) > 0L) {
  stop(
    "Expected RTI components are missing from the hospital file: ",
    paste(
      missing_rti_components,
      collapse = ", "
    )
  )
}

coinfection_rows <- hospital_raw %>%
  filter(week_start %in% analysis_calendar$week_start,
         age_group %in% age_group_lookup$source_age_group,
         pathogen %in% c("Influenza_coinfection", "RSV_coinfection", "HMPV_coinfection", "PIV_coinfection"))
if (any(coinfection_rows$number > 0)) {
  message("Retaining the documented minor coinfection overlap in the hospital category sum; ",
          "see data/README.md.")
}


missing_source_age_groups <- setdiff(
  age_group_lookup$source_age_group,
  unique(hospital_raw$age_group)
)

if (length(missing_source_age_groups) > 0L) {
  stop(
    "Expected source age groups are missing from the hospital file: ",
    paste(
      missing_source_age_groups,
      collapse = ", "
    )
  )
}

# Hospital-admission data in long format
hospital_data_long <- hospital_raw %>%
  rename(
    source_age_group = age_group
  ) %>%
  filter(
    pathogen %in% rti_components
  ) %>%
  inner_join(
    age_group_lookup,
    by = "source_age_group"
  ) %>%
  group_by(
    week_start,
    age_group
  ) %>%
  summarise(
    number = sum(
      number,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  inner_join(
    analysis_calendar %>%
      select(
        week_start,
        year,
        week
      ),
    by = "week_start"
  ) %>%
  select(
    year,
    week,
    age_group,
    number
  ) %>%
  arrange(
    year,
    week,
    factor(
      age_group,
      levels = age_groups
    )
  )

if (nrow(hospital_data_long) == 0L) {
  stop(
    "No hospital admissions remain after applying the selected filters."
  )
}

missing_hospital_cells <- anti_join(
  tidyr::crossing(analysis_calendar %>% distinct(year, week), age_group = age_groups),
  hospital_data_long, by = c("year", "week", "age_group"))
if (nrow(missing_hospital_cells) && !hospital_absent_is_zero) {
  stop(nrow(missing_hospital_cells), " hospital week/age combinations are absent from the extract.")
}
if (nrow(missing_hospital_cells)) {
  message("Zero-filling ", nrow(missing_hospital_cells), " absent hospital week/age combinations.")
}

# Explicitly add zero counts for week-age combinations that are absent from
# the extract. This keeps the hospital time series aligned with the complete
# ISO-week calendar used by the virological data and by Stan.
hospital_data <- tidyr::crossing(
  analysis_calendar %>%
    distinct(
      year,
      week
    ),
  age_group = age_groups
) %>%
  left_join(
    hospital_data_long,
    by = c(
      "year",
      "week",
      "age_group"
    )
  ) %>%
  mutate(
    number = replace_na(
      number,
      0L
    )
  ) %>%
  pivot_wider(
    names_from = age_group,
    values_from = number
  ) %>%
  select(
    year,
    week,
    all_of(age_groups)
  ) %>%
  arrange(
    year,
    week
  )

# observed totals used as denominators in the attribution table
observed_totals_by_age <- hospital_data_long %>%
  group_by(
    age_group
  ) %>%
  summarise(
    observed_total = sum(
      number
    ),
    .groups = "drop"
  ) %>%
  mutate(
    age_group = factor(
      age_group,
      levels = age_groups
    )
  ) %>%
  arrange(
    age_group
  )

# checks of the prepared hospital data
stopifnot(
  nrow(hospital_data) == 365L,
  all(age_groups %in% names(hospital_data)),
  !anyNA(hospital_data[, age_groups]),
  all(
    unlist(
      hospital_data[, age_groups]
    ) >= 0L
  )
)

# load virological data
# The virological workbook is stored with the other required inputs.
if (!file.exists(virology_file)) stop("Virological file not found: ", virology_file)
virologische_weekstaten_20240322 <- read_excel(virology_file)
required_virology <- c("Jaar", "Week", "VIRINFA", "VIRINFB", "VIRRS", "VIRHMPV",
                       "VIRPAR1", "VIRPAR2", "VIRPAR3", "VIRPAR4", "VIRPARANGETYP", "VIRPARANTYP")
check_columns(virologische_weekstaten_20240322, required_virology, "Virological data")

# weekly virological data
vir_data <- virologische_weekstaten_20240322 %>%
  rename(
    year = Jaar,
    week = Week
  ) %>%
  mutate(
    year = as.integer(year),
    week = as.integer(week)
  ) %>%
  filter(
    year >= first_year,
    year <= last_year
  )

check_unique(vir_data, c("year", "week"), "Virological data")
if (anyNA(vir_data$year) || anyNA(vir_data$week)) stop("Invalid virological year/week values.")
unused_virology_weeks <- anti_join(vir_data, analysis_calendar, by = c("year", "week"))
if (nrow(unused_virology_weeks)) {
  message("Ignoring ", nrow(unused_virology_weeks), " virological rows outside the valid study ISO calendar.")
}
vir_data <- semi_join(vir_data, analysis_calendar, by = c("year", "week"))
missing_virus_weeks <- anti_join(analysis_calendar, vir_data, by = c("year", "week"))
if (nrow(missing_virus_weeks)) stop("Virological data are missing ", nrow(missing_virus_weeks), " study weeks.")

# merge both datasets 
data <- analysis_calendar %>%
  left_join(
    hospital_data,
    by = c("year", "week")
  ) %>%
  left_join(
    vir_data,
    by = c("year", "week")
  ) %>%
  arrange(week_start)

# Combine subcategories of PIV and influenza
data <- data %>%
  mutate(
    VIRPIV = VIRPAR1 + VIRPAR2 + VIRPAR3 +
      VIRPAR4 + VIRPARANGETYP + VIRPARANTYP,
    
    VIRINF = VIRINFA + VIRINFB
  )


check_unique(data, c("year", "week"), "Merged study data")
data_diagnosis_scope <- diagnosis_scope
run_signature <- NULL  # Reprepare Stan inputs after every successful data read.
export_table(observed_totals_by_age, "observed_totals_by_age", analysis_output_dir)
message("Aligned ", nrow(data), " ISO weeks; ", sum(observed_totals_by_age$observed_total),
        " summed RTI diagnosis-category counts across the eight age groups.")

# Export the minimal weekly input for the public attribution analysis.
# Run this script once for each diagnosis_scope to produce separate files.
# Only the study weeks, model age groups and relevant virology counts are kept.
# Influenza A and B are retained separately, so either A alone or their sum can
# be analysed. PIV subtypes have already been combined into VIRPIV above.
# year/week identify ISO weeks; week_start and VIRINF can be reconstructed.
# The hospital counts retain the original suppression conversion, aggregation
# and zero-filling. Original hospital categories are not part of this export.
#
# Optionally set compact_data_dir before sourcing either reader to change the
# default location. Run from the repository root when using this default.
if (!exists("compact_data_dir", inherits = FALSE)) {
  compact_data_dir <- file.path("data", "compact")
}
compact_data_file <- file.path(
  compact_data_dir,
  paste0("attribution_", diagnosis_scope, ".csv")
)

compact_virology_columns <- c(
  "VIRINFA", "VIRINFB", "VIRRS", "VIRHMPV", "VIRPIV"
)
compact_data <- data %>%
  select(year, week, all_of(age_groups), all_of(compact_virology_columns))

compact_count_columns <- c(age_groups, compact_virology_columns)
if (!all(vapply(compact_data[compact_count_columns], is.numeric, logical(1))) ||
    any(!is.finite(as.matrix(compact_data[compact_count_columns]))) ||
    any(as.matrix(compact_data[compact_count_columns]) < 0)) {
  stop("Compact-data counts must be numeric, finite and non-negative.")
}

dir.create(compact_data_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(compact_data, compact_data_file)
message("Wrote ", nrow(compact_data), " weeks to ", compact_data_file, ".")
