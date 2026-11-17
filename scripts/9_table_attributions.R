# Overall, winter, and complementary other-weeks attribution tables.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
#
# Within each posterior draw:
# allocated admissions = observed admissions * pathogen intensity / total intensity.
#
# Incidence for every period:
# allocated admissions in selected weeks / population person-years * 100,000.
#
# Percentages for every period:
# allocated admissions in selected weeks / observed admissions in those weeks * 100.
#
# Winter and other-weeks incidence are contributions to annual incidence.
# They are NOT annualised by dividing by the fraction of weeks selected.

require_fit()

message(
  "Calculating attribution tables, credible intervals, ",
  "and observed-admission denominators."
)

# -------------------------------------------------------------------------
# Display order and labels
# -------------------------------------------------------------------------

# Display order only. Stan indices remain mapped through pathogen_lookup.
pathogen_table_order <- c(
  "Influenza",
  "RSV",
  "HMPV",
  "Parainfluenza"
)

manuscript_pathogen_labels <- c(
  "Influenza" = "Influenza",
  "RSV" = "RSV",
  "HMPV" = "Human metapneumovirus",
  "Parainfluenza" = "Parainfluenza virus"
)

manuscript_age_labels <- c(
  "0-2m"   = "0\u20132 months",
  "3-5m"   = "3\u20135 months",
  "6-11m"  = "6\u201311 months",
  "1-2y"   = "1\u20132 years",
  "3-4y"   = "3\u20134 years",
  "65-74y" = "65\u201374 years",
  "75-84y" = "75\u201384 years",
  "85+y"   = "\u226585 years"
)

stopifnot(
  first_year == 2013L,
  last_year == 2019L,
  length(age_groups) == num_types,
  !anyDuplicated(age_groups),
  all(age_groups %in% names(manuscript_age_labels)),
  nrow(pathogen_lookup) == num_path,
  !anyDuplicated(as.character(pathogen_lookup$pathogen)),
  setequal(
    as.character(pathogen_lookup$pathogen),
    pathogen_table_order
  ),
  nrow(y_data) == n_data,
  ncol(y_data) == num_types,
  nrow(data) == n_data,
  all(is.finite(y_data)),
  all(y_data >= 0)
)

# -------------------------------------------------------------------------
# Population denominators
# -------------------------------------------------------------------------

# CBS mean populations for 2013-2019.
# Each year's population was approximated by the mean of 1 January that year
# and 1 January the next year. Infants aged 0 years are divided 3/12, 3/12, 6/12.
# Sources: CBS StatLine 7461eng and 85510ENG.
# Match values by age label.

population_reference <- c(
  "0-2m"   = 42814.785714,
  "3-5m"   = 42814.785714,
  "6-11m"  = 85629.571429,
  "1-2y"   = 349337.000000,
  "3-4y"   = 358851.285714,
  "65-74y" = 1793692.500000,
  "75-84y" = 968578.000000,
  "85+y"   = 356764.785714
)

num_analysis_years <- last_year - first_year + 1L

population_denominators <- tibble(
  age_group = age_groups,
  mean_population = unname(population_reference[age_groups]),
  person_years = mean_population * num_analysis_years
)

stopifnot(
  num_analysis_years == 7L,
  all(is.finite(population_denominators$person_years)),
  all(population_denominators$person_years > 0)
)

# -------------------------------------------------------------------------
# Period definitions
# -------------------------------------------------------------------------

periods <- make_periods(
  data$week,
  winter_start_week,
  winter_end_week
)

stopifnot(
  all(c("overall", "winter", "other_weeks") %in% names(periods))
)

periods <- periods[c("overall", "winter", "other_weeks")]

stopifnot(
  all(lengths(periods) > 0L),
  all(vapply(
    periods,
    function(indices) {
      all(is.finite(indices)) &&
        all(indices == as.integer(indices)) &&
        all(indices >= 1L & indices <= n_data) &&
        !anyDuplicated(indices)
    },
    logical(1)
  )),
  setequal(periods$overall, seq_len(n_data)),
  length(intersect(periods$winter, periods$other_weeks)) == 0L,
  setequal(
    c(periods$winter, periods$other_weeks),
    periods$overall
  )
)

period_calendar <- data %>%
  transmute(
    week_start,
    ISOweek,
    year,
    week,
    period = if_else(
      seq_len(n_data) %in% periods$winter,
      "winter",
      "other_weeks"
    )
  )

export_table(
  period_calendar,
  "attribution_period_calendar",
  analysis_output_dir
)

# Exact observed-admission denominators, calculated directly from the data.
attribution_denominators <- bind_rows(
  lapply(names(periods), function(selected_period) {
    tibble(
      period = selected_period,
      age_group = age_groups,
      n_weeks = length(periods[[selected_period]]),
      observed_total = as.numeric(
        colSums(
          y_data[periods[[selected_period]], , drop = FALSE]
        )
      )
    )
  })
) %>%
  mutate(
    period = factor(period, levels = names(periods)),
    age_group = factor(age_group, levels = age_groups)
  ) %>%
  arrange(period, age_group)

observed_partition_check <- attribution_denominators %>%
  select(age_group, period, observed_total) %>%
  pivot_wider(
    names_from = period,
    values_from = observed_total
  )

stopifnot(
  all(
    observed_partition_check$overall ==
      observed_partition_check$winter +
      observed_partition_check$other_weeks
  )
)

# -------------------------------------------------------------------------
# Observed-count allocations within each posterior draw
# -------------------------------------------------------------------------

lambda_variables <- unlist(
  lapply(seq_len(num_types), function(a) {
    sprintf("lambda[%d,%d]", a, seq_len(n_data))
  }),
  use.names = FALSE
)

lambda_draw_matrix <- get_draw_matrix(
  fit_hmc,
  lambda_variables
)

stopifnot(
  all(lambda_variables %in% colnames(lambda_draw_matrix))
)

lambda_draw_matrix <- lambda_draw_matrix[
  ,
  lambda_variables,
  drop = FALSE
]

if (
  any(!is.finite(lambda_draw_matrix)) ||
  any(lambda_draw_matrix <= 0)
) {
  stop("Total fitted intensities must be finite and positive.")
}

combined_contribution <- matrix(
  0,
  nrow = nrow(lambda_draw_matrix),
  ncol = ncol(lambda_draw_matrix)
)

pathogen_total_list <- vector("list", num_path)

for (p in seq_len(num_path)) {
  
  contribution_names <- unlist(
    lapply(seq_len(num_types), function(a) {
      sprintf(
        "lambda_contributions[%d,%d,%d]",
        p,
        a,
        seq_len(n_data)
      )
    }),
    use.names = FALSE
  )
  
  contribution <- get_draw_matrix(
    fit_hmc,
    contribution_names
  )
  
  stopifnot(
    all(contribution_names %in% colnames(contribution))
  )
  
  contribution <- contribution[
    ,
    contribution_names,
    drop = FALSE
  ]
  
  if (!identical(dim(contribution), dim(lambda_draw_matrix))) {
    stop("Posterior draw or variable counts differ across variables.")
  }
  
  if (
    any(!is.finite(contribution)) ||
    any(contribution < 0)
  ) {
    stop("Pathogen intensities must be finite and non-negative.")
  }
  
  combined_contribution <- combined_contribution + contribution
  
  pathogen_total_list[[p]] <- bind_rows(
    lapply(seq_len(num_types), function(a) {
      
      columns <- seq.int(
        from = (a - 1L) * n_data + 1L,
        to = a * n_data
      )
      
      allocated <- allocate_observed(
        contribution[, columns, drop = FALSE],
        lambda_draw_matrix[, columns, drop = FALSE],
        y_data[, a]
      )
      
      bind_rows(
        lapply(names(periods), function(selected_period) {
          
          indices <- periods[[selected_period]]
          observed <- sum(y_data[indices, a])
          
          total <- rowSums(
            allocated[, indices, drop = FALSE]
          )
          
          tibble(
            draw = seq_along(total),
            period = selected_period,
            pathogen_index = p,
            pathogen = as.character(pathogen_lookup$pathogen[p]),
            age_index = a,
            age_group = age_groups[a],
            n_weeks = length(indices),
            observed_total = observed,
            attributed_total = total,
            
            fraction_observed = if (observed > 0) {
              100 * total / observed
            } else {
              NA_real_
            },
            
            attributed_per_year = total / num_analysis_years,
            attributed_per_week = total / length(indices),
            
            # Same population person-years for all three periods.
            incidence_per_100000 = (
              100000 * total /
                population_denominators$person_years[a]
            )
          )
        })
      )
    })
  )
  
  rm(contribution)
}

if (
  any(combined_contribution > lambda_draw_matrix * (1 + 1e-8))
) {
  stop(
    "The sum of pathogen contributions exceeds the fitted total ",
    "in some draw/week/age cells."
  )
}

attribution_draws <- bind_rows(pathogen_total_list)

# Combined attribution to the four pathogens.
# Uses existing posterior allocations; no additional fitting is needed.

four_pathogen_draws <- attribution_draws %>%
  group_by(draw, period, age_index, age_group) %>%
  summarise(
    n_pathogens = n_distinct(pathogen_index),
    n_rows = n(),
    n_denominators = n_distinct(observed_total),
    observed_total = first(observed_total),
    attributed_total = sum(attributed_total),
    .groups = "drop"
  )

# Each draw/period/age group must contain exactly four pathogens
# with the same observed admission denominator.
stopifnot(
  all(four_pathogen_draws$n_pathogens == 4L),
  all(four_pathogen_draws$n_rows == 4L),
  all(four_pathogen_draws$n_denominators == 1L),
  all(is.finite(four_pathogen_draws$observed_total)),
  all(four_pathogen_draws$observed_total > 0),
  all(is.finite(four_pathogen_draws$attributed_total)),
  all(four_pathogen_draws$attributed_total >= 0),
  all(
    four_pathogen_draws$attributed_total <=
      four_pathogen_draws$observed_total * (1 + 1e-8)
  )
)

four_pathogen_draws <- four_pathogen_draws %>%
  mutate(
    four_pathogen_pct = 100 * attributed_total / observed_total,
    background_pct = 100 - four_pathogen_pct
  )

# Summarise the combined fraction, rather than adding
# separately calculated pathogen medians or interval limits.
four_pathogen_summary <- four_pathogen_draws %>%
  group_by(period, age_index, age_group) %>%
  summarise(
    observed_total = first(observed_total),
    
    four_pathogen_pct_median = median(four_pathogen_pct),
    four_pathogen_pct_lower = unname(
      quantile(four_pathogen_pct, 0.025)
    ),
    four_pathogen_pct_upper = unname(
      quantile(four_pathogen_pct, 0.975)
    ),
    
    background_pct_median = median(background_pct),
    background_pct_lower = unname(
      quantile(background_pct, 0.025)
    ),
    background_pct_upper = unname(
      quantile(background_pct, 0.975)
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    diagnosis_scope = .env$diagnosis_scope,
    period = factor(
      period,
      levels = c("overall", "winter", "other_weeks")
    )
  ) %>%
  arrange(period, age_index) %>%
  select(diagnosis_scope, everything())

# Export numeric estimates to the current scope-specific output directory.
export_table(
  four_pathogen_summary,
  "four_pathogen_attributable_fractions",
  analysis_output_dir
)

# Print tab-separated results for inspection and copying into Word.
period_titles <- c(
  overall = "ALL STUDY WEEKS",
  winter = "WINTER WEEKS",
  other_weeks = "OTHER WEEKS"
)

for (selected_period in names(period_titles)) {
  
  printed_table <- four_pathogen_summary %>%
    filter(as.character(period) == selected_period) %>%
    transmute(
      `Age group` = as.character(age_group),
      `Observed RTI admissions` = observed_total,
      `Four pathogens combined, % (95% CrI)` = sprintf(
        "%.1f (%.1f–%.1f)",
        four_pathogen_pct_median,
        four_pathogen_pct_lower,
        four_pathogen_pct_upper
      ),
      `Background, % (95% CrI)` = sprintf(
        "%.1f (%.1f–%.1f)",
        background_pct_median,
        background_pct_lower,
        background_pct_upper
      )
    )
  
  cat(
    "\n", period_titles[[selected_period]],
    " — ", diagnosis_scope, "\n\n",
    sep = ""
  )
  
  utils::write.table(
    printed_table,
    file = "",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

message(
  "Saved combined four-pathogen fractions in: ",
  analysis_output_dir
)

rm(
  pathogen_total_list,
  lambda_draw_matrix,
  combined_contribution
)

if (anyNA(attribution_draws$fraction_observed)) {
  warning(
    "Some period/age groups have zero observed admissions; ",
    "their percentage is undefined."
  )
}

# -------------------------------------------------------------------------
# Verify the partition in every posterior draw
# -------------------------------------------------------------------------

partition_check <- attribution_draws %>%
  select(
    draw,
    pathogen,
    age_group,
    period,
    attributed_total
  ) %>%
  pivot_wider(
    names_from = period,
    values_from = attributed_total
  )

stopifnot(
  all(
    abs(
      partition_check$overall -
        partition_check$winter -
        partition_check$other_weeks
    ) <
      1e-8 * pmax(1, partition_check$overall)
  )
)

# -------------------------------------------------------------------------
# Posterior summaries
# -------------------------------------------------------------------------

attribution_summary_functions <- list(
  median = function(x) {
    if (all(is.na(x))) return(NA_real_)
    median(x)
  },
  lower = function(x) {
    if (all(is.na(x))) return(NA_real_)
    unname(quantile(x, probs = 0.025))
  },
  upper = function(x) {
    if (all(is.na(x))) return(NA_real_)
    unname(quantile(x, probs = 0.975))
  }
)

attribution_summary <- attribution_draws %>%
  group_by(period, pathogen, age_group) %>%
  summarise(
    n_weeks = first(n_weeks),
    observed_total = first(observed_total),
    across(
      .cols = c(
        attributed_total,
        attributed_per_year,
        attributed_per_week,
        incidence_per_100000,
        fraction_observed
      ),
      .fns = attribution_summary_functions,
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    observed_per_week = observed_total / n_weeks,
    period = factor(period, levels = names(periods)),
    pathogen = factor(pathogen, levels = pathogen_table_order),
    age_group = factor(age_group, levels = age_groups)
  ) %>%
  arrange(period, pathogen, age_group)

# Confirm that summary denominators match the direct data calculation.
denominator_check <- attribution_summary %>%
  select(period, age_group, observed_total) %>%
  distinct() %>%
  left_join(
    attribution_denominators %>%
      select(
        period,
        age_group,
        observed_from_data = observed_total
      ),
    by = c("period", "age_group")
  )

stopifnot(
  !anyNA(denominator_check$observed_from_data),
  all(
    denominator_check$observed_total ==
      denominator_check$observed_from_data
  )
)

# Joint attribution: sum across pathogens within each draw first.
all_pathogens_draws <- attribution_draws %>%
  group_by(draw, period, age_group) %>%
  summarise(
    observed_total = first(observed_total),
    all_pathogens_total = sum(attributed_total),
    .groups = "drop"
  ) %>%
  mutate(
    all_pathogens_fraction = if_else(
      observed_total > 0,
      100 * all_pathogens_total / observed_total,
      NA_real_
    ),
    background_total = observed_total - all_pathogens_total
  )

stopifnot(
  all(
    all_pathogens_draws$background_total >=
      -1e-8 * pmax(1, all_pathogens_draws$observed_total)
  )
)

all_pathogens_summary <- all_pathogens_draws %>%
  group_by(period, age_group) %>%
  summarise(
    observed_total = first(observed_total),
    across(
      .cols = c(
        all_pathogens_total,
        all_pathogens_fraction,
        background_total
      ),
      .fns = attribution_summary_functions,
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    period = factor(period, levels = names(periods)),
    age_group = factor(age_group, levels = age_groups)
  ) %>%
  arrange(period, age_group)

# -------------------------------------------------------------------------
# Formatting helpers
# -------------------------------------------------------------------------

format_manuscript_interval <- function(
    estimate,
    lower,
    upper,
    percent = FALSE,
    digits = 1L) {
  
  number_format <- paste0("%.", digits, "f")
  
  interval_format <- paste0(
    number_format,
    if (percent) "%%" else "",
    " (",
    number_format,
    "\u2013",
    number_format,
    ")"
  )
  
  result <- sprintf(
    interval_format,
    estimate,
    lower,
    upper
  )
  
  result[
    is.na(estimate) | is.na(lower) | is.na(upper)
  ] <- "NA"
  
  result
}

table_labels <- function(summary_data) {
  summary_data %>%
    transmute(
      Pathogen = unname(
        manuscript_pathogen_labels[as.character(pathogen)]
      ),
      `Age group` = unname(
        manuscript_age_labels[as.character(age_group)]
      )
    )
}

print_tab_separated <- function(table_data) {
  write.table(
    as.data.frame(table_data),
    file = "",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

# -------------------------------------------------------------------------
# Incidence tables for all three periods
# -------------------------------------------------------------------------

incidence_column_labels <- c(
  "overall" =
    "Average annual admissions per 100,000 population (95% CrI)",
  "winter" =
    "Average annual winter-attributed admissions per 100,000 population (95% CrI)",
  "other_weeks" =
    "Average annual other-weeks-attributed admissions per 100,000 population (95% CrI)"
)

percentage_column_labels <- c(
  "overall" =
    "Percentage of RTI admissions explained (95% CrI)",
  "winter" =
    "Percentage of winter RTI admissions explained (95% CrI)",
  "other_weeks" =
    "Percentage of other-weeks RTI admissions explained (95% CrI)"
)

make_incidence_table <- function(selected_period) {
  
  summary_i <- attribution_summary %>%
    filter(as.character(period) == selected_period) %>%
    arrange(pathogen, age_group)
  
  estimates <- tibble(
    incidence = format_manuscript_interval(
      summary_i$incidence_per_100000_median,
      summary_i$incidence_per_100000_lower,
      summary_i$incidence_per_100000_upper
    ),
    percentage = format_manuscript_interval(
      summary_i$fraction_observed_median,
      summary_i$fraction_observed_lower,
      summary_i$fraction_observed_upper,
      percent = TRUE
    )
  )
  
  names(estimates) <- c(
    incidence_column_labels[[selected_period]],
    percentage_column_labels[[selected_period]]
  )
  
  bind_cols(
    table_labels(summary_i),
    estimates
  )
}

incidence_tables <- setNames(
  lapply(names(periods), make_incidence_table),
  names(periods)
)

overall_attribution_summary <- attribution_summary %>%
  filter(period == "overall")

attribution_table <- incidence_tables$overall
winter_attribution_incidence_table <- incidence_tables$winter
other_weeks_attribution_incidence_table <- incidence_tables$other_weeks

export_table(
  attribution_table,
  "attribution_table",
  analysis_output_dir,
  word = TRUE
)

export_table(
  winter_attribution_incidence_table,
  "winter_attribution_incidence_table",
  analysis_output_dir,
  word = TRUE
)

export_table(
  other_weeks_attribution_incidence_table,
  "other_weeks_attribution_incidence_table",
  analysis_output_dir,
  word = TRUE
)

# -------------------------------------------------------------------------
# Retain separate cumulative-count tables for winter and other weeks
# -------------------------------------------------------------------------

for (selected_period in c("winter", "other_weeks")) {
  
  period_data <- attribution_summary %>%
    filter(as.character(period) == selected_period)
  
  period_table <- bind_cols(
    table_labels(period_data),
    tibble(
      `Observed RTI admissions in selected weeks` =
        period_data$observed_total,
      
      `Virus-attributed admissions in selected weeks (95% CrI)` =
        format_manuscript_interval(
          period_data$attributed_total_median,
          period_data$attributed_total_lower,
          period_data$attributed_total_upper,
          digits = 0L
        ),
      
      `Percentage of RTI admissions in selected weeks explained (95% CrI)` =
        format_manuscript_interval(
          period_data$fraction_observed_median,
          period_data$fraction_observed_lower,
          period_data$fraction_observed_upper,
          percent = TRUE
        )
    )
  )
  
  export_table(
    period_table,
    paste0(selected_period, "_attribution_table"),
    analysis_output_dir,
    word = TRUE
  )
  
  export_table(
    period_data,
    paste0(selected_period, "_attribution_summary_numeric"),
    analysis_output_dir
  )
  
  if (selected_period == "winter") {
    winter_attribution_table <- period_table
  } else {
    other_weeks_attribution_table <- period_table
  }
}

# -------------------------------------------------------------------------
# Overall totals table
# -------------------------------------------------------------------------

totals_table <- bind_cols(
  table_labels(overall_attribution_summary),
  tibble(
    `Observed RTI admissions, 2013-2019` =
      overall_attribution_summary$observed_total,
    
    `Virus-attributed admissions, 2013-2019 (95% CrI)` =
      format_manuscript_interval(
        overall_attribution_summary$attributed_total_median,
        overall_attribution_summary$attributed_total_lower,
        overall_attribution_summary$attributed_total_upper,
        digits = 0L
      ),
    
    `Average virus-attributed admissions per year (95% CrI)` =
      format_manuscript_interval(
        overall_attribution_summary$attributed_per_year_median,
        overall_attribution_summary$attributed_per_year_lower,
        overall_attribution_summary$attributed_per_year_upper
      ),
    
    `Average attributed admissions per week` = round(
      overall_attribution_summary$attributed_per_week_median,
      2
    ),
    
    `Average observed admissions per week` = round(
      overall_attribution_summary$observed_per_week,
      2
    )
  )
)

export_table(
  totals_table,
  "attribution_totals_table",
  analysis_output_dir,
  word = TRUE
)

# -------------------------------------------------------------------------
# Numeric exports
# -------------------------------------------------------------------------

export_table(
  overall_attribution_summary,
  "attribution_summary_numeric",
  analysis_output_dir
)

export_table(
  attribution_summary,
  "attribution_all_periods_numeric",
  analysis_output_dir
)

export_table(
  all_pathogens_summary,
  "all_pathogens_summary",
  analysis_output_dir
)

export_table(
  population_denominators,
  "population_denominators",
  analysis_output_dir
)

export_table(
  attribution_denominators,
  "attribution_observed_denominators",
  analysis_output_dir
)

period_summary <- attribution_denominators %>%
  group_by(period) %>%
  summarise(
    n_weeks = first(n_weeks),
    observed_total = sum(observed_total),
    .groups = "drop"
  )

export_table(
  period_summary,
  "attribution_period_summary",
  analysis_output_dir
)

# Wide table of exact percentage denominators.
denominator_table <- attribution_denominators %>%
  select(age_group, period, observed_total) %>%
  pivot_wider(
    names_from = period,
    values_from = observed_total
  ) %>%
  arrange(age_group) %>%
  transmute(
    `Age group` = unname(
      manuscript_age_labels[as.character(age_group)]
    ),
    `Observed RTI admissions: all weeks` = overall,
    `Observed RTI admissions: winter` = winter,
    `Observed RTI admissions: other weeks` = other_weeks
  )

export_table(
  denominator_table,
  "attribution_denominators_table",
  analysis_output_dir,
  word = TRUE
)

# -------------------------------------------------------------------------
# Three-column manuscript incidence tables
# -------------------------------------------------------------------------

# Pathogen heading row followed by its age-group rows.
make_manuscript_attribution_table <- function(selected_period) {
  
  source_table <- incidence_tables[[selected_period]]
  
  result <- bind_rows(
    lapply(pathogen_table_order, function(selected_pathogen) {
      
      pathogen_label <- unname(
        manuscript_pathogen_labels[selected_pathogen]
      )
      
      pathogen_rows <- source_table %>%
        filter(Pathogen == pathogen_label) %>%
        select(-Pathogen)
      
      heading <- tibble(
        age = pathogen_label,
        incidence = "",
        percentage = ""
      )
      
      names(heading) <- names(pathogen_rows)
      
      bind_rows(heading, pathogen_rows)
    })
  )
  
  result
}

manuscript_attribution_tables <- setNames(
  lapply(names(periods), make_manuscript_attribution_table),
  names(periods)
)

# -------------------------------------------------------------------------
# Exact denominator text for manuscript footnotes
# -------------------------------------------------------------------------

winter_week_label <- paste0(
  winter_start_week,
  "\u2013",
  winter_end_week
)

period_descriptions <- c(
  "overall" = "all study weeks",
  "winter" = paste0(
    "winter weeks (weeks ",
    winter_week_label,
    ")"
  ),
  "other_weeks" = paste0(
    "weeks outside the winter definition (weeks ",
    winter_week_label,
    ")"
  )
)

make_denominator_text <- function(selected_period) {
  
  denominator_data <- attribution_denominators %>%
    filter(as.character(period) == selected_period) %>%
    arrange(age_group)
  
  entries <- paste0(
    format(
      denominator_data$observed_total,
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    ),
    " for ",
    unname(
      manuscript_age_labels[
        as.character(denominator_data$age_group)
      ]
    )
  )
  
  entry_text <- if (length(entries) == 1L) {
    entries
  } else {
    paste0(
      paste(entries[-length(entries)], collapse = ", "),
      ", and ",
      entries[length(entries)]
    )
  }
  
  paste0(
    "The observed numbers of RTI hospital admissions in ",
    period_descriptions[[selected_period]],
    " during 2013\u201319, used as denominators for the percentages, were ",
    entry_text,
    "."
  )
}

attribution_denominator_notes <- tibble(
  period = names(periods),
  denominator_text = vapply(
    names(periods),
    make_denominator_text,
    character(1)
  )
)

export_table(
  attribution_denominator_notes,
  "attribution_denominator_notes",
  analysis_output_dir
)

# -------------------------------------------------------------------------
# Printed output for inspection and copying into Word
# -------------------------------------------------------------------------

cat(
  "\nATTRIBUTION RESULTS FOR THE CURRENT FITTED MODEL\n",
  "Output directory: ", analysis_output_dir, "\n",
  "Estimates are posterior medians with 95% credible intervals.\n",
  "All incidence estimates use seven years of population person-time.\n",
  "Winter and other-weeks estimates are contributions to annual incidence.\n",
  sep = ""
)

cat("\nEXACT OBSERVED-ADMISSION DENOMINATORS\n\n")
print_tab_separated(denominator_table)

cat("\nPOPULATION DENOMINATORS USED FOR INCIDENCE\n\n")

print_tab_separated(
  population_denominators %>%
    transmute(
      `Age group` = unname(
        manuscript_age_labels[as.character(age_group)]
      ),
      `Mean population` = round(mean_population, 2),
      `Population person-years, 2013-2019` = round(person_years, 2)
    )
)

manuscript_period_titles <- c(
  "overall" = "ALL STUDY WEEKS",
  "winter" = "WINTER WEEKS",
  "other_weeks" = "OTHER WEEKS"
)

for (selected_period in names(periods)) {
  
  manuscript_table <-
    manuscript_attribution_tables[[selected_period]]
  
  cat(
    "\n\n",
    manuscript_period_titles[[selected_period]],
    " (",
    length(periods[[selected_period]]),
    " weeks)\n\n",
    sep = ""
  )
  
  cat(
    make_denominator_text(selected_period),
    "\n\n",
    sep = ""
  )
  
  # No row numbers, quotes, or tibble truncation.
  print_tab_separated(manuscript_table)
  
  export_table(
    manuscript_table,
    paste0(selected_period, "_attribution_manuscript"),
    analysis_output_dir,
    word = TRUE
  )
}

message(
  "\nSaved overall, winter, and other-weeks incidence tables, ",
  "separate count tables, numeric summaries, and exact denominators.\n",
  "The printed manuscript tables are tab-separated for copying into Word."
)