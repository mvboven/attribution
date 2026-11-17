# Gamma beta prior and posterior intervals.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
if (make_parameter_figures) {
  require_fit()
  # Prior shape/rate are plotting settings; they match Gamma(2,2) in the Stan file.
  beta_draws <- fit_hmc$draws(
    variables = "beta",
    format = "draws_df"
  ) %>%
    posterior::as_draws_df() %>%
    as_tibble() %>%
    select(
      starts_with("beta[")
    ) %>%
    mutate(
      draw = row_number()
    ) %>%
    pivot_longer(
      cols = -draw,
      names_to = "variable",
      values_to = "beta"
    ) %>%
    tidyr::extract(
      variable,
      into = c(
        "age_index",
        "pathogen_index"
      ),
      regex =
        "^beta\\[(\\d+),(\\d+)\\]$",
      convert = TRUE
    ) %>%
    mutate(
      age_group =
        age_groups[age_index],

      pathogen =
        names(pathogen_colors)[pathogen_index],

      age_group = factor(
        age_group,
        levels = rev(age_groups)
      ),

      pathogen = factor(
        pathogen,
        levels = names(pathogen_colors)
      )
    )

  stopifnot(
    !anyNA(beta_draws$beta),
    all(beta_draws$beta >= 0),
    n_distinct(beta_draws$age_group) == num_types,
    n_distinct(beta_draws$pathogen) == num_path
  )


  ### Summarise posterior distributions ###

  beta_posterior_summary <- beta_draws %>%
    group_by(
      pathogen,
      age_group
    ) %>%
    summarise(
      posterior_lower_95 =
        quantile(
          beta,
          0.025
        ),

      posterior_lower_50 =
        quantile(
          beta,
          0.25
        ),

      posterior_median =
        median(beta),

      posterior_upper_50 =
        quantile(
          beta,
          0.75
        ),

      posterior_upper_95 =
        quantile(
          beta,
          0.975
        ),

      .groups = "drop"
    )


  ### Prior intervals ###

  # Both prior intervals are one-sided:
  #
  # 0 to the 50th percentile
  # 0 to the 95th percentile

  prior_upper_50 <- qgamma(
    0.50,
    shape = beta_prior_shape,
    rate = beta_prior_rate
  )

  prior_upper_95 <- qgamma(
    0.95,
    shape = beta_prior_shape,
    rate = beta_prior_rate
  )


  ### Plotting limits ###

  # The lower plotting limit must be strictly positive because
  # beta is shown on a logarithmic axis.

  smallest_posterior_value <- min(
    beta_posterior_summary$posterior_lower_95[
      beta_posterior_summary$posterior_lower_95 > 0
    ],
    na.rm = TRUE
  )

  if (!is.finite(smallest_posterior_value) || smallest_posterior_value <= 0) {
    stop("No positive lower posterior beta interval is available for a logarithmic plot.")
  }
  beta_plot_limits <- c(
    smallest_posterior_value / 1.5,

    max(
      prior_upper_95,
      beta_posterior_summary$posterior_upper_95,
      na.rm = TRUE
    ) * 1.5
  )


  ### Construct prior summary ###

  beta_prior_summary <- tidyr::expand_grid(
    pathogen = factor(
      names(pathogen_colors),
      levels = names(pathogen_colors)
    ),

    age_group = factor(
      rev(age_groups),
      levels = rev(age_groups)
    )
  ) %>%
    mutate(
      # Zero cannot be displayed on a logarithmic axis.
      # Both prior intervals therefore start at the
      # left plotting boundary.
      prior_lower_display =
        beta_plot_limits[1],

      prior_upper_50 =
        prior_upper_50,

      prior_upper_95 =
        prior_upper_95
    )


  ### Axis breaks ###

  beta_axis_breaks <- 10^(-8:3)

  beta_axis_breaks <- beta_axis_breaks[
    beta_axis_breaks >= beta_plot_limits[1] &
      beta_axis_breaks <= beta_plot_limits[2]
  ]


  ### Make figure ###

  plot_beta_prior_posterior <- ggplot() +

    # One-sided 95% prior range
    geom_segment(
      data = beta_prior_summary,
      aes(
        x = prior_lower_display,
        xend = prior_upper_95,
        y = age_group,
        yend = age_group
      ),
      color = "grey82",
      linewidth = 4.5,
      lineend = "round"
    ) +

    # One-sided 50% prior range
    geom_segment(
      data = beta_prior_summary,
      aes(
        x = prior_lower_display,
        xend = prior_upper_50,
        y = age_group,
        yend = age_group
      ),
      color = "grey62",
      linewidth = 4.5,
      lineend = "round"
    ) +

    # Prior median
    geom_point(
      data = beta_prior_summary,
      aes(
        x = prior_upper_50,
        y = age_group
      ),
      shape = 21,
      fill = "white",
      color = "grey35",
      size = 2
    ) +

    # Posterior 95% credible interval
    geom_segment(
      data = beta_posterior_summary,
      aes(
        x = posterior_lower_95,
        xend = posterior_upper_95,
        y = age_group,
        yend = age_group,
        color = pathogen
      ),
      linewidth = 0.8,
      lineend = "round"
    ) +

    # Posterior 50% credible interval
    geom_segment(
      data = beta_posterior_summary,
      aes(
        x = posterior_lower_50,
        xend = posterior_upper_50,
        y = age_group,
        yend = age_group,
        color = pathogen
      ),
      linewidth = 2.5,
      lineend = "round"
    ) +

    # Posterior median
    geom_point(
      data = beta_posterior_summary,
      aes(
        x = posterior_median,
        y = age_group,
        color = pathogen
      ),
      size = 2.2
    ) +

    facet_wrap(
      ~pathogen,
      ncol = 2
    ) +

    scale_color_manual(
      values = pathogen_colors,
      guide = "none"
    ) +

    scale_x_log10(
      limits = beta_plot_limits,
      breaks = beta_axis_breaks,
      labels = format_beta_axis
    ) +

    labs(
      title =
        "Prior and posterior distributions of pathogen coefficients",

      x =
        expression(beta ~ "(logarithmic scale)"),

      y =
        "Age group"
    ) +

    theme_minimal() +

    theme(
      panel.grid.minor =
        element_blank(),

      panel.grid.major.y =
        element_blank(),

      plot.title =
        element_text(
          size = 12
        ),

      strip.text =
        element_text(
          size = 11
        )
    )

  if (show_figures) print(plot_beta_prior_posterior)


  ### Save figure ###

  ggsave(
    filename = file.path(
      analysis_figure_dir,
      "beta_prior_posterior_intervals.png"
    ),
    plot =
      plot_beta_prior_posterior,
    dpi = 300,
    height = 7,
    width = 10,
    bg = "white"
  )

  export_table(beta_posterior_summary, "beta_posterior_intervals_numeric", analysis_output_dir)
  message("Saved the beta prior/posterior figure (Gamma prior shape ", beta_prior_shape,
          ", rate ", beta_prior_rate, ").")
} else message("Beta prior/posterior figure skipped.")
