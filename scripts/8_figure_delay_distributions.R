# Posterior mean discrete delay distributions and 95% credible intervals.
# Lancet-inspired presentation.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.

if (make_parameter_figures) {
  
  require_fit()
  
  require_objects(
    c(
      "min_lag",
      "max_lag",
      "age_groups",
      "pathogen_colors",
      "analysis_figure_dir",
      "analysis_output_dir",
      "show_figures"
    ),
    instruction = "Run the data preparation scripts first."
  )
  
  delay_figure_font <- "sans"
  delay_figure_width_mm <- 180
  delay_figure_height_mm <- 150
  
  lag_values <- min_lag:max_lag
  
  # Preserve the existing mapping from Stan indices to pathogens.
  # This order must match the order used to prepare the Stan data.
  delay_stan_pathogen_order <- names(pathogen_colors)
  
  # Display order is specified separately.
  delay_facet_labels <- c(
    "Influenza" = "A  Influenza",
    "RSV" = "B  RSV",
    "HMPV" = "C  HMPV",
    "Parainfluenza" = "D  PIV"
  )
  
  delay_component_order <- names(delay_facet_labels)
  
  age_alpha <- c(
    "0-2m"   = 0.25,
    "3-5m"   = 0.35,
    "6-11m"  = 0.45,
    "1-2y"   = 0.55,
    "3-4y"   = 0.65,
    "65-74y" = 0.75,
    "75-84y" = 0.87,
    "85+y"   = 1.00
  )
  
  stopifnot(
    all(delay_component_order %in% delay_stan_pathogen_order),
    all(age_groups %in% names(age_alpha)),
    !anyDuplicated(delay_stan_pathogen_order),
    !anyDuplicated(age_groups)
  )
  
  # Extract posterior draws of the discrete lag probabilities.
  # Stan dimensions: age group x pathogen x lag.
  delay_draws <- fit_hmc$draws(
    variables = "pmf",
    format = "draws_df"
  ) %>%
    posterior::as_draws_df() %>%
    as_tibble() %>%
    select(starts_with("pmf[")) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(
      cols = -draw,
      names_to = "variable",
      values_to = "probability"
    ) %>%
    tidyr::extract(
      variable,
      into = c(
        "age_index",
        "pathogen_index",
        "lag_index"
      ),
      regex = "^pmf\\[(\\d+),(\\d+),(\\d+)\\]$",
      convert = TRUE
    ) %>%
    mutate(
      age_group = factor(
        age_groups[age_index],
        levels = age_groups
      ),
      pathogen = factor(
        delay_stan_pathogen_order[pathogen_index],
        levels = delay_component_order
      ),
      lag = lag_values[lag_index]
    )
  
  stopifnot(
    !anyNA(delay_draws$age_group),
    !anyNA(delay_draws$pathogen),
    !anyNA(delay_draws$lag),
    all(is.finite(delay_draws$probability))
  )
  
  # Posterior means sum to one across lags for each age/pathogen pair.
  # Credible intervals are marginal intervals for each lag probability.
  delay_summary <- delay_draws %>%
    group_by(pathogen, age_group, lag) %>%
    summarise(
      mean = mean(probability),
      lower = quantile(probability, 0.025),
      upper = quantile(probability, 0.975),
      .groups = "drop"
    ) %>%
    arrange(pathogen, age_group, lag)
  
  delay_probability_sums <- delay_summary %>%
    group_by(pathogen, age_group) %>%
    summarise(
      total = sum(mean),
      .groups = "drop"
    )
  
  stopifnot(
    all(abs(delay_probability_sums$total - 1) < 1e-8),
    all(is.finite(delay_summary$mean)),
    all(is.finite(delay_summary$lower)),
    all(is.finite(delay_summary$upper)),
    all(delay_summary$lower >= 0),
    all(delay_summary$upper <= 1 + 1e-8)
  )
  
  delay_dodge <- position_dodge(width = 0.85)
  
  plot_delay <- ggplot(
    delay_summary,
    aes(
      x = factor(lag, levels = lag_values),
      y = mean,
      group = age_group
    )
  ) +
    
    # Age is represented by increasing colour opacity.
    geom_col(
      aes(
        fill = pathogen,
        alpha = age_group
      ),
      position = delay_dodge,
      width = 0.78,
      colour = NA
    ) +
    
    # Keep intervals fully opaque for readability.
    geom_errorbar(
      aes(
        ymin = lower,
        ymax = upper,
        colour = pathogen
      ),
      position = delay_dodge,
      width = 0.10,
      linewidth = 0.25,
      alpha = 1,
      show.legend = FALSE
    ) +
    
    facet_wrap(
      ~pathogen,
      ncol = 2,
      scales = "fixed",
      labeller = as_labeller(delay_facet_labels)
    ) +
    
    scale_fill_manual(
      values = pathogen_colors,
      guide = "none"
    ) +
    
    scale_colour_manual(
      values = pathogen_colors,
      guide = "none"
    ) +
    
    scale_alpha_manual(
      values = age_alpha,
      breaks = age_groups,
      limits = age_groups,
      drop = FALSE,
      name = "Age group",
      guide = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          fill = "grey20",
          colour = NA
        )
      )
    ) +
    
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(add = 0.55)
    ) +
    
    # Common probability scale for comparison across pathogens.
    scale_y_continuous(
      limits = c(0, 1 + 1e-8),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::label_number(accuracy = 0.01),
      expand = expansion(mult = c(0, 0.03))
    ) +
    
    labs(
      x = "Lag (weeks)",
      y = "Probability"
    ) +
    
    theme_minimal(
      base_size = 9,
      base_family = delay_figure_font
    ) +
    
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(
        colour = "grey92",
        linewidth = 0.2
      ),
      
      axis.line = element_line(
        colour = "grey40",
        linewidth = 0.25
      ),
      axis.ticks = element_line(
        colour = "grey40",
        linewidth = 0.25
      ),
      axis.ticks.length = grid::unit(1, "mm"),
      
      axis.text = element_text(
        size = 8,
        colour = "grey20"
      ),
      axis.title = element_text(
        size = 9,
        colour = "grey20"
      ),
      axis.title.x = element_text(
        margin = margin(t = 6)
      ),
      axis.title.y = element_text(
        margin = margin(r = 6)
      ),
      
      strip.background = element_blank(),
      strip.text = element_text(
        size = 10,
        face = "bold",
        colour = "grey15",
        hjust = 0,
        margin = margin(t = 3, b = 6)
      ),
      
      panel.spacing = grid::unit(5, "mm"),
      plot.margin = margin(
        t = 5,
        r = 4,
        b = 3,
        l = 4
      ),
      
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(
        size = 9,
        colour = "grey20"
      ),
      legend.text = element_text(
        size = 8,
        colour = "grey20"
      ),
      legend.key.width = grid::unit(5, "mm"),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.spacing.x = grid::unit(2, "mm"),
      legend.margin = margin(
        t = 4,
        r = 0,
        b = 0,
        l = 0
      )
    )
  
  if (show_figures) {
    print(plot_delay)
  }
  
  dir.create(
    analysis_figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Vector PDF.
  ggsave(
    filename = file.path(
      analysis_figure_dir,
      "delay_distributions.pdf"
    ),
    plot = plot_delay,
    device = grDevices::cairo_pdf,
    width = delay_figure_width_mm,
    height = delay_figure_height_mm,
    units = "mm",
    bg = "white"
  )
  
  # High-resolution PNG at the same dimensions.
  ggsave(
    filename = file.path(
      analysis_figure_dir,
      "delay_distributions.png"
    ),
    plot = plot_delay,
    dpi = 600,
    width = delay_figure_width_mm,
    height = delay_figure_height_mm,
    units = "mm",
    bg = "white"
  )
  
  export_table(
    delay_summary,
    "delay_distributions_numeric",
    analysis_output_dir
  )
  
  message(
    "Saved delay-distribution figure (PDF and PNG) and numeric summary."
  )
  
} else {
  message("Delay-distribution figure skipped.")
}