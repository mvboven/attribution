# Linked-admissions correlation figure: Lancet-inspired presentation.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
# Use as the replacement for scripts/14_validation_figure_correlation.R.
# Run after 13_prepare_validation_data.R; no refitting is required.

if (run_validation) {
  require_objects(
    c(
      "validation_weekly",
      "validation_statistics",
      "validation_run_signature",
      "run_signature",
      "pathogen_colors"
    ),
    instruction = "Source 13_prepare_validation_data.R first."
  )
  
  if (!identical(validation_run_signature, run_signature)) {
    stop("Validation objects do not match the current model/data.")
  }
  
  # Match the time-series figure's typography and export width.
  # These are practical presentation settings, not journal specifications.
  validation_figure_font <- "sans"
  validation_figure_width_mm <- 180
  validation_figure_height_mm <- 140
  
  validation_facet_labels <- c(
    "Influenza" = "A  Influenza",
    "RSV" = "B  RSV",
    "HMPV" = "C  HMPV",
    "Parainfluenza" = "D  Parainfluenza"
  )
  validation_component_order <- names(validation_facet_labels)
  
  stopifnot(
    all(validation_component_order %in% names(pathogen_colors)),
    all(as.character(validation_weekly$component) %in%
          validation_component_order),
    all(as.character(validation_statistics$component) %in%
          validation_component_order)
  )
  
  validation_plot_data <- validation_weekly %>%
    mutate(
      component = factor(component, levels = validation_component_order)
    )
  
  validation_labels <- validation_statistics %>%
    mutate(
      component = factor(component, levels = validation_component_order),
      label = if_else(
        is.finite(spearman_rho) &
          is.finite(rho_ci_low) &
          is.finite(rho_ci_high),
        sprintf(
          "Spearman \u03c1 = %.2f\n95%% block CI %.2f\u2013%.2f",
          spearman_rho,
          rho_ci_low,
          rho_ci_high
        ),
        "Spearman \u03c1 / 95% block CI\nundefined"
      )
    )
  
  plot_validation_correlation <- ggplot(
    validation_plot_data,
    aes(
      x = confirmed,
      y = attributed,
      colour = component
    )
  ) +
    geom_point(alpha = 0.35, size = 0.7, shape = 16) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      span = 0.8,
      linewidth = 0.55
    ) +
    geom_text(
      data = validation_labels,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.03,
      vjust = 1.15,
      size = 2.8,  # Approximately 8 points.
      lineheight = 1.1,
      family = validation_figure_font,
      colour = "grey20"
    ) +
    facet_wrap(
      ~component,
      scales = "free",
      ncol = 2,
      labeller = as_labeller(validation_facet_labels)
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    # Reserve space above the observations for the two-line annotation.
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.32))
    ) +
    scale_colour_manual(values = pathogen_colors, guide = "none") +
    labs(
      x = "Confirmed pathogen-linked admissions per week",
      # Stage 13 sums age-specific median fitted pathogen intensities.
      y = "Estimated pathogen-attributed admissions per week"
    ) +
    theme_minimal(
      base_size = 9,
      base_family = validation_figure_font
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(
        colour = "grey92",
        linewidth = 0.2
      ),
      axis.line = element_line(colour = "grey50", linewidth = 0.2),
      axis.ticks = element_line(colour = "grey50", linewidth = 0.2),
      axis.ticks.length = grid::unit(1, "mm"),
      axis.text = element_text(size = 8, colour = "grey20"),
      axis.title = element_text(size = 9, colour = "grey20"),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6)),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 10,
        colour = "grey20",
        hjust = 0,
        margin = margin(t = 3, b = 5)
      ),
      panel.spacing = grid::unit(4, "mm"),
      plot.margin = margin(t = 5, r = 4, b = 3, l = 4),
      legend.position = "none"
    )
  
  if (show_figures) print(plot_validation_correlation)
  
  # Vector PDF: Cairo preserves Unicode labels and embeds the resolved font.
  ggsave(
    filename = file.path(
      analysis_figure_dir,
      "validation_confirmed_pathogens_correlation.pdf"
    ),
    plot = plot_validation_correlation,
    device = grDevices::cairo_pdf,
    width = validation_figure_width_mm,
    height = validation_figure_height_mm,
    units = "mm",
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      analysis_figure_dir,
      "validation_confirmed_pathogens_correlation.png"
    ),
    plot = plot_validation_correlation,
    dpi = 600,
    width = validation_figure_width_mm,
    height = validation_figure_height_mm,
    units = "mm",
    bg = "white"
  )
  
  message("Saved linked-admissions correlation figure (PDF and PNG).")
} else {
  message("Validation figure skipped.")
}

