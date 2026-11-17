# Confirmed and scaled fitted pathogen time series.
# Lancet-inspired presentation; practical settings, not journal specifications.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
# Replacement for scripts/15_validation_figure_time_series.R.
# Run after 13_prepare_validation_data.R; no refitting is required.

if (run_validation) {
  
  require_objects(
    c(
      "validation_weekly",
      "validation_statistics",
      "validation_run_signature",
      "run_signature",
      "week_seq",
      "pathogen_colors",
      "analysis_figure_dir",
      "show_figures"
    ),
    instruction = "Source 13_prepare_validation_data.R first."
  )
  
  if (!identical(validation_run_signature, run_signature)) {
    stop("Validation objects do not match the current model/data.")
  }
  
  # Figure dimensions and typography.
  validation_figure_font <- "sans"
  validation_figure_width_mm <- 180
  validation_figure_height_mm <- 140
  
  # Panel order and labels; scaling factors are not included.
  validation_facet_labels <- c(
    "Influenza" = "A  Influenza",
    "RSV" = "B  RSV",
    "HMPV" = "C  HMPV",
    "Parainfluenza" = "D  PIV"
  )
  
  validation_component_order <- names(validation_facet_labels)
  
  stopifnot(
    all(validation_component_order %in% names(pathogen_colors)),
    setequal(
      as.character(validation_weekly$component),
      validation_component_order
    ),
    setequal(
      as.character(validation_statistics$component),
      validation_component_order
    ),
    !anyDuplicated(as.character(validation_statistics$component)),
    all(is.finite(validation_statistics$scaling_factor)),
    all(validation_statistics$scaling_factor > 0)
  )
  
  # Stage 13 supplies the multiplier from attributed to confirmed.
  # Its reciprocal gives the multiplier from confirmed to attributed.
  validation_scaling_output <- validation_statistics %>%
    transmute(
      component = factor(
        as.character(component),
        levels = validation_component_order
      ),
      confirmed_to_attributed = 1 / scaling_factor
    ) %>%
    arrange(component)
  
  stopifnot(
    all(is.finite(validation_scaling_output$confirmed_to_attributed))
  )
  
  cat(
    "\nConfirmed-to-attributed scaling factors\n",
    "Multiply confirmed admissions by these factors to match ",
    "the attributed totals over the validation period.\n",
    "These are overall multipliers, not week-specific ratios.\n\n",
    sep = ""
  )
  
  print(
    validation_scaling_output %>%
      transmute(
        Pathogen = ifelse(
          as.character(component) == "Parainfluenza",
          "PIV",
          as.character(component)
        ),
        `Confirmed-to-attributed multiplier` = sprintf(
          "%.2f",
          confirmed_to_attributed
        )
      ),
    n = Inf,
    width = Inf
  )
  
  # Scale fitted attribution onto the confirmed-count scale.
  validation_time_plot_data <- validation_weekly %>%
    mutate(component = as.character(component)) %>%
    left_join(
      validation_statistics %>%
        transmute(
          component = as.character(component),
          scaling_factor
        ),
      by = "component"
    ) %>%
    mutate(
      week = as.Date(week),
      component = factor(
        component,
        levels = validation_component_order
      ),
      scaled_attributed = scaling_factor * attributed
    ) %>%
    arrange(component, week)
  
  # Use the study dates and show a year label every two years.
  validation_date_limits <- range(as.Date(week_seq), na.rm = TRUE)
  
  validation_first_year <- as.integer(
    format(validation_date_limits[1], "%Y")
  )
  validation_last_year <- as.integer(
    format(validation_date_limits[2], "%Y")
  )
  
  validation_year_breaks <- as.Date(
    paste0(
      seq.int(validation_first_year, validation_last_year),
      "-01-01"
    )
  )
  
  validation_year_breaks <- validation_year_breaks[
    validation_year_breaks >= validation_date_limits[1] &
      validation_year_breaks <= validation_date_limits[2]
  ]
  
  validation_year_breaks <- validation_year_breaks[
    seq_along(validation_year_breaks) %% 2L == 1L
  ]
  
  plot_validation_time <- ggplot() +
    
    # Fine dark-grey reference curve.
    geom_line(
      data = validation_time_plot_data,
      aes(
        x = week,
        y = confirmed,
        group = component,
        linetype = "Confirmed"
      ),
      colour = "grey25",
      linewidth = 0.3,
      alpha = 0.85
    ) +
    
    # Coloured fitted curve drawn on top.
    geom_line(
      data = validation_time_plot_data,
      aes(
        x = week,
        y = scaled_attributed,
        group = component,
        colour = component,
        linetype = "Scaled fitted attribution"
      ),
      linewidth = 0.55,
      alpha = 1
    ) +
    
    facet_wrap(
      ~component,
      scales = "free_y",
      ncol = 2,
      labeller = as_labeller(validation_facet_labels)
    ) +
    
    scale_colour_manual(
      values = pathogen_colors,
      guide = "none"
    ) +
    
    scale_linetype_manual(
      values = c(
        "Confirmed" = "solid",
        "Scaled fitted attribution" = "solid"
      ),
      breaks = c(
        "Confirmed",
        "Scaled fitted attribution"
      ),
      name = NULL,
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          colour = c("grey25", "grey40"),
          linewidth = c(0.3, 0.55),
          alpha = c(0.85, 1)
        )
      )
    ) +
    
    scale_x_date(
      breaks = validation_year_breaks,
      date_labels = "%Y",
      limits = validation_date_limits,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    
    scale_y_continuous(
      limits = function(x) c(0, max(x, na.rm = TRUE)),
      labels = scales::label_number(
        accuracy = 1,
        big.mark = ","
      ),
      expand = expansion(mult = c(0, 0.06))
    ) +
    
    labs(
      x = "Year",
      y = "Pathogen-linked admissions per week"
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
      legend.spacing.x = grid::unit(3, "mm"),
      legend.text = element_text(
        size = 9,
        colour = "grey20"
      ),
      legend.key.width = grid::unit(8, "mm"),
      legend.key.height = grid::unit(4, "mm"),
      legend.margin = margin(
        t = 4,
        r = 0,
        b = 0,
        l = 0
      )
    )
  
  if (show_figures) {
    print(plot_validation_time)
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
      "validation_confirmed_pathogens_time.pdf"
    ),
    plot = plot_validation_time,
    device = grDevices::cairo_pdf,
    width = validation_figure_width_mm,
    height = validation_figure_height_mm,
    units = "mm",
    bg = "white"
  )
  
  # High-resolution PNG with identical dimensions.
  ggsave(
    filename = file.path(
      analysis_figure_dir,
      "validation_confirmed_pathogens_time.png"
    ),
    plot = plot_validation_time,
    dpi = 600,
    width = validation_figure_width_mm,
    height = validation_figure_height_mm,
    units = "mm",
    bg = "white"
  )
  
  message(
    "Saved confirmed and scaled fitted pathogen time series (PDF and PNG)."
  )
  
} else {
  message("Validation figure skipped.")
}