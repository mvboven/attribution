# Posterior seasonal multipliers: Lancet-inspired presentation.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.

if (make_parameter_figures) {
  
  require_fit()
  
  require_objects(
    c(
      "pathogen_colors",
      "season_labels",
      "round_log_down",
      "round_log_up",
      "format_multiplier",
      "analysis_figure_dir",
      "analysis_output_dir",
      "show_figures"
    ),
    instruction = "Run the preparation scripts first."
  )
  
  season_figure_font <- "sans"
  season_figure_width_mm <- 180
  season_figure_height_mm <- 145
  
  # Preserve the original mapping from Stan indices to pathogens.
  # This must match the pathogen order used to prepare the Stan data.
  season_stan_pathogen_order <- names(pathogen_colors)
  
  # Display order is specified separately.
  season_panel_labels <- c(
    "Influenza" = "A  Influenza",
    "RSV" = "B  RSV",
    "HMPV" = "C  HMPV",
    "Parainfluenza" = "D  PIV"
  )
  
  season_component_order <- names(season_panel_labels)
  
  stopifnot(
    all(season_component_order %in% season_stan_pathogen_order),
    !anyDuplicated(season_stan_pathogen_order),
    !anyDuplicated(season_labels)
  )
  
  # Extract draws of the normalised seasonal multipliers.
  # Stan dimensions: pathogen x season.
  season_draws <- fit_hmc$draws(
    variables = "seasoneffect",
    format = "draws_df"
  ) %>%
    posterior::as_draws_df() %>%
    as_tibble() %>%
    select(starts_with("seasoneffect[")) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(
      cols = -draw,
      names_to = "variable",
      values_to = "seasonal_multiplier"
    ) %>%
    tidyr::extract(
      variable,
      into = c(
        "pathogen_index",
        "season_index"
      ),
      regex = "^seasoneffect\\[(\\d+),(\\d+)\\]$",
      convert = TRUE
    ) %>%
    mutate(
      pathogen = factor(
        season_stan_pathogen_order[pathogen_index],
        levels = season_component_order
      ),
      season = factor(
        season_labels[season_index],
        levels = season_labels
      )
    )
  
  stopifnot(
    !anyNA(season_draws$pathogen),
    !anyNA(season_draws$season),
    all(is.finite(season_draws$seasonal_multiplier)),
    all(season_draws$seasonal_multiplier > 0)
  )
  
  # Posterior medians; retain the existing export structure.
  season_medians <- season_draws %>%
    group_by(pathogen, season) %>%
    summarise(
      median = median(seasonal_multiplier),
      .groups = "drop"
    ) %>%
    arrange(pathogen, season)
  
  # Common visual settings for all four panels.
  season_theme <- theme_minimal(
    base_size = 9,
    base_family = season_figure_font
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
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
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
      
      plot.title = element_text(
        size = 10,
        face = "bold",
        colour = "grey15",
        hjust = 0,
        margin = margin(b = 7)
      ),
      plot.title.position = "plot",
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 4,
        l = 4
      ),
      
      legend.position = "none"
    )
  
  make_season_panel <- function(
    selected_pathogen,
    show_x_title = TRUE,
    show_y_title = TRUE) {
    
    data_i <- season_draws %>%
      filter(pathogen == selected_pathogen)
    
    medians_i <- season_medians %>%
      filter(pathogen == selected_pathogen)
    
    panel_colour <- unname(
      pathogen_colors[selected_pathogen]
    )
    
    # Include the full posterior range and the reference value of one.
    y_min <- min(data_i$seasonal_multiplier, 1)
    y_max <- max(data_i$seasonal_multiplier, 1)
    
    y_min_plot <- round_log_down(y_min)
    y_max_plot <- round_log_up(y_max)
    
    # Guard against rounding functions returning an inward limit.
    y_min_plot <- min(y_min_plot, y_min)
    y_max_plot <- max(y_max_plot, y_max)
    
    stopifnot(
      is.finite(y_min_plot),
      is.finite(y_max_plot),
      y_min_plot > 0
    )
    
    if (y_min_plot == y_max_plot) {
      y_min_plot <- y_min_plot / 1.1
      y_max_plot <- y_max_plot * 1.1
    }
    
    y_breaks <- sort(
      unique(c(y_min_plot, 1, y_max_plot))
    )
    
    ggplot(
      data_i,
      aes(
        x = season,
        y = seasonal_multiplier
      )
    ) +
      
      # Reference line behind the posterior distributions.
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        linewidth = 0.3,
        colour = "grey50"
      ) +
      
      geom_violin(
        fill = panel_colour,
        colour = panel_colour,
        scale = "width",
        trim = TRUE,
        width = 0.8,
        linewidth = 0.35,
        alpha = 0.25
      ) +
      
      geom_point(
        data = medians_i,
        aes(
          x = season,
          y = median
        ),
        inherit.aes = FALSE,
        colour = panel_colour,
        size = 1.4
      ) +
      
      scale_x_discrete(
        drop = FALSE,
        expand = expansion(add = 0.55)
      ) +
      
      scale_y_log10(
        breaks = y_breaks,
        labels = format_multiplier,
        expand = expansion(mult = c(0.04, 0.06))
      ) +
      
      # Apply display limits without dropping draws before density estimation.
      coord_cartesian(
        ylim = c(y_min_plot, y_max_plot)
      ) +
      
      labs(
        title = unname(
          season_panel_labels[selected_pathogen]
        ),
        x = if (show_x_title) "Season" else NULL,
        y = if (show_y_title) "Seasonal multiplier" else NULL
      ) +
      
      season_theme
  }
  
  p_influenza <- make_season_panel(
    "Influenza",
    show_x_title = FALSE,
    show_y_title = TRUE
  )
  
  p_rsv <- make_season_panel(
    "RSV",
    show_x_title = FALSE,
    show_y_title = FALSE
  )
  
  p_hmpv <- make_season_panel(
    "HMPV",
    show_x_title = TRUE,
    show_y_title = TRUE
  )
  
  p_piv <- make_season_panel(
    "Parainfluenza",
    show_x_title = TRUE,
    show_y_title = FALSE
  )
  
  plot_season <- cowplot::plot_grid(
    p_influenza,
    p_rsv,
    p_hmpv,
    p_piv,
    ncol = 2,
    align = "hv",
    axis = "tblr"
  )
  
  if (show_figures) {
    print(plot_season)
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
      "seasonal_multipliers.pdf"
    ),
    plot = plot_season,
    device = grDevices::cairo_pdf,
    width = season_figure_width_mm,
    height = season_figure_height_mm,
    units = "mm",
    bg = "white"
  )
  
  # High-resolution PNG at the same dimensions.
  ggsave(
    filename = file.path(
      analysis_figure_dir,
      "seasonal_multipliers.png"
    ),
    plot = plot_season,
    dpi = 600,
    width = season_figure_width_mm,
    height = season_figure_height_mm,
    units = "mm",
    bg = "white"
  )
  
  export_table(
    season_medians,
    "seasonal_multiplier_medians",
    analysis_output_dir
  )
  
  message(
    "Saved seasonal-multiplier figure (PDF and PNG) and numeric medians."
  )
  
} else {
  message("Seasonal-multiplier figure skipped.")
}