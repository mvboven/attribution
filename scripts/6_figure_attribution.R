# Observed admissions, fitted totals, and pathogen contributions.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.

if (make_hmc_figure) {
  require_fit()
  
  # Use the actual dates of the ISO weeks.
  week_seq <- as.Date(data$week_start)
  
  stopifnot(
    length(week_seq) == n_data,
    nrow(y_data) == n_data,
    all(diff(week_seq) == 7),
    "2015-W53" %in% data$ISOweek
  )
  
  # Overall fitted incidence: posterior median only.
  overall_summary <- fit_hmc$summary(
    variables = "lambda",
    "median"
  ) %>%
    tidyr::extract(
      variable,
      into = c("age_index", "week_index"),
      regex = "^lambda\\[(\\d+),(\\d+)\\]$",
      convert = TRUE
    ) %>%
    transmute(
      week = week_seq[week_index],
      y_data = median,
      y_lo = NA_real_,
      y_hi = NA_real_,
      age_group = age_groups[age_index],
      component = "Overall"
    )
  
  # Shared summary also used by validation; no figure dependency.
  pathogen_summary <- summarise_pathogen_intensities(
    fit_hmc, week_seq, age_groups, pathogen_lookup$pathogen
  )
  
  # Observed admissions.
  observed_data <- tidyr::expand_grid(
    age_index = seq_len(num_types),
    week_index = seq_len(n_data)
  ) %>%
    transmute(
      week = week_seq[week_index],
      y_data = y_data[cbind(week_index, age_index)],
      y_lo = NA_real_,
      y_hi = NA_real_,
      age_group = age_groups[age_index],
      component = "Data"
    )
  
  # Combine plotting data. Keep internal component and age-group names.
  plot_data_hmc <- bind_rows(
    observed_data,
    overall_summary,
    pathogen_summary
  ) %>%
    mutate(
      age_group = factor(age_group, levels = age_groups),
      component = factor(
        component,
        levels = c("Data", "Overall", names(pathogen_colors))
      )
    )
  
  # Consistency checks.
  stopifnot(
    nrow(overall_summary) == num_types * n_data,
    nrow(pathogen_summary) == num_path * num_types * n_data,
    !anyNA(pathogen_summary$y_data),
    !anyNA(pathogen_summary$y_lo),
    !anyNA(pathogen_summary$y_hi)
  )
  
  # Separate data used by the different plot layers.
  plot_data_observed <- plot_data_hmc %>%
    filter(component == "Data")
  
  plot_data_overall <- plot_data_hmc %>%
    filter(component == "Overall")
  
  plot_data_pathogens <- plot_data_hmc %>%
    filter(component %in% names(pathogen_colors))
  
  # Use the shared pathogen palette for both lines and ribbons.
  # Change pathogen_colors upstream to update the pathogen colours.
  # Add the light-grey Total locally; cbPalette is not used here.
  # "Overall" remains the data value; its legend label becomes "Total".
  figure_colours <- c(Overall = "grey70", pathogen_colors)
  legend_breaks <- c("Overall", names(pathogen_colors))
  legend_labels <- c("Total", names(pathogen_colors))
  
  # Figure settings. Font sizes are points at the exported dimensions.
  figure_font <- "sans"
  figure_width_mm <- 180
  figure_height_mm <- 125
  
  # Display labels only: preserve the internal age-group values.
  age_labels <- c(
    "0-2m" = "0\u20132 months",
    "3-5m" = "3\u20135 months",
    "6-11m" = "6\u201311 months",
    "1-2y" = "1\u20132 years",
    "3-4y" = "3\u20134 years",
    "65-74y" = "65\u201374 years",
    "75-84y" = "75\u201384 years",
    "85+y" = "\u226585 years"
  )
  
  # Helper function for one age-group panel.
  make_age_panel <- function(selected_age_group,
                             show_legend = FALSE,
                             show_x_title = FALSE,
                             show_y_title = FALSE) {
    
    data_observed_i <- plot_data_observed %>%
      filter(age_group == selected_age_group)
    
    data_overall_i <- plot_data_overall %>%
      filter(age_group == selected_age_group)
    
    data_pathogens_i <- plot_data_pathogens %>%
      filter(age_group == selected_age_group)
    
    panel_title <- if (selected_age_group %in% names(age_labels)) {
      unname(age_labels[[selected_age_group]])
    } else {
      selected_age_group
    }
    
    ggplot() +
      # Draw uncertainty behind the fitted curves.
      geom_ribbon(
        data = data_pathogens_i,
        aes(
          x = week,
          ymin = y_lo,
          ymax = y_hi,
          fill = component,
          group = component
        ),
        alpha = 0.12,
        colour = NA,
        show.legend = FALSE
      ) +
      # Fitted total: light grey and slightly thinner.
      geom_line(
        data = data_overall_i,
        aes(
          x = week,
          y = y_data,
          group = component,
          colour = component
        ),
        linewidth = 0.4,
        alpha = 1
      ) +
      # Observed admissions.
      geom_point(
        data = data_observed_i,
        aes(x = week, y = y_data),
        colour = "grey20",
        alpha = 0.6,
        size = 0.3,
        show.legend = FALSE
      ) +
      # Draw pathogen curves last, using fully opaque colours.
      geom_line(
        data = data_pathogens_i,
        aes(
          x = week,
          y = y_data,
          group = component,
          colour = component
        ),
        linewidth = 0.45,
        alpha = 1
      ) +
      scale_color_manual(
        values = figure_colours,
        breaks = legend_breaks,
        labels = legend_labels
      ) +
      scale_fill_manual(
        values = pathogen_colors,
        breaks = names(pathogen_colors),
        guide = "none"
      ) +
      scale_x_date(
        breaks = seq(
          as.Date("2013-01-01"),
          as.Date("2020-01-01"),
          by = "2 years"
        ),
        date_labels = "%Y",
        limits = range(week_seq)
      ) +
      scale_y_continuous(
        limits = function(x) c(0, max(x, na.rm = TRUE)),
        expand = expansion(mult = c(0, 0.05))
      ) +
      labs(
        title = panel_title,
        x = if (show_x_title) "Year" else NULL,
        y = if (show_y_title) "Weekly RTI hospital admissions" else NULL,
        colour = NULL
      ) +
      theme_minimal(base_size = 9, base_family = figure_font) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(
          colour = "grey92",
          linewidth = 0.2
        ),
        axis.line = element_line(
          colour = "grey50",
          linewidth = 0.2
        ),
        axis.ticks = element_line(
          colour = "grey50",
          linewidth = 0.2
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
        plot.title = element_text(
          size = 10,
          hjust = 0.5,
          margin = margin(b = 4)
        ),
        plot.margin = margin(t = 5, r = 4, b = 3, l = 4),
        legend.position = if (show_legend) "right" else "none",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.text = element_text(
          size = 9,
          colour = "grey20"
        ),
        legend.key.width = grid::unit(7, "mm"),
        legend.key.height = grid::unit(4, "mm"),
        legend.margin = margin(0, 0, 0, 0)
      ) +
      guides(
        colour = guide_legend(
          ncol = 1,
          override.aes = list(
            alpha = 1,
            linewidth = 0.6
          )
        )
      )
  }
  
  # Make the eight panels. Add shared axis titles below.
  p_02m <- make_age_panel("0-2m")
  p_35m <- make_age_panel("3-5m")
  p_611m <- make_age_panel("6-11m")
  p_12y <- make_age_panel("1-2y")
  p_34y <- make_age_panel("3-4y")
  p_6574 <- make_age_panel("65-74y")
  p_7584 <- make_age_panel("75-84y")
  p_85p <- make_age_panel("85+y")
  
  # Extract the shared five-entry legend.
  p_for_legend <- make_age_panel(
    "3-4y",
    show_legend = TRUE
  )
  
  legend_hmc <- cowplot::get_legend(p_for_legend)
  
  # Empty panel containing only the legend.
  p_legend <- cowplot::ggdraw() +
    cowplot::draw_grob(
      legend_hmc,
      x = 0,
      y = 0,
      width = 1,
      height = 1
    )
  
  # Assemble the 3 x 3 panel layout.
  figure_core <- cowplot::plot_grid(
    p_02m,  p_35m,  p_611m,
    p_12y,  p_34y,  p_legend,
    p_6574, p_7584, p_85p,
    ncol = 3,
    align = "hv"
  )
  
  # One common x-axis title beneath the panels.
  figure_with_x_title <- cowplot::plot_grid(
    figure_core,
    cowplot::ggdraw() +
      cowplot::draw_label(
        "Year",
        size = 9,
        fontfamily = figure_font
      ),
    ncol = 1,
    rel_heights = c(1, 0.04)
  )
  
  # One common y-axis title on the left.
  plot_hmc <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_label(
        "Weekly RTI hospital admissions",
        angle = 90,
        size = 9,
        fontfamily = figure_font
      ),
    figure_with_x_title,
    ncol = 2,
    rel_widths = c(0.045, 1)
  )
  
  if (show_figures) print(plot_hmc)
  
  # Editable vector PDF for submission.
  ggsave(
    filename = file.path(
      analysis_figure_dir,
      "attribution_hmc.pdf"
    ),
    plot = plot_hmc,
    device = grDevices::cairo_pdf,
    width = figure_width_mm,
    height = figure_height_mm,
    units = "mm",
    bg = "white"
  )
  
  # PNG for convenient viewing.
  ggsave(
    filename = file.path(
      analysis_figure_dir,
      "attribution_hmc.png"
    ),
    plot = plot_hmc,
    dpi = 600,
    width = figure_width_mm,
    height = figure_height_mm,
    units = "mm",
    bg = "white"
  )
  
  message("Saved the attribution time-series figure (PDF and PNG).")
} else {
  message("Attribution figure skipped.")
}