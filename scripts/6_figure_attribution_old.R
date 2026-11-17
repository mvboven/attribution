# Observed admissions, fitted totals, and pathogen contributions.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
if (make_hmc_figure) {
  require_fit()
  # Use the actual dates of the ISO weeks
  week_seq <- as.Date(data$week_start)

  stopifnot(
    length(week_seq) == n_data,
    nrow(y_data) == n_data,
    all(diff(week_seq) == 7),
    "2015-W53" %in% data$ISOweek
  )

  # Overall fitted incidence: posterior median only
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
    fit_hmc, week_seq, age_groups, pathogen_lookup$pathogen)

  # Observed admissions
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

  # Combine plotting data
  plot_data_hmc <- bind_rows(
    observed_data,
    overall_summary,
    pathogen_summary
  ) %>%
    mutate(
      age_group = factor(
        age_group,
        levels = age_groups
      ),
      component = factor(
        component,
        levels = c(
          "Data",
          "Overall",
          names(pathogen_colors)
        )
      )
    )

  # Consistency checks
  stopifnot(
    nrow(overall_summary) ==
      num_types * n_data,

    nrow(pathogen_summary) ==
      num_path * num_types * n_data,

    !anyNA(pathogen_summary$y_data),
    !anyNA(pathogen_summary$y_lo),
    !anyNA(pathogen_summary$y_hi)
  )

  # Separate data used by the different plot layers
  plot_data_observed <- plot_data_hmc %>%
    filter(component == "Data")

  plot_data_overall <- plot_data_hmc %>%
    filter(component == "Overall")

  plot_data_pathogens <- plot_data_hmc %>%
    filter(component %in% names(pathogen_colors))

  # helper function for one age-group panel
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

    ggplot() +
      geom_line(
        data = data_overall_i,
        aes(
          x = week,
          y = y_data,
          group = component,
          color = component
        ),
        linewidth = 0.9,
        alpha = 0.9
      ) +
      geom_point(
        data = data_observed_i,
        aes(
          x = week,
          y = y_data
        ),
        alpha = 0.9,
        size = 0.3
      ) +
      geom_ribbon(
        data = data_pathogens_i,
        aes(
          x = week,
          ymin = y_lo,
          ymax = y_hi,
          fill = component,
          group = component
        ),
        alpha = 0.15,
        colour = NA
      ) +
      geom_line(
        data = data_pathogens_i,
        aes(
          x = week,
          y = y_data,
          group = component,
          color = component
        ),
        linewidth = 0.5,
        alpha = 0.6
      ) +
      scale_color_manual(
        values = cbPalette,
        breaks = names(pathogen_colors)
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
        limits = function(x) c(0, max(x, na.rm = TRUE))
      ) +
      labs(
        title = selected_age_group,
        x = if (show_x_title) "Year" else NULL,
        y = if (show_y_title) "Weekly RTI hospital admissions" else NULL,
        color = NULL
      ) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = if (show_legend) "right" else "none",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.text = element_text(size = 10),
        legend.key.size = grid::unit(0.4, "cm"),
        axis.title.x = element_text(),
        axis.title.y = element_text()
      ) +
      guides(
        color = guide_legend(ncol = 1)
      )
  }

  # make the eight age-group panels
  p_02m <- make_age_panel(
    "0-2m",
    show_x_title = FALSE,
    show_y_title = FALSE
  )

  p_35m <- make_age_panel(
    "3-5m",
    show_x_title = FALSE,
    show_y_title = FALSE
  )

  p_611m <- make_age_panel(
    "6-11m",
    show_x_title = FALSE,
    show_y_title = FALSE
  )

  p_12y <- make_age_panel(
    "1-2y",
    show_x_title = FALSE,
    show_y_title = FALSE
  )

  p_34y <- make_age_panel(
    "3-4y",
    show_x_title = FALSE,
    show_y_title = FALSE
  )

  p_6574 <- make_age_panel(
    "65-74y",
    show_x_title = TRUE,
    show_y_title = FALSE
  )

  p_7584 <- make_age_panel(
    "75-84y",
    show_x_title = TRUE,
    show_y_title = FALSE
  )

  p_85p <- make_age_panel(
    "85+y",
    show_x_title = TRUE,
    show_y_title = FALSE
  )

  # make a temporary plot containing the legend
  p_for_legend <- make_age_panel(
    "3-4y",
    show_legend = TRUE,
    show_x_title = FALSE,
    show_y_title = FALSE
  )

  # extract legend only
  legend_hmc <- cowplot::get_legend(
    p_for_legend +
      theme(
        legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(size = 14),
        legend.key.size = grid::unit(0.7, "cm"),
        legend.spacing.y = grid::unit(0.25, "cm")
      )
  )

  # make an empty panel containing only the legend
  p_legend <- cowplot::ggdraw() +
    cowplot::draw_grob(
      legend_hmc,
      x = 0,
      y = 0,
      width = 1,
      height = 1
    )

  # assemble the 3 x 3 panel layout
  figure_core <- cowplot::plot_grid(
    p_02m,  p_35m,  p_611m,
    p_12y,  p_34y,  p_legend,
    p_6574, p_7584, p_85p,
    ncol = 3,
    align = "hv"
  )

  # add one common y-axis title on the left
  plot_hmc <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_label(
        "Weekly RTI hospital admissions",
        angle = 90,
        size = 11
      ),
    figure_core,
    ncol = 2,
    rel_widths = c(0.04, 1)
  )

  if (show_figures) print(plot_hmc)

  # save figure
  ggsave(
    filename = file.path(
      analysis_figure_dir,
      "attribution_hmc.png"
    ),
    plot = plot_hmc,
    dpi = 300,
    height = 6.8,
    width = 10,
    bg = "white"
  )


  message("Saved the attribution time-series figure.")
} else message("Attribution figure skipped.")
