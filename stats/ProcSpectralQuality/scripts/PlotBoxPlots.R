# Box plots for each  (with jittered points)

PlotBoxPlots <- function(data, out_dir, y_vars, x_vars) {

  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  } else {
    unlink(out_dir, recursive = TRUE)
    dir.create(out_dir)
  }

  # Box plot function ---------------------------------------------------------

  make_box_plot <- function(data, x_var, y_var, y_label, group_var, out_dir, file_name) {
    every_nth <- function(n) {
      function(x) x[seq(1, length(x), by = n)]
    }

    plot_data <- data %>%
      dplyr::filter(
        !is.na(.data[[x_var]]),
        !is.na(.data[[y_var]]),
        !is.na(.data[[group_var]])
      ) %>%
      dplyr::arrange(.data[[group_var]], .data[[x_var]]) %>%
      dplyr::mutate(
        !!x_var := factor(.data[[x_var]], levels = unique(.data[[x_var]]))
      )
    
    n_levels <- nlevels(plot_data[[x_var]])
    nth <- max(1, round(n_levels / 10))
    
    if (n_levels > 10) {
      a <- 45
      j <- 1
    }  else {
      a <- 0
      j <- 0.5
    }
    
    p <- ggplot(
      plot_data,
      aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
      )
    ) +
      geom_boxplot(
        fill = "#EA7E2D",
        width = 0.65,
        alpha = 0.9,
        outlier.shape = 21,
        outlier.size = 1.8,
        outlier.stroke = 0.3,
        color = "black",
        linewidth = 0.4
      ) +
      geom_point(
        alpha = 0.85,
        color = "black",
        shape = 16
      ) +
      labs(
        title = paste(y_label, "by", x_var),
        x = x_var,
        y = y_label
      ) +
      scale_x_discrete(
        breaks = every_nth(nth)
      ) +
      theme_comp() +
      theme(
        axis.text.x = element_text(angle = a, hjust = j)
      )

    print(p)

    ggsave(
      file.path(out_dir, file_name),
      plot = p,
      width = 11,
      height = 5,
      units = "in",
      dpi = 300
    )
  }


  # Make box plots ---------------------------------------------------------

  # Create standard boxplots
  all_boxplots <- list()

  for (y in y_vars) {
    for (x in x_vars) {
      file_name <- paste0("boxplot_", y$var, "_by_", x, ".png")
      key <- paste(y$var, x, sep = "_by_")

      all_boxplots[[key]] <- make_box_plot(
        data = data,
        x_var = x,
        y_var = y$var,
        y_label = y$label,
        group_var = "DP",
        out_dir = out_dir,
        file_name = file_name
      )
    }
  }
  
  return(all_boxplots)

}
