# R script to visualize and statistically analyze data from CoMP-MRS
#
# Authors:
#   Mark Mikkelsen, Ph.D. (mam4041@med.cornell.edu)
#   Diana G. Rotaru, Ph.D. (diana.rotaru@meduniwien.ac.at)
#
# Last updated: 2026-04-08


# Plots -----------------------------------------------------------------------

### Faceted boxplots ----------------------------------------------------------

make_boxplot_vendor_species <- function(data, y_var, y_label, file_name) {
  
  plot_data <- data %>%
    dplyr::filter(
      !is.na(MRvendor),
      !is.na(AnimalSpecies),
      !is.na(.data[[y_var]])
    ) %>%
    dplyr::mutate(
      MRvendor = as.factor(MRvendor),
      AnimalSpecies = as.factor(AnimalSpecies)
    )
  
  p <- ggplot(plot_data, aes(x = MRvendor, y = .data[[y_var]], fill = MRvendor)) +
    geom_boxplot(
      width = 0.65,
      alpha = 0.9,
      outlier.shape = 21,
      outlier.size = 2.0,
      outlier.stroke = 0.3,
      color = "black",
      linewidth = 0.4
    ) +
    geom_jitter(
      width = 0.12,
      alpha = 0.4,
      size = 1.6,
      color = "black"
    ) +
    facet_wrap(~ AnimalSpecies, scales = "free_x") +
    labs(
      title = paste(y_label, "by vendor and animal species"),
      x = "Vendor",
      y = y_label
    ) +
    scale_fill_manual(values = get_discrete_colors(nlevels(plot_data$MRvendor))) +
    theme_comp()
  
  print(p)
  
  ggsave(
    filename = file.path(plots_dir, file_name),
    plot = p,
    width = 8,
    height = 5,
    units = "in",
    dpi = 300
  )
  
  return(p)
}

if (show_facet_plots) {
  
  facet_specs <- list(
    list(
      var = "LW_norm",
      label = "Normalized LW",
      file = "boxplot_LW_norm_by_MRvendor_AnimalSpecies.png"
    ),
    list(
      var = "SNR_norm",
      label = "Normalized SNR",
      file = "boxplot_SNR_norm_by_MRvendor_AnimalSpecies.png"
    ),
    list(
      var = "SNR_LW_Product_norm",
      label = "Normalized SNR×LW product",
      file = "boxplot_SNR_LW_Product_norm_by_MRvendor_AnimalSpecies.png"
    ),
    list(
      var = "SNR_LW_Ratio_norm",
      label = "Normalized SNR/LW ratio",
      file = "boxplot_SNR_LW_Ratio_norm_by_MRvendor_AnimalSpecies.png"
    )
  )
  
  facet_plots <- lapply(facet_specs, function(x) {
    make_boxplot_vendor_species(
      data = DATA_DP,
      y_var = x$var,
      y_label = x$label,
      file_name = x$file
    )
  })
  
  names(facet_plots) <- sapply(facet_specs, `[[`, "var")
}



# End -------------------------------------------------------------------------
