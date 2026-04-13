# Descriptive statistics for each grouping variable and their combinations.

DescripStats <- function(data) {
  
  STATS <- list()
  
  ### Descriptive stats function ----------------------------------------------
  
  make_stats <- function(data, group_var, value_var, prefix) {
    data %>%
      dplyr::group_by(across(all_of(group_var))) %>%
      dplyr::summarise(
        mean = mean(.data[[value_var]], na.rm = TRUE),
        sd   = sd(.data[[value_var]], na.rm = TRUE),
        cv   = cv(.data[[value_var]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      rename(
        !!paste0("mean", prefix) := mean,
        !!paste0("sd", prefix)   := sd,
        !!paste0("cv", prefix)   := cv
      )
  }
  
  ### DP ----------------------------------------------------------------------
  
  STATS$DP <- list(
    LW      = make_stats(data, "DP", "LW_norm", "LW"),
    SNR     = make_stats(data, "DP", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "DP", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "DP", "SNR_LW_Product_norm", "Product")
  )
  
  ### SiteID ------------------------------------------------------------------
  
  STATS$SiteID <- list(
    LW      = make_stats(data, "SiteID", "LW_norm", "LW"),
    SNR     = make_stats(data, "SiteID", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "SiteID", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "SiteID", "SNR_LW_Product_norm", "Product")
  )
  
  ### AnimalSpecies -----------------------------------------------------------
  
  STATS$AnimalSpecies <- list(
    LW      = make_stats(data, "AnimalSpecies", "LW_norm", "LW"),
    SNR     = make_stats(data, "AnimalSpecies", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "AnimalSpecies", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "AnimalSpecies", "SNR_LW_Product_norm", "Product")
  )
  
  ### MRvendor ----------------------------------------------------------------
  
  STATS$MRvendor <- list(
    LW      = make_stats(data, "MRvendor", "LW_norm", "LW"),
    SNR     = make_stats(data, "MRvendor", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "MRvendor", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "MRvendor", "SNR_LW_Product_norm", "Product")
  )
  
  ### MRfield -----------------------------------------------------------------
  
  STATS$MRfield <- list(
    LW      = make_stats(data, "MRfield", "LW_norm", "LW"),
    SNR     = make_stats(data, "MRfield", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "MRfield", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "MRfield", "SNR_LW_Product_norm", "Product")
  )
  
  ### MRsequence --------------------------------------------------------------
  
  STATS$MRsequence <- list(
    LW      = make_stats(data, "MRsequence", "LW_norm", "LW"),
    SNR     = make_stats(data, "MRsequence", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "MRsequence", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "MRsequence", "SNR_LW_Product_norm", "Product")
  )
  
  ### MRbrainregion -----------------------------------------------------------
  
  STATS$MRbrainregion <- list(
    LW      = make_stats(data, "MRbrainregion", "LW_norm", "LW"),
    SNR     = make_stats(data, "MRbrainregion", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "MRbrainregion", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "MRbrainregion", "SNR_LW_Product_norm", "Product")
  )
  
  ### AnimalSpecies x MRvendor ------------------------------------------------
  
  STATS$AnimalSpecies_MRvendor <- data %>%
    dplyr::group_by(MRvendor, AnimalSpecies) %>%
    dplyr::summarise(
      meanLW      = mean(LW_norm, na.rm = TRUE),
      sdLW        = sd(LW_norm, na.rm = TRUE),
      cvLW        = cv(LW_norm, na.rm = TRUE),
      meanSNR     = mean(SNR_norm, na.rm = TRUE),
      sdSNR       = sd(SNR_norm, na.rm = TRUE),
      cvSNR       = cv(SNR_norm, na.rm = TRUE),
      meanRatio   = mean(SNR_LW_Ratio_norm, na.rm = TRUE),
      sdRatio     = sd(SNR_LW_Ratio_norm, na.rm = TRUE),
      cvRatio     = cv(SNR_LW_Ratio_norm, na.rm = TRUE),
      meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE),
      sdProduct   = sd(SNR_LW_Product_norm, na.rm = TRUE),
      cvProduct   = cv(SNR_LW_Product_norm, na.rm = TRUE),
      .groups = "drop"
    )
  
  ### AnimalSpecies x MRbrainregion -------------------------------------------
  
  STATS$AnimalSpecies_MRbrainregion <- data %>%
    dplyr::group_by(AnimalSpecies, MRbrainregion) %>%
    dplyr::summarise(
      meanLW      = mean(LW_norm, na.rm = TRUE),
      sdLW        = sd(LW_norm, na.rm = TRUE),
      cvLW        = cv(LW_norm, na.rm = TRUE),
      meanSNR     = mean(SNR_norm, na.rm = TRUE),
      sdSNR       = sd(SNR_norm, na.rm = TRUE),
      cvSNR       = cv(SNR_norm, na.rm = TRUE),
      meanRatio   = mean(SNR_LW_Ratio_norm, na.rm = TRUE),
      sdRatio     = sd(SNR_LW_Ratio_norm, na.rm = TRUE),
      cvRatio     = cv(SNR_LW_Ratio_norm, na.rm = TRUE),
      meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE),
      sdProduct   = sd(SNR_LW_Product_norm, na.rm = TRUE),
      cvProduct   = cv(SNR_LW_Product_norm, na.rm = TRUE),
      .groups = "drop"
    )
  
  ### All ---------------------------------------------------------------------
  
  STATS$all <- list(
    meanLW      = mean(data$LW_norm, na.rm = TRUE),
    sdLW        = sd(data$LW_norm, na.rm = TRUE),
    cvLW        = cv(data$LW_norm, na.rm = TRUE),
    meanSNR     = mean(data$SNR_norm, na.rm = TRUE),
    sdSNR       = sd(data$SNR_norm, na.rm = TRUE),
    cvSNR       = cv(data$SNR_norm, na.rm = TRUE),
    meanRatio   = mean(data$SNR_LW_Ratio_norm, na.rm = TRUE),
    sdRatio     = sd(data$SNR_LW_Ratio_norm, na.rm = TRUE),
    cvRatio     = cv(data$SNR_LW_Ratio_norm, na.rm = TRUE),
    meanProduct = mean(data$SNR_LW_Product_norm, na.rm = TRUE),
    sdProduct   = sd(data$SNR_LW_Product_norm, na.rm = TRUE),
    cvProduct   = cv(data$SNR_LW_Product_norm, na.rm = TRUE)
  )
  
  return(STATS)
  
}
