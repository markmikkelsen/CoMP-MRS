# Load CoMP-MRS data and perform initial cleaning and outlier removal

LoadData <- function(csv_file = "CoMP_MRS_Rstats_input.csv",  verbose = TRUE) {
  
  # Load the participants.csv file with appropriate column types
  DATA <- read_csv(
    file = csv_file,
    col_types = list(
      CompID = col_factor(),
      SiteID = col_factor(),
      DP = col_factor(),
      AnimalID = col_factor(),
      AnimalSpecies = col_factor(),
      AnimalStrain = col_factor(),
      AnimalAge = col_double(),
      AnimalSex = col_factor(),
      AnimalWeight = col_double(),
      MRvendor = col_factor(),
      MRfield = col_double(),
      MRsequence = col_factor(),
      MRbrainregion = col_factor(),
      MRvoxelvolume = col_double(),
      MRaverages = col_double(),
      MRsoftwareversion = col_factor(),
      MRcoildetail = col_factor(),
      MRSsw = col_double(),
      MRSnpts = col_double(),
      MRSTE = col_double(),
      MRSTR = col_double(),
      MRSshimmethod = col_factor(),
      LW = col_double(),
      SNR = col_double(),
      SNR_LW_Ratio = col_double(),
      CompCheck = col_factor()
    )
  )
  
  # Remove "compMR" prefix from CompID to make it more concise (optional)
  DATA <- DATA %>%
    mutate(CompID = str_remove_all(CompID, "compMR"))
  
  
  # Spectral quality metrics and normalization --------------------------------
  
  DATA <- DATA %>%
    mutate(
      SNR_LW_Ratio_norm   = SNR_LW_Ratio / (sqrt(MRaverages) * MRvoxelvolume),
      SNR_LW_Product_norm = (LW * SNR) / (sqrt(MRaverages) * MRvoxelvolume),
      LW_norm             = LW / (sqrt(MRaverages) * MRvoxelvolume),
      SNR_norm            = SNR / (sqrt(MRaverages) * MRvoxelvolume)
    )
  
  
  # Clean up data -------------------------------------------------------------
  
  # Ensure that all missing values are represented as NA (in case there are any
  # other representations of missing data in the original CSV)
  DATA[is.na(DATA)] <- NA
  
  ### Remove datasets that failed QC ------------------------------------------
  
  DATA <- DATA %>%
    filter(CompCheck != "fail") %>%
    select(-CompCheck) # Remove the CompCheck column as it's no longer needed
  
  ### Remove outliers based on the MAD method ---------------------------------
  # This will also remove any rows where LW or SNR is NA, which is appropriate for outlier removal
  
  DATA_clean <- mad_outlier_removal(
    data      = DATA,
    outcome   = "SNR_LW_Ratio_norm",
    groups    = "DP",  # outer → inner
    strategy  = "group",  # options: "global", "group", "multilevel"
    threshold = 2.5,
    verbose = verbose
  )
  DATA_orig <- DATA
  DATA <- DATA_clean$data_clean # Update DATA to the cleaned dataset
  
  # Optional but recommended:
  # aggregate to one row per DP for cleaner plots
  DATA_DP <- DATA %>%
    dplyr::group_by(
      DP, SiteID, AnimalSpecies, AnimalSex,
      MRvendor, MRfield, MRsequence, MRbrainregion
    ) %>%
    dplyr::summarise(
      LW_norm = mean(LW_norm, na.rm = TRUE),
      SNR_norm = mean(SNR_norm, na.rm = TRUE),
      SNR_LW_Product_norm = mean(SNR_LW_Product_norm, na.rm = TRUE),
      SNR_LW_Ratio_norm = mean(SNR_LW_Ratio_norm, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(
    list(
      data       = DATA,
      data_orig  = DATA_orig,
      data_by_dp = DATA_DP
    )
  )
  
}
