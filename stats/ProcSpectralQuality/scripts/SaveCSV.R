### Save descriptive stats tables

SaveCSV <- function(data, out_dir) {
  
  write_csv(data$DP$LW, file.path(out_dir, "stats_DP_LW.csv"))
  write_csv(data$DP$SNR, file.path(out_dir, "stats_DP_SNR.csv"))
  write_csv(data$DP$Ratio, file.path(out_dir, "stats_DP_Ratio.csv"))
  write_csv(data$DP$Product, file.path(out_dir, "stats_DP_Product.csv"))
  
  write_csv(data$SiteID$LW, file.path(out_dir, "stats_SiteID_LW.csv"))
  write_csv(data$SiteID$SNR, file.path(out_dir, "stats_SiteID_SNR.csv"))
  write_csv(data$SiteID$Ratio, file.path(out_dir, "stats_SiteID_Ratio.csv"))
  write_csv(data$SiteID$Product, file.path(out_dir, "stats_SiteID_Product.csv"))
  
  write_csv(data$AnimalSpecies$LW, file.path(out_dir, "stats_AnimalSpecies_LW.csv"))
  write_csv(data$AnimalSpecies$SNR, file.path(out_dir, "stats_AnimalSpecies_SNR.csv"))
  write_csv(data$AnimalSpecies$Ratio, file.path(out_dir, "stats_AnimalSpecies_Ratio.csv"))
  write_csv(data$AnimalSpecies$Product, file.path(out_dir, "stats_AnimalSpecies_Product.csv"))
  
  write_csv(data$MRvendor$LW, file.path(out_dir, "stats_MRvendor_LW.csv"))
  write_csv(data$MRvendor$SNR, file.path(out_dir, "stats_MRvendor_SNR.csv"))
  write_csv(data$MRvendor$Ratio, file.path(out_dir, "stats_MRvendor_Ratio.csv"))
  write_csv(data$MRvendor$Product, file.path(out_dir, "stats_MRvendor_Product.csv"))
  
  write_csv(data$MRfield$LW, file.path(out_dir, "stats_MRfield_LW.csv"))
  write_csv(data$MRfield$SNR, file.path(out_dir, "stats_MRfield_SNR.csv"))
  write_csv(data$MRfield$Ratio, file.path(out_dir, "stats_MRfield_Ratio.csv"))
  write_csv(data$MRfield$Product, file.path(out_dir, "stats_MRfield_Product.csv"))
  
  write_csv(data$MRsequence$LW, file.path(out_dir, "stats_MRsequence_LW.csv"))
  write_csv(data$MRsequence$SNR, file.path(out_dir, "stats_MRsequence_SNR.csv"))
  write_csv(data$MRsequence$Ratio, file.path(out_dir, "stats_MRsequence_Ratio.csv"))
  write_csv(data$MRsequence$Product, file.path(out_dir, "stats_MRsequence_Product.csv"))
  
  write_csv(data$MRbrainregion$LW, file.path(out_dir, "stats_MRbrainregion_LW.csv"))
  write_csv(data$MRbrainregion$SNR, file.path(out_dir, "stats_MRbrainregion_SNR.csv"))
  write_csv(data$MRbrainregion$Ratio, file.path(out_dir, "stats_MRbrainregion_Ratio.csv"))
  write_csv(data$MRbrainregion$Product, file.path(out_dir, "stats_MRbrainregion_Product.csv"))
  
  write_csv(
    data$AnimalSpecies_MRvendor,
    file.path(out_dir, "stats_AnimalSpecies_MRvendor.csv")
  )
  
  write_csv(
    data$AnimalSpecies_MRbrainregion,
    file.path(out_dir, "stats_AnimalSpecies_MRbrainregion.csv")
  )
  
  write_csv(
    tibble(
      meanLW = data$all$meanLW,
      sdLW = data$all$sdLW,
      cvLW = data$all$cvLW,
      meanSNR = data$all$meanSNR,
      sdSNR = data$all$sdSNR,
      cvSNR = data$all$cvSNR,
      meanRatio = data$all$meanRatio,
      sdRatio = data$all$sdRatio,
      cvRatio = data$all$cvRatio,
      meanProduct = data$all$meanProduct,
      sdProduct = data$all$sdProduct,
      cvProduct = data$all$cvProduct
    ),
    file.path(out_dir, "stats_all.csv")
  )
  
}
