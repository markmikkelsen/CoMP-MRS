# Wrapper R script to analyze data from CoMP-MRS
#
# Authors:
#   Mark Mikkelsen, Ph.D. (mam4041@med.cornell.edu)
#   Diana G. Rotaru, Ph.D. (diana.rotaru@meduniwien.ac.at)
#
# Last updated: 2026-04-08

# Initialize ------------------------------------------------------------------

rm(list = ls()) # Clear the current environment
try(dev.off(dev.list()["RStudioGD"]), silent = TRUE) # If using RStudio, this clears all plots
cat("\014") # CTRL+L (clear console)

# Set data directories --------------------------------------------------------

##### CHANGE THESE DIRECTORIES AS NEEDED #####

base_dir  <- dirname(rstudioapi::getSourceEditorContext()$path)
data_dir  <- file.path(base_dir, "data")
deriv_dir <- file.path(data_dir, "derivatives")
plots_dir <- file.path(deriv_dir, "plots")
csv_dir   <- file.path(deriv_dir, "csv")

if (!dir.exists(deriv_dir)) {
  dir.create(deriv_dir)
}

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
} else {
  unlink(plots_dir, recursive = TRUE)
  dir.create(plots_dir)
}

if (!dir.exists(csv_dir)) {
  dir.create(csv_dir)
} else {
  unlink(csv_dir, recursive = TRUE)
  dir.create(csv_dir)
}


# Pre-load code ---------------------------------------------------------------

source(file.path(base_dir, "fun", "mad_outlier_removal.R"))
source(file.path(base_dir, "scripts", "LoadPackages.R"))
source(file.path(base_dir, "scripts", "InlineFun.R"))
source(file.path(base_dir, "scripts", "LoadData.R"))
source(file.path(base_dir, "scripts", "DescripStats.R"))
source(file.path(base_dir, "scripts", "SaveCSV.R"))
source(file.path(base_dir, "scripts", "PlotPieCharts.R"))
source(file.path(base_dir, "scripts", "PlotDotPlots.R"))
source(file.path(base_dir, "scripts", "PlotBoxPlots.R"))
source(file.path(base_dir, "scripts", "RunLMEM.R"))
source(file.path(base_dir, "scripts", "ExtractVPCs.R"))


# Analysis options ------------------------------------------------------------

save_csv               <- FALSE # Set to TRUE to save descriptive statistics tables as CSV files in the derivatives directory
show_pie_charts        <- FALSE # Set to TRUE to create pie charts of categorical variables (e.g., species
show_amcharts          <- FALSE # Set to TRUE to create interactive 3D pie charts using amCharts4
show_dot_plots         <- FALSE  # Set to TRUE to create dot plots of spectral quality metrics by different grouping variables
show_box_plots         <- TRUE  # Set to TRUE to create box plots of spectral quality metrics by different grouping variables
# show_facet_plots       <- FALSE # Set to TRUE to create facet plots of spectral quality metrics by different grouping variables
show_model_diagnostics <- TRUE # Set to TRUE to show model diagnostic plots (e.g., residuals, Q-Q plots) for linear mixed-effects models
calc_VPCs              <- TRUE # Set to TRUE to calculate variance partition coefficients (VPCs) from linear mixed-effects models to assess the proportion of variance explained by each random effect
run_pbkrtest           <- TRUE # Set to TRUE to run parametric bootstrapping using the pbkrtest package 
                               # to compare linear mixed-effects models with different random effects structures
                               # and derive p-values for the added random effects (can be time-consuming with larger datasets)


# Load data -------------------------------------------------------------------
# Also clean up data (incl. outlier removal) and create new variables (e.g., normalized SNR/LW ratio)

DATA <- LoadData(
  csv_file = file.path(data_dir, "CoMP_MRS_Rstats_input.csv"),
  outl_rm_strategy = "group",
  verbose = TRUE
)


# Run descriptive statistics --------------------------------------------------

STATS <- DescripStats(data = DATA$data)


# Save descriptive statistics tables as CSV files -----------------------------

if (save_csv) {
  SaveCSV(data = STATS, out_dir = csv_dir)
}


# Plot pie charts -------------------------------------------------------------

if (show_pie_charts) {
  PlotPieCharts(
    data = DATA$data,
    out_dir = file.path(plots_dir, "pie_charts"),
    plots_dir = plots_dir,
    show_pie_charts = show_pie_charts,
    show_amcharts = show_amcharts
  )
}


# Plot dot plots --------------------------------------------------------------

if (show_dot_plots) {
  PlotDotPlots(
    data = DATA$data,
    data_orig = DATA$data_orig,
    out_dir = file.path(plots_dir, "dot_plots")
    )
}


# Plot box plots --------------------------------------------------------------

if (show_box_plots) {
  
  y_vars <- list(
    #list(var = "LW_norm",             label = "Normalized LW"),
    #list(var = "SNR_norm",            label = "Normalized SNR"),
    #list(var = "SNR_LW_Product_norm", label = "Normalized SNR×LW product"),
    list(var = "SNR_LW_Ratio_norm",   label = "Normalized SNR/LW ratio")
  )
  
  x_vars <- c(
    "DP",
    "SiteID",
    "AnimalSpecies",
    "AnimalSex",
    "MRvendor",
    "MRfield",
    "MRsequence",
    "MRbrainregion"
  )
  
  all_boxplots <- PlotBoxPlots(
    data = DATA$data,
    out_dir = file.path(plots_dir, "box_plots"),
    y_vars = y_vars,
    x_vars = x_vars
  )
  
}


# Run linear mixed-effects modeling -------------------------------------------

### Variables -----------------------------------------------------------------

dv <- "SNR_LW_Ratio_norm"

random_effects <- list(
  M.SNRLWrationorm.0.1 = list(
    DP = "1"
  ),
  M.SNRLWrationorm.0.2 = list(
    SiteID = "1"
  ),
  M.SNRLWrationorm.0.3 = list(
    MRvendor = "1"
  ),
  M.SNRLWrationorm.0.4 = list(
    MRsequence = "1"
  ),
  M.SNRLWrationorm.0.5 = list(
    DP = "1",
    MRsequence = "1"
  )
)

fixed_effects <- list(
  M.SNRLWrationorm.0.1 = c(
    ""
  ),
  M.SNRLWrationorm.0.2 = c(
    ""
  ),
  M.SNRLWrationorm.0.3 = c(
    ""
  ),
  M.SNRLWrationorm.0.4 = c(
    ""
  ),
  M.SNRLWrationorm.0.5 = c(
    ""
  )
)

### Run LMEM models -----------------------------------------------------------

LMEM_MODELS <- lapply(names(random_effects), function(model_name) {
  RunLMEM(
    data = DATA$data,
    dv = dv,
    rand_ef = random_effects[[model_name]],
    fix_ef = fixed_effects[[model_name]]
  )
})
names(LMEM_MODELS) <- names(random_effects)

### Model diagnostics ---------------------------------------------------------

# model <- M.SNRLWrationorm
# 
# if (show_model_diagnostics) {
#   
#   cat("\n── Random-effects diagnostics ──\n")
#   re <- ranef(model, condVar = TRUE)
#   vapply(re, function(g) {
#     print(dotplot(re, scales = "free"))
#   }, list(1))
#   # print(dotplot(ranef(model, condVar = TRUE), scales = "free"))
#   
#   cat("\n── Model performance indices ──\n")
#   m_pef <- model_performance(model)
#   
#   cat("\n── Check model assumptions (performance) ──\n")
#   chk <- check_model(model)
#   print(chk)
#   
# }

### Variance partitioning -----------------------------------------------------

if (calc_VPCs) {
  VPCs <- ExtractVPCs(LMEM_MODELS, verbose = TRUE)
}

### Inference by LRT with parametric bootstrapping ----------------------------

# Compare large model with smaller model to derive p-value for added random effect (e.g., MRsequence)
# Note, models have to be fitted with REML = FALSE for valid comparison by LRT,
# and this can be time-consuming with larger datasets
if (run_pbkrtest) {
  large_model <- LMEM_MODELS$M.SNRLWrationorm.0.5
  small_model <- LMEM_MODELS$M.SNRLWrationorm.0.1
  KR_results  <- pbkrtest::PBmodcomp(large_model, small_model, nsim = 1e3)
  print(KR_results)
}

# Confidence intervals for fixed effects (intercept)
# confint(Y.M1.3, parm=c(3,4), level = 0.95, method = "boot", nsim = 1e3, boot.type = "perc")
# Bootstrapping
# b.par1 <- bootMer(Y.M0.3, fixef, nsim=1e4)
# b.par2 <- bootMer(Y.M1.3, fixef, nsim=1e4)
# boot.ci(b.par1, conf=0.95, type="perc", index=1)
# boot.ci(b.par2, conf=0.95, type="perc", index=1)





