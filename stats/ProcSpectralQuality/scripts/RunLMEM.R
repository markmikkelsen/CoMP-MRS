# Linear mixed-effects modeling

RunLMEM <- function(data, dv = "", rand_ef = list(), fix_ef = "") {
  
  cat("\n── Running linear mixed-effects model (LMEM) ──\n")
  
  # Fit dynamic lmer models ---------------------------------------------------
  
  fit_lmer <- function(response, random_effects, fixed_effects, data, ...) {
    #' Dynamically build and fit a linear mixed model using lmer
    
    #'
    #' @param response      Character string: name of the response variable
    #' @param fixed_effects Character vector: names of fixed-effect predictors
    #'                      (e.g., c("group", "age", "sex"))
    #' @param random_effects Named list: each element defines a random-effect term.
    #'                       The name is the grouping factor; the value is a
    #'                       character vector of random slopes (use "1" or NULL
    #'                       for intercept-only).
    #'                       Examples:
    #'                         list(subject = "1")
    #'                         list(subject = c("1", "time"), site = "1")
    #' @param data           Data frame
    #' @param ...            Additional arguments passed to lmer()
    #' @return An lmerMod object
    
    # --- validate inputs ---
    stopifnot(
      is.character(response), length(response) == 1,
      is.list(random_effects), length(random_effects) >= 1,
      !is.null(names(random_effects)),
      is.character(fixed_effects), length(fixed_effects) >= 1,
      is.data.frame(data)
    )
    
    # --- fixed-effects part ---
    if (length(fixed_effects) == 1 && nchar(fixed_effects) == 0) {
      fixed_part <- "1"
    } else {
      fixed_part <- paste(fixed_effects, collapse = " + ")
    }
    
    # --- random-effects part ---
    random_part <- vapply(names(random_effects), function(grp) {
      slopes <- random_effects[[grp]]
      if (is.null(slopes)) slopes <- "1"
      terms <- paste(slopes, collapse = " + ")
      sprintf("(%s | %s)", terms, grp)
    }, character(1))
    
    random_part <- paste(random_part, collapse = " + ")
    
    # --- assemble and fit ---
    formula_str <- paste(response, "~", fixed_part, "+", random_part)
    f <- as.formula(formula_str)
    
    message("Fitting model: ", deparse(f))
    
    m <- lme4::lmer(f, data = data, ...)
    return(m)
  }
  
  
  # Selection of optimizers to use --------------------------------------------
  
  optWrap.nloptwrap     <- lmerControl(optimizer = "nloptwrap")
  optWrap.optimx.nlminb <- lmerControl(optimizer = "optimx", optCtrl = list(method = "nlminb", eval.max = 1e5))
  optWrap.minqa.bobyqa  <- lmerControl(optimizer = "bobyqa")
  
  optimToUse <- optWrap.minqa.bobyqa
  
  
  # Remove rows with missing values -------------------------------------------
  
  data <- drop_na(data, .data[[dv]])
  
  
  # Models --------------------------------------------------------------------
  
  M <- fit_lmer(
    response = dv,
    random_effects = rand_ef,
    fixed_effects = fix_ef,
    data = data,
    REML = TRUE,
    control = optimToUse
  )
  
  
  # M.SNRLWrationorm.0.1.a <- lme4::lmer(dv ~ (1 | MRbrainregion),
  #                                      data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.0.1.b <- lme4::lmer(dv ~ (1 | MRsequence),
  #                                      data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.0.2 <- lme4::lmer(dv ~ (1 | MRvendor), data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.0.3 <- lme4::lmer(dv ~ (1 | SiteID), data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.0.4 <- lme4::lmer(dv ~ (1 | DP) + (1 | SiteID), data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.0.full <- lme4::lmer(dv ~ (1 | DP),     # (1 | MRvendor) + (1 | SiteID) only account for negligible variance
  #                                       data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.1.0 <- lme4::lmer(dv ~ (1 | AnimalSpecies) + (1 | DP), # (1 | AnimalSpecies) only accounts for negligible variance
  #                                    data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.2.0 <- lme4::lmer(dv ~ (1 | MRbrainregion) + (1 | DP), # (1 | MRbrainregion) only accounts for negligible variance
  #                                    data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.3.0 <- lme4::lmer(dv ~ (1 | MRsequence) + (1 | DP),
  #                                    data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.4.0 <- lme4::lmer(dv ~ MRfield + (1 | DP),
  #                                    data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.4.1 <- lme4::lmer(dv ~ MRfield + (MRfield | DP),
  #                                    data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.5.0 <- lme4::lmer(dv ~  AnimalAge + (1 | DP),
  #                                    data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.5.1 <- lme4::lmer(dv ~  AnimalAge + (AnimalAge | DP),
  #                                    data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.6.0 <- lme4::lmer(dv ~  AnimalSex + (1 | DP),
  #                                    data = data, REML = FALSE, control = optimToUse)
  # M.SNRLWrationorm.6.1 <- lme4::lmer(dv ~  AnimalSex + (AnimalSex | DP),
  #                                    data = data, REML = FALSE, control = optimToUse)
  
  

  # Inferential tests ---------------------------------------------------------
  
  # Compare null model to model with DP random effect
  # LRT.1 <- anova(M.SNRLWrationorm.0.full, M.SNRLWrationorm.3.0)
  
  # Bootstrapping confidence intervals for fixed effects (intercept)
  # b.par1 <- bootMer(Y.M0.3, fixef, nsim = 1e4) # bootstrap fixed effects
  # b.par2 <- bootMer(Y.M1.3, fixef, nsim = 1e4)
  # boot.ci(b.par1, conf = 0.95, type = "bca", index = 1)
  # boot.ci(b.par2, conf = 0.95, type = "bca", index = 1)
  
  # Alternative bootstrapping approach
  # confint(M.SNRLWrationorm.4.0, parm = c(3,4), level = 0.95, method = "boot",
  #         nsim = 500, boot.type = "perc")
  
  # if (run_pbkrtest) {
  #   
  #   # Set up parallel computation for bootstrapping
  #   nc <- detectCores()
  #   clus <- makeCluster(rep("localhost", nc))
  #   clusterEvalQ(clus, {
  #     library(lme4)
  #     library(pbkrtest)
  #   })
  #   clusterExport(clus, c("nloptr", "defaultControl", "nloptFun",
  #                         "data", "M.SNRLWrationorm.0.full", "M.SNRLWrationorm.3.0"))
  #   
  #   # Bootstrap likelihood ratio tests
  #   PB_LRT.1 <- PBmodcomp(M.SNRLWrationorm.3.0, M.SNRLWrationorm.0.full, nsim = 2e3, cl = clus)
  #   
  #   stopCluster(clus)
  #   
  # }
  
  return(M)
  
}
