# Set up
setwd("YOUR PATH")
.libPaths("YOUR PATH")

time1 <- Sys.time()

# Load packages
library(mice)
library(tidyverse)
library(bkmr)
library(bkmrhat)

# load manipulated data----------------------------------------------------------

load("PATH TO IMPUTED DATA SETS")
# df_afimp below is data after imputation

# set date
systime <- Sys.Date()

#####################################
##                                 ##
##              BKMR               ##-------------------------------------------
##                                 ##
#####################################


#################################################################################
# Function to perform BKMR analysis for a given outcome
# Function to perform BKMR analysis for a given outcome
run_bkmr_analysis <- function(outcome, Z_matrix, X_matrix, chains, iterations, group, sel) {
  # Set random seed for reproducibility
  set.seed(2024)
  
  # Run BKMR model
  fit <- kmbayes_parallel(nchains = chains, y = outcome, Z = Z_matrix, X = X_matrix, 
                          iter = iterations, verbose = FALSE, varsel = TRUE, groups = group)
  
  # Combine multiple chains
  fit_combined <- kmbayes_combine(fit)
  
  
  # Summarize model outputs
  univar_pred <- PredictorResponseUnivar(fit = fit_combined, sel = sel)
  bivar_pred <- PredictorResponseBivar(fit = fit_combined, sel = sel, min.plot.dist = 1)
  bivar_levels <- PredictorResponseBivarLevels(pred.resp.df = bivar_pred, 
                                               Z = Z_matrix, qs = c(0.25, 0.5, 0.75))
  overall_risk <- OverallRiskSummaries(fit= fit_combined, y = outcome, Z = Z_matrix, 
                                       X = X_matrix, sel = sel, qs = seq(0.25, 0.75, by = 0.05), 
                                       q.fixed = 0.5, method = "exact")
  singvar_risk <- SingVarRiskSummaries(fit = fit_combined, y = outcome, Z = Z_matrix, 
                                       X = X_matrix, qs.diff = c(0.25, 0.75), 
                                       q.fixed = c(0.25, 0.50, 0.75), method = "exact")
  risk_int <- SingVarIntSummaries(fit = fit_combined, qs.diff = c(0.25, 0.75))
  
  
  # Return all results as a list
  list(
    fit_combined = fit_combined,
    univar_pred = univar_pred,
    bivar_pred = bivar_pred,
    bivar_levels = bivar_levels,
    overall_risk = overall_risk,
    singvar_risk = singvar_risk,
    risk_int = risk_int
  )
}

#################################################################################
# Define variables

# extract the 30th imputed data set
# BKMR ONLY WORKS WITH COMPLETE DATASET, NO MISSING VALUES ALLOWED
dd_m <- df_afimp %>% mice::complete(30)

# Exposure matrix
Z_preg <- dd_m %>% 
  select(
    "YOUR LIST OF METABOLITES"
  ) %>% 
  as.matrix()

# Z-score exposure matrix to put on same scale 
Z_preg_scale <- scale(Z_preg)


# Covariates - mothers
cov <- c(
  "YOUR LIST"
)

# covariates need to be numeric
X_preg <- dd_m %>% select(all_of(cov)) %>% 
  mutate(across(where(is.factor), as.numeric)) %>% 
  mutate(across(where(is.numeric), scale)) %>% 
  as.matrix()

#################################################################################################
time1 <- Sys.time()

future::plan(strategy = future::multisession)

# List of cognitive outcomes
outcomes <- list(
  voc = dd_m$WISC13_Voc_Tscore
)


# Initialize list to store results
results <- list()

# Loop over each outcome and run analysis
for (outcome_name in names(outcomes)) {
  outcome <- outcomes[[outcome_name]]
  
  # Run the BKMR analysis and store results, 50000 iterations
  # CHANGE PARAMETERS: change "chains=" to how many cores you want to use, for testing run, recommend to try with "chain = 2, iterations = 1000". If you change this parameter, also change sel = seq(1500, 2000, by = 2)
  results[[outcome_name]] <- run_bkmr_analysis(outcome, Z_preg_scale, X_preg, chains = 8, iterations = 6250, 
                                               group = c(1,2,2,2,2,3,3,4,5,5,6,7,7,7,7,7),
                                               sel = seq(40000, 50000, by = 5))
}

time2 <- Sys.time()

# running time for one outcome
time2-time1

#################################################################################################
# Save image
save(BKMRfits, file = paste0("YOUR PATH/FILE NAME_", systime, ".RData"))
