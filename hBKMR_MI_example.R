# Set up
setwd("YOUR PATH")
.libPaths("YOUR PATH")

time1 <- Sys.time()

# Load packages
library(mice)
library(tidyverse)
library(bkmr)
library(bkmrhat) # parallel running extension of BKMR

# load manipulated data----------------------------------------------------------

load("/home/r044073/Chemx_IQ/data/DATA_afimp_manipulated_2024-11-06.RData")

# set date
systime <- Sys.Date()

#####################################
##                                 ##
##            BKMR-MI              ##-------------------------------------------
##                                 ##
#####################################

#################################################################################
# Manipulate imputed data for analysis 

# Covariates - mothers
cov <- c(
  "YOUR LIST"
)

dd_m0 <- df_afimp %>% complete("long", include = T) %>% 
  filter(mis_chemx_preg_iq == 1) %>% 
  select(.imp, .id,
         "YOUR_LIST_OF_METABOLITES", "YOUR_LIST_OF_OUTCOME",
         all_of(cov)
  ) %>% 
  mutate(across(all_of(cov), as.numeric))

dd_m0 <- as.mids(dd_m0)

# Create the 30 imputed data sets
# FOR TEST RUN: RECOMMEND TO USE LESS THAN 5 DATA SETS
num <- 0
repeat{
  num = num+1
  dd <- complete(dd_m0, num)
  assign(paste('dd.',num,sep=""), dd, envir = .GlobalEnv)
  if (num == 30){
    print("repeat loop ends");
    break
  }
}

# Create a list of data frames from dd.1 to dd.30
dd_ls <- lapply(1:30, function(i) get(paste0("dd.", i)))


#################################################################################################
# Analysis with multiple imputed datasets

time1 <- Sys.time()

# Initialize list to store BKMR fits
BKMRfits <- list()

dsnum <- 0
repeat{
  dsnum <- dsnum + 1
  dd_m <- dd_ls[[dsnum]]
  
  # Define variables
  # Exposure matrix
  Z_preg <- dd_m %>% 
    select(
      "YOUR_LIST_OF_METABOLITES"
    ) %>% 
    as.matrix()
  
  # Z-score exposure matrix to put on same scale 
  Z_preg_scale <- scale(Z_preg)
  
  # covariates need to be numeric
  X_preg <- dd_m %>% select(all_of(cov)) %>% 
    mutate(across(where(is.numeric), scale)) %>% 
    as.matrix()
  
  # Outcome
  outcome <- dd_m$WISC13_CD_Tscore
  
  # Set random seed for reproducibility
  set.seed(2024)
  
  future::plan(strategy = future::multisession)
  
  # Run BKMR model
  # CHANGE PARAMETERS: change "chains=" to how many cores you want to use, for testing run, recommend to try with "chain = 2, iterations = 1000". If you change this parameter, also change sel = seq(1500, 2000, by = 2)
  fit <- bkmrhat::kmbayes_parallel(outcome, Z_preg_scale, X_preg, nchains = 8, iter = 6250, 
                                   group = c(1,2,2,2,2,3,3,4,5,5,6,7,7,7,7,7),
                                   verbose = FALSE, varsel = TRUE)
  
  # Combine multiple chains
  BKMRfits[[dsnum]] <- kmbayes_combine(fit)
  
  if (dsnum == 30){
    print("repeat loop ends");
    break
  }
}

time2 <- Sys.time()

# running time for one outcome
time2-time1

#################################################################################################
# Save image
save(BKMRfits, file = paste0("/home/r044073/Chemx_IQ/results/BKMR_MI_EDC_CD_", systime, ".RData"))
