# Set up
setwd("/home/r044073/Chemx_IQ/")
.libPaths("/home/r044073/R/x86_64-pc-linux-gnu-library/3.6/")

time1 <- Sys.time()

# Load packages
library(mice)
library(tidyverse)
library(bkmr)
library(bkmrhat) # parallel running extension of BKMR

# load manipulated data----------------------------------------------------------

load("/home/r044073/Chemx_IQ/data/DATA_afimp_manipulated_2025-07-23.RData")

# set date
systime <- Sys.Date()

analysis = "TSH mother"
print(analysis)

#####################################
##                                 ##
##            BKMR-MI              ##-------------------------------------------
##                                 ##
#####################################


#################################################################################
# Manipulate imputed data for analysis 

# Covariates - mothers
cov <- c(
  "gestage_serum_thyroid_g1", "AGE_M_v2", "BMI_0", "PARITY",
  "GENDER", "edu", "ethm", 
  "SMOKE_ALL",
  "TPOAb"
)

dd_m0 <- df_afimp %>% complete("long", include = T) %>% 
  filter(tdis == 0 & med == 0 & ivf == 0) %>% 
  filter(mis_chemx_preg_thyroidm == 1 & mis_tsh_preg == 1 & tsh_mother_outlier == 0) %>% 
  select(.imp, .id,
         # BzBP, high molecular weight
         mBzBPngmL_AVGg1g2_imputed_bratio_log2, 
         # DEHP, high molecular weight
         mCMHPngmL_AVGg1g2_imputed_bratio_log2, mECPPngmL_AVGg1g2_imputed_bratio_log2, 
         mEHHPngmL_AVGg1g2_imputed_bratio_log2, mEOHPngmL_AVGg1g2_imputed_bratio_log2, 
         # DBP, low molecular weight
         mBPngmL_AVGg1g2_imputed_bratio_log2, mIBPngmL_AVGg1g2_imputed_bratio_log2,
         # DEP, low molecular weight
         mEPngmL_AVGg1g2_imputed_bratio_log2, 
         # DMP, low molecular weight
         mMPngmL_AVGg1g2_imputed_bratio_log2, mCPPngmL_AVGg1g2_imputed_bratio_log2,  
         # PA
         #PAngmL_AVGg1g2_imputed_bratio_log2,
         # BPA
         BPAngmL_AVGg1g2_imputed_bratio_log2, 
         # DE
         DEP_AVGg1g2_imputed_bratio_log2, DETP_AVGg1g2_imputed_bratio_log2, 
         # DM
         DMDTP_AVGg1g2_imputed_bratio_log2, DMP_AVGg1g2_imputed_bratio_log2, DMTP_AVGg1g2_imputed_bratio_log2,
         # outcome
         TSH_mother_log2, TT4_mother_log2, FT4_eci_mother_log2, 
         TSH_child_log2, FT4_eci_child_log2,
         TSHchildF5_log2, FT4childF5_log2, 
         all_of(cov)
  ) %>% 
  rename_with(~ str_remove(., "_AVGg1g2_imputed_bratio_log2"), contains("_AVGg1g2_imputed_bratio_log2")) %>% 
  rename_with(~ str_remove(., "ngmL"), contains("ngmL")) %>% 
  mutate(across(mBzBP:DMTP, ~ pmin(.x, quantile(.x, 0.99)))) %>% 
  mutate(across(all_of(cov), as.numeric))

dd_m0 <- as.mids(dd_m0)

# Create the 30 imputed datasets
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
      mBzBP, 
      mCMHP, mECPP, 
      mEHHP, mEOHP, 
      mBP, mIBP,
      mEP, 
      mMP, mCPP,  
      BPA, 
      DEP, DETP, 
      DMDTP, DMP, DMTP
    ) %>% 
    as.matrix()
  
  # Z-score exposure matrix to put on same scale 
  Z_preg_scale <- scale(Z_preg)
  
  # covariates need to be numeric
  X_preg <- dd_m %>% select(all_of(cov)) %>% 
    mutate(across(where(is.numeric), scale)) %>% 
    as.matrix()
  
  # Outcome
  outcome <- dd_m$TSH_mother_log2
  
  # Set random seed for reproducibility
  set.seed(2024)
  
  future::plan(strategy = future::multisession)
  
  # Run BKMR model
  fit <- bkmrhat::kmbayes_parallel(outcome, Z_preg_scale, X_preg, nchains = 8, iter = 6250, 
                                   group = c(1,2,2,2,2,3,3,4,5,6,7,8,8,8,8,8),
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
save(BKMRfits, file = paste0("/home/r044073/Chemx_IQ/results/thyroid/BKMR_MI_EDC_TSHm_", systime, ".RData"))

