## --------------------------------------------------------- ##
##
## Routines to reproduce the results of the paper:
## "When do parents bury a child? Quantifying uncertainty
## in the parental age at offspring loss"
## 
## STEP 3: Create matrix objects for 16 countries
##
## Code by Diego Alburez-Gutierrez (2022) unless otherwise stated
##
##  sessionInfo() details:
##
##  R version 4.0.2 (2020-06-22)
##  Platform: x86_64-w64-mingw32/x64 (64-bit)
##  Running under: Windows 10 x64 (build 19044)
##
##  attached base packages:
##  stats   graphics  grDevices utils   datasets 
##  methods   base
## 
##  other attached packages:
##  quadprog_1.5-8        httr_1.4.2            countrycode_1.2.0 
##  data.table_1.13.2     forcats_0.5.1         stringr_1.4.0        
##  dplyr_1.0.7           purrr_0.3.4           readr_1.4.0          
##  tidyr_1.1.4           tibble_3.1.5          ggplot2_3.3.5        
##  tidyverse_1.3.1
## --------------------------------------------------------- ##

# This script does the following using data for 16 countries using WPP data:

# A. Load all data 
# B. Approximates cohort data from period
# C. Estimate matrices for analysis

rm(list=ls())

library(tidyverse)
library(data.table)
library(countrycode)

# These two are only needed if reestimate = T inthe function 
# get_unwpp_cohort_rates() below
library(httr)
library(quadprog)

source("functions/analysisFUNS.R")

# 0. Parameters --------
# Countries to keep
country_key <- 
  read.csv("data/input/unwpp/country_key.csv", stringsAsFactors = F, encoding = "UTF-8") %>% 
  filter(!is.na(new_name))

lookup_c <-  c(
  "United States of America"
  # These two are estimated in separate scripts because they use data from 
  # elsewhere
  # ,  "Denmark" , "Sweden"
  , "China", "Guatemala", "Algeria"
  , "Mali", "Chad", "Burkina Faso", "Senegal", "Ghana", "El Salvador"
  , "Australia", "Japan",
  "Cuba", "India", "Myanmar", "Angola"
  )

# Create country lookup tables
names(lookup_c) <- country_key$LocID[match(lookup_c, country_key$new_name)]
lookup_c2 <- names(lookup_c)
names(lookup_c2) <-  lookup_c

# PARAMETERS FOR KINSHIP MODELS
# Estimates are stored in disk
run_matrix_models <- T

parallelise_kin_models <- F

# cores <- detectCores() - 3
print(paste0("Parallelise kinship models: ", parallelise_kin_models))

# This is the ages at death for children, seen from the 
# perspective of mothers. 
ages_mother <- 0:100
reprod_ages <- 15:50
ages_child <- 0:(max(ages_mother) - min(reprod_ages))

years_unwpp <- 1950:2100
cohorts_unwpp <- 1950:2000
# For running kins function
# For which cohorts should we get the entire kinship structure
# using period rates as input?
cohorts_demokin <- 1950:2050

# Which maternal cohorts do you want to get the child loss rates for?
# Requires some thinking, but ultimately this shoudl be true:
cohorts_kin <- cohorts_unwpp

# Smaller selection of cohorts for which rates will be actually estimated
# If you don't need estimates for all cohorts, chose 
# how many yeares you want in between cohort.  
by_cohs <- 5
coh_keep <- cohorts_kin[seq(1, length(cohorts_kin), by = by_cohs)]

# To read the demokin output files
base <- "data/demokin_unwpp/"
pat <- "tree"

# Make sure that directory exists
if(!dir.exists(base)) dir.create(base)

# Parameters to expand period life tables into the future
# assuming stability at the last observed value
expand_lt_to_future_stable <- T
years_to_expand <- 50
up_to_year <- max(years_unwpp) + years_to_expand

cohort_lt <- 
  switch(
    expand_lt_to_future_stable
    , min(cohorts_unwpp):(max(cohorts_unwpp)+years_to_expand)
    , cohorts_unwpp
  )

# A. Load all data ~~~~~~~~~ ----------- 

# We basically want the lx, LX, ex, asfr, and TFR for each country, 
# historical and projected data.

rates_l <- 
  get_unwpp_cohort_rates(
    lookup_c
    # reestimate = T downloads the raw data from the UNWPP
    # website and does the period-cohort conversion
    # at the end, it saves the products as csv files
    # reestimate = F assumes that these products already exist
    # The function is smart, so if the raw unwpp files are not available
    # it will download them.
    # So, the first time you run the function, it's better to run with 
    # reestimate = T
    , reestimate = T
    , expand_lt_to_future_stable
    , up_to_year = up_to_year
    , cohorts_unwpp = cohorts_unwpp
    , cohort_lt = cohort_lt
    , country_key
    , export = T
  )

# Cohort life tables
lt_coh <- rates_l[['lt_coh']]
# Period life tables
lt_1_1 <- rates_l[['lt_1_1']]
lt_1_1 <- lt_1_1 %>% filter(Year <= 2100)
# Period age-specific fertility rates
asfr_1_1 <- rates_l[['asfr_1_1']]
# Yearly population figures
pop_1_1 <- rates_l[['pop_1_1']]

rm("rates_l")

# B. Kinship matrices ~~~~~~~~~ ----------- 

# This runs DemoKin to get the matrices of maternal and 
# offspring availability in non-stable populations.

# It runs the model for all country/year combinations, so it
# may take a while to run...


if(run_matrix_models){
  get_kin_trees_unwpp(
    asfr_1_1 = asfr_1_1
    , lt_1_1 = lt_1_1
    , pop_1_1 = pop_1_1
    , cohorts = cohorts_demokin
    , age = 100
    , max_age = max(ages_mother)
    , years_unwpp = years_unwpp
    , base = base
    , pat = pat
    , parallel = parallelise_kin_models
    , cores = cores
  )  
}

# C. Estimate matrices for analysis ~~~~~~~~~ ----------- 
# rcode ks4dgf

# Here, I wrapped everything up in a function, but, essentially, we are doing
# the same we did for Sweden and Denmark.
# This function computes the different matrices needed for the analysis and
# save the to the disk, one for each country. 

for(con in lookup_c){
  get_matrices_wpp(country_keep = con, return_list = F) 
}
