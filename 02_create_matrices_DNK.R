# Note: wd should be one level up from the R dir - where the Rproj is stored

# This script does the following using data for Denmark combining HMD, HFD, and WPP:

# A. Load all data 
# B. Approximates cohort data from period
# C. Estimate matrices for analysis

rm(list=ls())

library(tidyverse)
source("functions/DiegoFUNS.R")

# 0. Parameters ~~~~~~~~~ ----------- 

# Data for Denmark: 
# HMD: 1835-2020
# HFD: 1916-2020
# HFD: 1878-2020


# This is the ages at death for children, seen from the 
# perspective of mothers. 
ages_mother <- 0:100
reprod_ages <- 15:50
ages_child <- 0:(max(ages_mother) - min(reprod_ages))

years_unwpp <- 1878:2100
cohorts_unwpp <- 1878:2050

# For running kins function
# For which cohorts should we get the entire kinship structure
# using period rates as input?
cohorts_demokin <- cohorts_unwpp

# Which maternal cohorts do you want to get the child loss rates for?
cohorts_kin <- 1878:2000

# To rescale life ables from hmd
max.age <- 100

# To read the demokin output files
base <- "data/demokin_denmark/"
pat <- "tree"

# Make sure that directory exists
if(!dir.exists(base)) dir.create(base)

# For identifying countries
lookup_c <- c("Denmark")
names(lookup_c) <- 1

lookup_c2 <- names(lookup_c)
names(lookup_c2) <-  lookup_c

# Parameters to expand period life tables into the future
# assuming stability at the last observed value
# Thhis is needed to account for children born in year 2020 
# dying at age 85+, for example, since WPPP projections just 
# get to year 2100
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

# This function loads period fertility and mortality data from 
# various sources and processes it. This includes (whenever applicable):
# - Renaming columns and changing radices of different sources to 
# make sure that they can be combined into a single object
# - 'Unpacking' 5-y age or calendar years measures into single-age
# and single year estimates
# - Changing the open-age interval of life tables and population totals to be 100
# - Converting period to cohort values by taking the values along the diagonal
# - Expand cohort life tables to the future (assuming demographic stability 
# after the 2000 cohort). This is needed because some offspring may be born 
# (and thus die) after 2000, which is the last LT we have (the last offspring)
# is born in cohort 2049 for woman born in year 2000.

rate_dat <- 
  get_denmark_hmd_hfd_data_clean2(
    max.age = max.age
    , years_unwpp = years_unwpp
    , cohorts_unwpp = cohorts_unwpp
    , expand_lt_to_future_stable = T
    , export = T
  )

lt_coh <- rate_dat[["lt_coh"]] 
lt_1_1 <- rate_dat[["lt_1_1"]] 
asfr_1_1 <- rate_dat[["asfr_1_1"]] 
pop_1_1 <- rate_dat[["pop_1_1"]] 

# C. Kinship matrices ~~~~~~~~~ ----------- 

# In this section, I take the UNWPP rates and run the matrix kinship models to 
# get arrays of the expected number of children and mothers by birth cohort and 
# ages in non-stable population. 
# I recommend checking the Methods sectino of our paper for more details. 
# Note that I use an optimised version of `DemoKin::kins()` that considers 
# only mothers and daughters (See the vignette of DemoKin and 
# Caswell 2019 for details). 
# This function saves the results to the disk and only needs to be run once!
# Takes about 5 sec per year

# Some parameters for this chunk
age <- max_age <- max.age
cohorts <- cohorts_demokin

# Iterate over all countries 
con <- 1
# cohort <- cohorts[1]
for(cohort in cohorts){
  print(paste0("Working on country: ", lookup_c, " - ",cohort))
  print(Sys.time())
  
  out_file <- paste0(base, pat, "_", lookup_c,"_", cohort, ".rds")
  
  # Get rates to input in DemoKin function
  if(file.exists(out_file)){
    print(paste0(out_file, " already exists. I'll just start with the next one."))
  } else {
    ages_v <- 0:max_age
    ages <- length(ages_v)
    w <- max_age
    
    L <- 
      lt_1_1 %>% 
      filter(age <= w) %>% 
      mutate(nLx = ifelse(age == w, Tx, nLx)) %>% 
      select(Year, age, nLx) %>% 
      spread(Year, nLx) %>% 
      select(-age)
    
    age_plus1 <- max_age + 1
    
    # The argument P in DemoKin::kins() is S(x)=L(x+1)/L(x)
    # The first row should be L(1)/L(0) 
    
    P <- 
      rbind(
        L[c(-1, -age_plus1), ]/L[-c(max_age:age_plus1), ]
        , L[age_plus1, ]/(L[max_age, ] + L[age_plus1, ])
        , L[age_plus1, ]/(L[max_age, ] + L[age_plus1, ])
      )
    
    rownames(P) <- ages_v
    P[is.na(P)] <- 0
    
    asfr_mat <- 
      asfr_1_1 %>% 
      select(Year, age, value = asfr) %>% 
      mutate(
        age = factor(age-1, levels = ages_v)
        , Year = factor(Year, levels = years_unwpp)
      ) %>% 
      tidyr::complete(Year, age, fill = list(value = 0)) %>% 
      pivot_wider(names_from = Year, values_from = value) %>% 
      select(-age) %>% 
      data.frame()
    
    colnames(asfr_mat) <- years_unwpp
    
    pop_mat <- 
      pop_1_1 %>% 
      select(Year, age, value = PopFemale) %>% 
      pivot_wider(names_from = Year, values_from = value) %>%
      select(-age) %>% 
      data.frame()
    
    colnames(pop_mat) <- years_unwpp
    
    # Actually implement the function to get the matrices of
    # kin availability for mothers and daughters
    tree <- 
      kins_fast(
        ego_age = age
        , year = cohort + age
        , P = P
        , asfr = asfr_mat
        , N = pop_mat
        , stable = FALSE
      )
    
    saveRDS(tree, out_file)        
  } # end get rates
}  # end cohort loop

# D. Estimate matrices for analysis ~~~~~~~~~ ----------- 

# Here, I apply the time-variant method described in the paper's Methods 
# section to get an array of the expected number of surviving children by child's
# age, maternal age, and maternal cohort. 
 
print("Reformating DemoKin outputs into arrays...")

# Get arrays for children
# The idea is to have the rows represent child's age
# columns represent mother's age, and the third dimension
# be mother's birth cohort
# The function `reformat_kin_trees` transforms this into a three-dimensional 
# array, where the third dimension is cohort.
arys <- reformat_kin_trees2(dataset = "kins", relatives = c("d"))

# 'd' stands for daughters in DemoKin
H <- arys[["d"]]

# Let's assign the cohorts as names to this array
dimnames(H)[[3]] <-
  (1:dim(H)[3]) + min(cohorts_unwpp) - 1

# View(H[,,1])

# Create matrices for the life-table objects that we need:

lt  <-
  lt_coh %>%
  group_by(LocID, cohort) %>%
  # lx should be normalised to be a density function
  mutate(lx = lx/first(lx)) %>%
  ungroup() %>%
  data.frame()

L <- 
  lt %>% 
  select(cohort, age, val = lx) %>% 
  pivot_wider(names_from = "cohort", values_from = "val") %>% 
  select(-age) %>% 
  as.matrix()

D <- 
  lt %>% 
  select(cohort, age, val = ndx) %>% 
  pivot_wider(names_from = "cohort", values_from = "val") %>% 
  select(-age) %>% 
  as.matrix()

Q <- 
  lt %>% 
  select(cohort, age, val = nqx) %>% 
  pivot_wider(names_from = "cohort", values_from = "val") %>% 
  select(-age) %>% 
  as.matrix()

# Let's assign the cohorts as names to this array

colnames(L) <- 
  colnames(D) <- 
  colnames(Q) <- 
  (1:dim(L)[2]) + min(cohorts_unwpp) - 1

# Matrix of Population ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# gET LX COLUMNS OF A LIFE TABLE GIVEN A RADIX
# AND VECTOR OF MORT PROBS
get_lx <- function(l0, nqx){ 
  lx <- cumprod(c(l0, 1 - nqx)) 
  lx[-length(lx)]
}

l0_df <- 
  pop_1_1 %>% 
  select(year = Year, age, value = PopFemale) %>% 
  filter(age == 0) %>% 
  filter(year %in% cohorts_unwpp)

l0 <- l0_df$value
names(l0) <- l0_df$year

POP <- mapply(
  get_lx
  , l0
  , as.data.frame(Q)
  , SIMPLIFY = T
)

dim(POP)

# Export arrays -----
# rcode lksdn437

save(H, D, L, Q, POP, file = "data/clean/DNK_arrays_clean.rdata")
