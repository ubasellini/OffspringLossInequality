## --------------------------------------------------------- ##
##
## Routines to reproduce the results of the paper:
## "When do parents bury a child? Quantifying uncertainty
## in the parental age at offspring loss"
## 
## STEP 6: Counterfactual analysis
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
##  ggrepel_0.8.2         countrycode_1.2.0     patchwork_1.1.1
##  forcats_0.5.1         stringr_1.4.0        
##  dplyr_1.0.7           purrr_0.3.4           readr_1.4.0          
##  tidyr_1.1.4           tibble_3.1.5          ggplot2_3.3.5        
##  tidyverse_1.3.1
## --------------------------------------------------------- ##

## cleaning the workspace
rm(list=ls())

# 0. Parameters ---------------

## packages
library(tidyverse) 
library(patchwork)
library(countrycode)
library(ggrepel)
source("functions/analysisFUNS.R")


## loading saved data
countries <- str_extract(list.files("data/clean", pattern = ".rdata$"), "^[A-Z]{3}")

last_cohort <- 2000

# Group countries by mortality levels
# I did this in a different script 
# rcode sljngs
ml <- c(MLI = "high", AGO = "high", TCD = "high", BFA = "high", SEN = "high", 
        GHA = "high", MMR = "medium", IND = "medium", SLV = "medium", 
        GTM = "medium", DZA = "medium", CHN = "medium", CUB = "low", 
        USA = "low", DNK = "low", AUS = "low", SWE = "low", JPN = "low"
)

ml_levs <- c("low", "medium", "high")

lookup_c <- countrycode(names(ml), origin = "iso3c", destination = "country.name")
lookup_c[lookup_c == "Myanmar (Burma)"] <- "Myanmar"
names(lookup_c) <- names(ml)

# For smaller plots
country_keep <- c("SWE", "GTM", "AGO")

sum_measures <- 
  c(
    "prop_died"
    , "M1"
    , "M2"
    , "M"
    , "w1"
    , "SD1"
    ,  "SD2"
    , "SD"
  )

sum_lookup <- c(
  "Prop offspring will die before mother"
  , "Mean age at loss (comp1)"
  , "Mean age at loss (comp2)"
  , "Mean age at loss"
  , "Weight of young-child deaths"
  , "SD age at loss (comp1)"
  , "SD age at loss (comp2)"
  , "SD age at loss"
)

names(sum_lookup) <- sum_measures

# Parameters for plotting

year_low <- 1850

pal <- "Dark2"

# Background of shape for historical area
alpha_hist <- 0.07

# A. Analysis +++++++++++++++++ ---------------

# 1. Create data frame of cohort ASFR =============

# Consolidate df of 1x1 cohort ASFR for all countries into single object
# This will be used in analysis below

if(!file.exists("data/input/asfr_coh.csv")){
  cn <- read.csv("data/input/unwpp/country_key.csv", stringsAsFactors = F)
  
  lookup_locid <- cn$iso
  names(lookup_locid) <- as.character(cn$LocID)
  
  asfr_1_1 <- 
    bind_rows(
      read.csv(paste0("data/input/asfr_1_1_SWE.csv"), stringsAsFactors = T) %>% 
        mutate(country = "SWE")
      , read.csv(paste0("data/input/asfr_1_1_DNK.csv"), stringsAsFactors = T) %>% 
        mutate(country = "DNK")
      , read.csv(paste0("data/input/asfr_1_1_unwpp.csv"), stringsAsFactors = T) %>% 
        mutate(country = lookup_locid[as.character(LocID)])
    ) %>% 
    select(-LocID) %>% 
    rename(LocID = country)
  
  asfr_coh <- 
    asfr_1_1 %>% 
    group_by(LocID, Year) %>%
    group_map( 
      ~ approximate_cohort_from_period(., add_all_combinations = F)
      , .keep = T
    ) %>%
    bind_rows() %>% 
    ungroup() %>% 
    select(-Year) %>% 
    rename(country = LocID)
  
  write.csv(asfr_coh, "data/input/asfr_coh.csv", row.names = F)  
} else{
  asfr_coh <- read.csv("data/input/asfr_coh.csv", stringsAsFactors = F)
}

# 2. Run analysis ----------

# The first time you run the scripts, you should use
# recompute <- T
# If, for any reason you want to run the scripts again but don't want to
# run the lengthy analysis scripts again, you can set recompute <- F
recompute <- T

if(!recompute){
  print("Reading analysis outcomes from disk...")
  
  # Just read the values without recomputing the measures
  
  age_labs <- c("[0,4]", "(4,25]", "(25,50]", "(50,75]", "(75,100]")
  
  cb <- 
    read.csv("results/summary/counts_by_child_age_counterfactual.csv", stringsAsFactors = F) %>% 
    mutate(
      child_age = factor(child_age, levels = age_labs)
      , country = factor(country, levels = lookup_c[country_keep])
    )
  
} else {
  
  print("Running analysis scripts...")
  
  # 2.5. Counterfactual ==============
  
  cba <-  
    lapply(country_keep, get_cd_by_age_batch_counterfactual) %>%
    bind_rows()
  
  child_age_br <- c(0,4, seq(25, 100, 25))
  
  cb <- 
    cba %>% 
    filter(country %in% country_keep) %>% 
    filter(mother_age >= 15) %>% 
    mutate(
      child_age2 = cut(child_age, child_age_br, include.lowest = T, right = T)
      , country = factor(lookup_c[country], levels = lookup_c[country_keep])
    ) %>% 
    group_by(country, cohort, mother_age, child_age = child_age2) %>% 
    summarise(value = sum(value, na.rm = T)) %>% 
    ungroup()
  
  # dput(levels(cb$child_age))
  
  write.csv(cb, "results/summary/counts_by_child_age_counterfactual.csv", row.names = F)
  
}

# B. Plots +++++++++++++++++ ----------

# rcode lknd83
# 1. Count of child deaths by child's age at death =========

# A dummy df to make sure that y axes are the same in the plots
# showing the distributinos
fix_axis <-
  data.frame(
    country = factor(lookup_c[country_keep], levels = lookup_c[country_keep])
    , mother_age = 50
    , value = c(100, 1000, 2100)
  )

# Line of overall number of deaths

cb2 <- 
  cb %>% 
  filter(cohort %in% c(1950))

counts <- 
  cb2 %>%
  group_by(country, cohort, mother_age) %>% 
  summarise(value = sum(value)) %>% 
  ungroup()

cb2 %>% 
  ggplot(aes(x = mother_age, y = value)) +
  geom_area(aes(fill = child_age), colour = "black") +
  geom_point(data = counts) +
  scale_fill_brewer("Offspring's age at death", palette = "Dark2") +
  scale_x_continuous("Age of mothers born in 1950") +
  scale_y_continuous("Offspring deaths (count)", breaks = scales::pretty_breaks(n=5)) +
  facet_wrap(~country, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "bottom"
    , strip.background = element_rect(fill = "honeydew2", colour = "black")
  )

ggsave("figures/summary/!6_counterfactual.pdf", width = 6, height = 3, units = "in")
