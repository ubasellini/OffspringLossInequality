## --------------------------------------------------------- ##
##
## Routines to reproduce the results of the paper:
## "When do parents bury a child? Quantifying uncertainty
## in the parental age at offspring loss"
## 
## STEP 5: Process results to produce results reported in paper
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

## main analysis or sensitivity??
# type <- "main"
type <- "sensitivity"

## if sensitiviy analysis, take Swedish weights
if (type=="sensitivity"){
  cou <- "SWE"
  load(paste0("data/clean/",cou,"_arrays_clean.rdata"))
  ## cohorts
  cohorts <- as.numeric(colnames(D))[1]:2000
  ## weights
  wei.swe <- POP[,which(cohorts==1900)]/sum(POP[,which(cohorts==1900)])
}

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
  
  if (type=="main"){
  cb <- 
    read.csv("results/summary/counts_by_child_age.csv", stringsAsFactors = F) %>% 
    mutate(
      child_age = factor(child_age, levels = age_labs)
      , country = factor(country, levels = lookup_c[country_keep])
    )
  
  all_measures <- 
    read.csv("results/summary/summary_measures.csv", stringsAsFactors = F) %>% 
    mutate(mort = factor(mort, levels = ml_levs))
  }else{
    cb <- 
      read.csv("results/summ_sens/counts_by_child_age.csv", stringsAsFactors = F) %>% 
      mutate(
        child_age = factor(child_age, levels = age_labs)
        , country = factor(country, levels = lookup_c[country_keep])
      )
    
    all_measures <- 
      read.csv("results/summ_sens/summary_measures.csv", stringsAsFactors = F) %>% 
      mutate(mort = factor(mort, levels = ml_levs))
  }
} else{

  print("Running analysis scripts...")
  cat(paste(type,"analysis"),"\n")
  
  # 2.1. Non-decomposed mean and SD ===============
  
  if (type=="main"){
  non_decomp <- 
    lapply(countries, get_maol_all_batch) %>% 
    bind_rows()
  }else{
    non_decomp <- 
      lapply(countries, get_maol_all_batch, type="sensitivity",wei.swe=wei.swe) %>% 
      bind_rows()
  }
  # 2.2. Decomposed mean and SD ===============
  
  if (type=="main"){
  decomp <- 
    lapply(countries, get_curves) %>% 
    bind_rows() 
  }else{
    decomp <- 
      lapply(countries, get_curves, type="sensitivity") %>% 
      bind_rows() 
  }
  
  # 2.3. Get Prop offspring died ===========
  
  prop_died <- 
    lapply(countries, get_share_offspring_died, asfr_coh) %>% 
    bind_rows() 
  
  # 2.4. Combine and export ==========
  
  all_measures <- 
    # Summary measures
    decomp %>% 
    bind_rows(
      # Proportion children die before mother
      prop_died %>% 
        mutate(cohort = as.character(cohort))    
    ) %>% 
    # Non-decomposed measures
    bind_rows(
      non_decomp %>% 
        mutate(cohort = as.character(cohort)) 
    ) %>% 
    mutate(
      mort = factor(ml[country], levels = ml_levs)
      , measure2 = factor(sum_lookup[measure], levels = sum_lookup)    
      , country_full = lookup_c[country]
    )
  
  if (type=="main"){
    write.csv(all_measures, "results/summary/summary_measures.csv", row.names = F)
    print("Saved to results/summary/summary_measures.csv!")
  }else{
    write.csv(all_measures, "results/summ_sens/summary_measures.csv", row.names = F)
    print("Saved to results/summ_sens/summary_measures.csv!")
  }
  # 2.5. DF of age-composition of children that died ==============
  
  cba <-  
    lapply(country_keep, get_cd_by_age_batch) %>%
    # lapply(countes, get_cd_by_age_batch) %>% 
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
  
  write.csv(cb, "results/summary/counts_by_child_age.csv", row.names = F)
  
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
  geom_point(data = fix_axis, colour = "white") +
  scale_fill_brewer("Offspring's age at death", palette = "Dark2") +
  scale_x_continuous("Age of mothers born in 1950") +
  scale_y_continuous("Offspring deaths (count)", breaks = scales::pretty_breaks(n=5)) +
  facet_wrap(~country, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "bottom"
    , strip.background = element_rect(fill = "honeydew2", colour = "black")
  )

ggsave("figures/summary/!1_counts_by_child_age_death.pdf", width = 6, height = 3, units = "in")

# 2. Summary measures (nondecomposed) =================

# rcode kfns6

m_keep <- c("M", "SD", "w1", "M1", "M2", "SD1", "SD2")
pf <- round(seq(year_low, last_cohort, length.out = 6))

# To show historical data
rect_df <-
  data.frame(
    xmin = -Inf
    , xmax= 1950
    , ymin=-Inf
    , ymax=Inf
    , mort = factor("low", levels = ml_levs)
  )


all_measures %>% 
  filter(
    measure %in% c("M", "SD")
    , cohort >= year_low
  ) %>% 
  mutate(
    cohort = as.numeric(cohort)
    , measure2 = factor(sum_lookup[measure], levels = sum_lookup[m_keep])
    , l1 = str_extract(measure, "[A-Z]+")
    , l2 = paste0("Component ", str_extract(measure, "[0-9]+"))
    , l1 = gsub("M", "Mean", l1)
    , l1 = gsub("SD", "Standard Deviation", l1)
  ) %>% 
  ggplot() +
  # backgprund lines
  geom_line(
    aes(x = cohort, y = value, colour = mort, group = country)
    , data = . %>% filter(!country %in% country_keep) 
    , size = 0.3
    , alpha = 0.4
    , show.legend = F
  ) +  
  # Highlighted lines
  geom_line(
    aes(x = cohort, y = value, colour = mort, group = country)
    , data = . %>%   filter(country %in% country_keep)  
    , size = 1
  ) +
  # Country labels
  geom_label_repel(
    aes(x = cohort, y = value, colour = mort, label = lookup_c[country])
    , data = . %>% filter(
      country %in% country_keep
      , cohort == 1950
      , measure == "M"
    )
    , nudge_x = -2
    , size = 3
    , show.legend = F
  ) +
  # Rectangle to show hist data
  geom_rect(
    aes(xmin = xmin, xmax= 1950, ymin=-Inf, ymax=Inf)
    , color = NA
    , alpha = alpha_hist
    , data = rect_df
  )  +
  # Shape to identify countries
  geom_point(
    aes(x = cohort, y = value, colour = mort, shape = mort)
    , data = . %>% filter(
      country %in% country_keep
      , cohort %in% pf
    )
    , size = 3
  ) +
  facet_wrap(~l1 , scales = "free") +
  # scale_color_discrete("Mortality group") +
  scale_color_brewer("Mortality group", palette = pal) +
  scale_shape_discrete("Mortality group") +
  scale_x_continuous("Mother's birth cohort") +
  scale_y_continuous("Maternal age (years)") +
  theme_bw() +
  theme(
    legend.position = "bottom"
    , strip.background = element_rect(fill = "honeydew2", colour = "black")
  ) 

ggsave("figures/summary/!2_summary_nondecomp.pdf", width = 5, height = 3, units = "in")

# 3. Count of child deaths decomposed ============
# rcode ldshih


# A dummy df to make sure that y axes are the same in the plots
# showing the distributinos
fix_axis2 <-
  fix_axis %>% 
  rename(age = mother_age, country_full = country) %>% 
  mutate(measure = "GAMMA1.hat")


# Plot raw curves

lookup_comp <- c("Component 1: Young child deaths", "Component 2: Adult child deaths")
names(lookup_comp) <- c("GAMMA1.hat", "GAMMA2.hat")

dec <- 
  all_measures %>%
  filter(country %in% country_keep) %>% 
  mutate(
    cohort = as.numeric(cohort)
    , country_full = factor(country_full, levels = lookup_c[country_keep])
  ) %>% 
  filter(cohort >= 1950) 

# DF to add weights as text

w_xy <-
dec %>% 
  filter(cohort %in% 1950) %>% 
  filter(measure %in% c("GAMMA1.hat", "GAMMA2.hat")) %>% 
  select(age, value, country_full, measure) %>% 
  # Get x,y position for labels
  group_by(country_full, measure) %>% 
  arrange(desc(value)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  # Add weights
  left_join(
    dec %>% 
      filter(cohort %in% 1950) %>% 
      filter(measure %in% c("w1")) %>% 
      select(label = value, country_full)
    , by = c("country_full")
  ) %>% 
  mutate(
    label = round(ifelse(measure == "GAMMA1.hat", label, 1 - label), 2)
    , value = value * 0.3
    , age = age + 2
    , label = paste0(round(label*100, 0), "%")
    )

w_lab <- 
  data.frame(
    age = 40, value = 80, measure = "GAMMA1.hat", country_full = "Sweden"
    , label = "relative\nimportance"
    ) %>% 
  mutate(country_full = factor(country_full, levels = lookup_c[country_keep]))

w2 <-
  w_xy %>% 
  filter(country_full == "Sweden") %>% 
  select(xend = age, yend = value, measure, country_full) %>% 
  left_join(
    w_lab %>% 
      select(x = age, y = value, country_full)
    , by = "country_full"
  ) %>% 
  mutate(country_full = factor(country_full, levels = lookup_c[country_keep]))

  
# Fitted values

dec %>% 
  filter(cohort %in% 1950) %>% 
  ggplot(aes(x = age, y = value)) +
  # Fitted curves
  geom_area(
    aes(fill = measure)
    , size = 1
    , alpha = 0.2
    , position = 'identity'
    , data = .  %>% filter(measure %in% c("GAMMA1.hat", "GAMMA2.hat")) 
    , show.legend = F
  ) +
  # Fitted lines
  geom_line(
    aes(colour = measure)
    , size = 1
    , data = .  %>% filter(measure %in% c("GAMMA1.hat", "GAMMA2.hat")) 
    , show.legend = T
  ) +
  # Arrows
  geom_curve(
    aes(x = x, y = y, xend = xend, yend = yend, colour = measure)
    , data = w2
    , curvature = 0
  ) +
  # Weight label
  geom_label(
    aes(label = label)
    , label.size = NA
    , size = 3
    , data = w_lab
    , show.legend = F
  ) +
  # Text to identify components
  geom_label(
    aes(label = label, colour = measure)
    , data = w_xy
    , size = 3
    , show.legend = F
  ) +
  # Real values
  geom_point(
    data = .  %>% filter(measure == "count") 
    , show.legend = F
  ) +
  # To fix axis
  geom_point(data = fix_axis2, colour = "white") + 
  scale_color_brewer("", labels = lookup_comp, palette = pal) +
  scale_fill_brewer("", labels = lookup_comp, palette = pal) +
  scale_x_continuous("Age of mothers born in 1950") +
  scale_y_continuous("Offspring deaths (count)", breaks = scales::pretty_breaks(n=5)) +
  facet_wrap(~country_full, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "bottom"
    , strip.background = element_rect(fill = "honeydew2", colour = "black")
  )

if (type=="main"){
  ggsave("figures/summary/!3_counts_components.pdf", width = 6, height = 3, units = "in")
}else{
  ggsave("figures/summ_sens/!3_counts_components.pdf", width = 6, height = 3, units = "in")
}


# 4. Summary measures (decomposed) ========

# rcode oduha8ij

lookup_comp2 <- c("Component 1: Young child deaths", "Component 2: Adult child deaths")
names(lookup_comp2) <- c("Component 1", "Component 2")

m_keep2 <- c("M1", "M2", "SD1", "SD2")
pf <- round(seq(year_low, last_cohort, length.out = 6))

dec_df <-
  all_measures %>% 
  filter(measure %in% m_keep2) %>% 
  mutate(
    cohort = as.numeric(cohort)
    , measure2 = factor(sum_lookup[measure], levels = sum_lookup[m_keep2])
    , l1 = str_extract(measure, "[A-Z]+")
    , l2 = paste0("Component ", str_extract(measure, "[0-9]+"))
    , l1 = gsub("M", "Mean", l1)
    , l1 = gsub("SD", "Standard Deviation", l1)
  )

# To show historical data
rect_df <-
  data.frame(
    xmin = -Inf
    , xmax= 1950
    , ymin=-Inf
    , ymax=Inf
    , mort = factor("low", levels = ml_levs)
  )

p <- 
  dec_df %>% 
  mutate(l2 = lookup_comp2[l2]) %>% 
  filter(cohort >= year_low) %>% 
  ggplot() +
  # background lines
  geom_line(
    aes(x = cohort, y = value, colour = mort, group = country)
    , data = . %>% filter(!country %in% country_keep) 
    , size = 0.3
    , alpha = 0.4
    , show.legend = F
  ) +  
  # Highlighted lines
  geom_line(
    aes(x = cohort, y = value, colour = mort, group = country)
    , data = . %>%   filter(country %in% country_keep)  
    , size = 1
  ) +
  # Country labels
  geom_label_repel(
    aes(x = cohort, y = value, colour = mort, label = lookup_c[country])
    , data = . %>% filter(
      country %in% country_keep
      , cohort == 1950
      , measure == "SD1"
    )
    , nudge_x = -2
    , size = 3
    , show.legend = F
  ) +
  # Rectangle to show hist data
  geom_rect(
    aes(xmin = xmin, xmax= 1950, ymin=-Inf, ymax=Inf)
    , color = NA
    , alpha = alpha_hist
    , data = rect_df
  )  +
  # Shape to identify countries
  geom_point(
    aes(x = cohort, y = value, colour = mort, shape = mort)
    , data = . %>% filter(
      country %in% country_keep
      , cohort %in% pf
    )
    , size = 3
  ) +
  facet_grid(l1~ l2 , scales = "free") +
  # scale_color_discrete("Mortality group") +
  scale_color_brewer("Mortality group", palette = pal) +
  scale_shape_discrete("Mortality group") +
  scale_x_continuous("Mother's birth cohort") +
  scale_y_continuous("Maternal age (years)") +
  theme_bw() +
  theme(
    legend.position = "bottom"
    , plot.margin = margin(0,0, 0, 0, "cm")
    , strip.background = element_rect(fill = "honeydew2", colour = "black")
  ) 

# DF of weights

w_df <- 
  all_measures %>% 
  filter(measure == "w1") %>%
  filter(cohort >= year_low) %>% 
  mutate(
    cohort = as.numeric(cohort)
    , l2 = "Component 1"
  ) %>% 
  select(country, cohort, mort, value, l2)

w_df <- 
  bind_rows(
    w_df
    , w_df %>% 
      mutate(
        value = 1 - value
        , l2 = "Component 2"
      )
  ) %>% 
  mutate(l1 = "Weight") %>% 
  mutate(l2 = lookup_comp2[l2])

p_w <- 
  w_df %>%
  ggplot() +
  # backgprund lines
  geom_line(
    aes(x = cohort, y = value, colour = mort, group = country)
    , data = . %>% filter(!country %in% country_keep) 
    , size = 0.3
    , alpha = 0.4
    , show.legend = F
  ) +  
  # HIghlighted lines
  geom_line(
    aes(x = cohort, y = value, colour = mort)
    , data =  . %>%  
      filter(country %in% country_keep) 
    
    , size = 1
    , show.legend = F
  ) +
  # Shape to identify countries
  geom_point(
    aes(x = cohort, y = value, colour = mort, shape = mort)
    , data = . %>% filter(
      country %in% country_keep
      , cohort %in% pf
    )
    , size = 3
    , show.legend = F
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_color_brewer("", palette = pal) +
  facet_grid(. ~ l2 , scales = "fixed") +
  labs(x = "", y = "Relative\nimportance") +
  coord_cartesian(ylim = c(0,1))+
  theme_minimal() +
  theme(
    legend.position = "bottom"
    , axis.title.x=element_blank()
    , axis.text.x=element_blank()
    , axis.ticks.x=element_blank()
    , strip.text.x = element_blank()
    , plot.margin = margin(0,0, 0, 0, "cm")
  ) 

p_w / p +
  plot_layout(widths = c(2, 1), heights = unit(c(1.75, 5), c('cm', 'null')))

if (type=="main"){
  ggsave("figures/summary/!4_summary_decomp.pdf", width = 5, height = 5, units = "in")  
}else{
  ggsave("figures/summ_sens/!4_summary_decomp.pdf", width = 5, height = 5, units = "in")  
}

# 5. Share of children outlive mothers =========
# rcode kjgw76

all_measures %>% 
  filter(
    measure %in% "prop_died"
    , cohort >= year_low
  ) %>% 
  mutate(
    cohort = as.numeric(cohort)
    , l2 = paste0("Component ", str_extract(measure, "[0-9]+"))
  ) %>% 
  ggplot() +
  # background lines
  geom_line(
    aes(x = cohort, y = value, colour = mort, group = country)
    , data = . %>% filter(!country %in% country_keep) 
    , size = 0.3
    , alpha = 0.4
    , show.legend = F
  ) +  
  # Highlighted lines
  geom_line(
    aes(x = cohort, y = value, colour = mort, group = country)
    , data = . %>%   filter(country %in% country_keep)  
    , size = 1
  ) +
  # Country labels
  geom_label_repel(
    aes(x = cohort, y = value, colour = mort, label = lookup_c[country])
    , data = . %>% filter(
      country %in% country_keep
      , cohort == 1950
    )
    , nudge_x = 1
    , nudge_y = 0.005
    , size = 3
    , show.legend = F
  ) +
  # REct to show hist data
  geom_rect(
    aes(xmin = xmin, xmax= 1950, ymin=-Inf, ymax=Inf)
    , color = NA
    , alpha = alpha_hist
    , data = rect_df
  )  +
  # Shape to identify country groups
  geom_point(
    aes(x = cohort, y = value, colour = mort, shape = mort)
    , data = . %>% filter(
      country %in% country_keep
      , cohort %in% pf
    )
    , size = 3
  ) +
  scale_color_brewer("Mortality group", palette = pal) +
  scale_shape_discrete("Mortality group") +
  scale_x_continuous("Mother's birth cohort") +
  scale_y_continuous("Offspring die before mother", labels = function(x) paste0(round(x*100), "%")) +
  theme_bw() +
  theme(
    legend.position = "bottom"
    , strip.background = element_rect(fill = "honeydew2", colour = "black")
  ) 

ggsave("figures/summary/!5_share_child_die.pdf", width = 4, height = 3, units = "in")