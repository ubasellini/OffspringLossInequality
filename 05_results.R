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
library(viridis)
source("functions/analysisFUNS.R")

windowsFonts(Times=windowsFont("Times New Roman"))


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
recompute <- F

if(!recompute){
  print("Reading analysis outcomes from disk...")

  # Just read the values without recomputing the measures
  
  age_labs <- c("[0,4]", "(4,25]", "(25,50]", "(50,75]", "(75,100]")
  
  cb <- 
    read.csv("results/summary/counts_by_child_age.csv", stringsAsFactors = F) %>% 
    mutate(
      child_age = factor(child_age, levels = age_labs)
      , country = factor(country, levels = lookup_c[country_keep])
    )
  
  all_measures <- 
    read.csv("results/summary/summary_measures.csv", stringsAsFactors = F) %>% 
    mutate(mort = factor(mort, levels = ml_levs))
  
} else {

  print("Running analysis scripts...")
  # 2.1. Non-decomposed mean and SD ===============
  
  non_decomp <- 
    lapply(countries, get_maol_all_batch) %>% 
    bind_rows()
  
  # 2.2. Decomposed mean and SD ===============
  
  decomp <- 
    lapply(countries, get_curves) %>% 
    bind_rows() 
  
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
  
  write.csv(all_measures, "results/summary/summary_measures.csv", row.names = F)
  print("Saved to results/summary/summary_measures.csv!")
  
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
    , name = "Count"
  )

# Line of overall number of deaths

cb2 <- 
  cb %>% 
  filter(cohort %in% c(1950))

counts <- 
  cb2 %>%
  group_by(country, cohort, mother_age) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  mutate(name = "Count")

# Add two panels below with deaths per 1,000 women
# First, get pop data per birth cohort
# Sweden
load("data/clean/SWE_arrays_clean.rdata")

pop1950coh_swe <- POP[,"1950"]
rm("H", 'D', "L", "Q", "POP")

# Guatemala
load("data/clean/GTM_arrays_clean.rdata")

pop1950coh_gtm <- POP[,"1950"]
rm("H", 'D', "L", "Q", "POP")

# Guatemala
load("data/clean/AGO_arrays_clean.rdata")

pop1950coh_ago <- POP[,"1950"]
rm("H", 'D', "L", "Q", "POP")

# Put together
pop_coh <- 
  data.frame(
    Sweden = pop1950coh_swe
    , Guatemala = pop1950coh_gtm
    , Angola = pop1950coh_ago
    , mother_age = 0:100
  ) %>% 
  pivot_longer(-mother_age, names_to = "country", values_to = "pop") %>% 
  mutate(cohort = 1950)

cb3 <- 
  cb2 %>% 
  left_join(pop_coh, by = c("country", "cohort", "mother_age")) %>% 
  mutate(`Rate per 1000 women` = (value / pop) * 1000) %>% 
  select(-pop) %>% 
  rename(Count = value) %>% 
  pivot_longer(Count:`Rate per 1000 women`) %>% 
  mutate(country = factor(country, levels = c("Sweden", "Guatemala", "Angola")))

lab_df <-
  cb3 %>% 
  group_by(country, name, mother_age) %>% 
  summarise(value = sum(value)) %>% 
  arrange(desc(value)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  arrange(name, country) %>% 
  mutate(
    mother_age = 50
    # , label = LETTERS[1:6]
    , label = paste0("(",letters[1:6],")")
    )

cb3 %>% 
  ggplot(aes(x = mother_age, y = value)) +
  geom_area(aes(fill = child_age), colour = "black") +
  geom_point(data = counts) +
  geom_point(data = fix_axis, colour = "white") +
  # Panel labels
  geom_text(
    aes(label = label)
    , family = "Times"
    , data = lab_df 
  ) +
  # scale_fill_brewer("Offspring's age at death", palette = "Dark2") +
  scale_fill_viridis_d("Offspring's age at death") +
  scale_x_continuous("Age of mothers born in 1950") +
  scale_y_continuous(
    "Offspring deaths", breaks = scales::pretty_breaks(n=5)
    , labels = function(x) format(x, big.mark = ",")
    ) +
  facet_wrap(name~country, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "bottom"
    , strip.background = element_rect(fill = "honeydew2", colour = "black")
    , text=element_text(family="Times")
  )

# ggsave("figures/summary/!1_counts_by_child_age_death.pdf", width = 6, height = 6, units = "in")
ggsave("figures/summary/Fig3.pdf", width = 6, height = 6, units = "in")

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
    , l1 = gsub("SD", "Standard deviation", l1)
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
    , family = "Times"
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
  scale_y_continuous("Maternal age at offspring loss (years)") +
  theme_bw() +
  theme(
    legend.position = "bottom"
    , strip.background = element_rect(fill = "honeydew2", colour = "black")
    , text=element_text(family="Times")
  ) 

# ggsave("figures/summary/!2_summary_nondecomp.pdf", width = 5, height = 3.5, units = "in")
ggsave("figures/summary/Fig4.pdf", width = 5, height = 3.5, units = "in")

# 3. Count of child deaths decomposed ============
# rcode ldshih

# A dummy df to make sure that y axes are the same in the plots
# showing the distributinos
fix_axis2 <-
  fix_axis %>% 
  rename(age = mother_age, country_full = country) %>% 
  mutate(measure = "GAMMA1.hat")

# Plot raw curves

lookup_comp <- c("Component 1: Young offspring deaths", "Component 2: Adult offspring deaths")
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
    , label = "Relative\nimportance"
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
    , family = "Times"
    , data = w_lab
    , show.legend = F
  ) +
  # Text to identify components
  geom_label(
    aes(label = label, colour = measure)
    , data = w_xy
    , size = 3
    , family = "Times"
    , show.legend = F
  ) +
  # Real values
  geom_point(
    data = .  %>% filter(measure == "count") 
    , show.legend = F
  ) +
  # To fix axis
  geom_point(data = fix_axis2, colour = "white") + 
  scale_color_viridis_d("", labels = lookup_comp) +
  scale_fill_viridis_d("", labels = lookup_comp) +
  scale_x_continuous("Age of mothers born in 1950") +
  scale_y_continuous(
    "Offspring deaths (counts)", breaks = scales::pretty_breaks(n=5)
    , labels = function(x) format(x, big.mark = ",")
    ) +
  facet_wrap(~country_full, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "bottom"
    , strip.background = element_rect(fill = "honeydew2", colour = "black")
    , text=element_text(family="Times")
  )

# ggsave("figures/summary/!3_counts_components.pdf", width = 6, height = 3, units = "in")
ggsave("figures/summary/Fig5.pdf", width = 6, height = 3, units = "in")

# 4. Summary measures (decomposed) ========

# rcode oduha8ij

lookup_comp2 <- c("Component 1: Young offspring deaths", "Component 2: Adult offspring deaths")
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
    , l1 = gsub("SD", "Standard deviation", l1)
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

# plot without weights

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
    , family = "Times"
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
  scale_y_continuous("Maternal age at offspring loss (years)") +
  theme_bw() +
  theme(
    legend.position = "bottom"
    , plot.margin = margin(0,0, 0, 0, "cm")
    , strip.background = element_rect(fill = "honeydew2", colour = "black")
    , text=element_text(family="Times")
     ) 

# ggsave("figures/summary/!4_summary_decomp1.pdf", width = 6, height = 5, units = "in")
ggsave("figures/summary/Fig7.pdf", width = 6, height = 5, units = "in")

# plot of weights ====

# DF of weights

w_df <- 
  all_measures %>% 
  filter(measure == "w1") %>%
  filter(cohort >= year_low) %>% 
  mutate(
    cohort = as.numeric(cohort)
    , l2 = "Component 1"
  ) %>% 
  select(country, country_full, cohort, mort, value, l2)

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

w_df %>%
  filter(l2  =='Component 1: Young offspring deaths') %>% 
  ggplot() +
  # backgprund lines
  geom_line(
    aes(x = cohort, y = value, colour = mort, group = country)
    , data = . %>% filter(!country %in% country_keep) 
    , size = 0.3
    , alpha = 0.4
    # , show.legend = F
  ) +  
  # HIghlighted lines
  geom_line(
    aes(x = cohort, y = value, colour = mort)
    , data =  . %>%  
      filter(country %in% country_keep) 
    
    , size = 1
    , show.legend = F
  ) +
  # Country labels
  geom_label_repel(
    aes(x = cohort, y = value, colour = mort, label = lookup_c[country])
    , data = . %>% filter(
      country %in% country_keep
      , cohort == 1950
    )
    , nudge_x = -2
    , size = 3
    , family = "Times"
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
    # , show.legend = F
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
  scale_color_brewer("Mortality group", palette = pal) +
  scale_shape("Mortality group") +
  labs(x = "Mother's birth cohort", y = "Relative importance\n(young component)") +
  coord_cartesian(ylim = c(0,1))+
  theme_bw() +
  theme(
    legend.position = "bottom"
    , text=element_text(family="Times")
  ) 

# ggsave("figures/summary/!5_rel_imp.pdf", width = 4.5, height = 3, units = "in")
ggsave("figures/summary/Fig6.pdf", width = 4.5, height = 3, units = "in")

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
    , family = "Times"
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
  scale_y_continuous("Percentage of offspring expected\nto die before their mother", labels = function(x) paste0(round(x*100), "")) +
  theme_bw() +
  theme(
    legend.position = "bottom"
    , strip.background = element_rect(fill = "honeydew2", colour = "black")
    , text=element_text(family="Times")
  ) 

# ggsave("figures/summary/!5_share_child_die.pdf", width = 4, height = 3, units = "in")
ggsave("figures/summary/Fig2.pdf", width = 4, height = 3, units = "in")

# 6. Create Lexis diagram ----

r <- 0:50
mu <- 30
shift <- 5
j <- mu + shift
max_j <- max(r) - shift - mu
max_i <- max(r) - shift
line_size <- 1

br <- c(shift, j, max(r))
labs <- c("t - i", "t - j", "t")
labs <- paste0("Year\n", labs)

br_y <- c(max_j, max_i)
labs_y <- c("j", "i")
labs_y <- paste0("Age\n   ", labs_y)

lex <- 
  data.frame(x = r) %>% 
  mutate(
    mother = x - shift
    , child = mother - mu
  ) %>% 
  pivot_longer(-x) 

labs1 <- data.frame(x = 48, value = 17, label = "Child deaths aged j\nto mothers aged i")

labs2 <- data.frame(x = 12, value = 10, label = "Mother's birth cohort")

labs3 <- data.frame(
  x = max(r) - (max(r) - mu) + 4
  , value = (max_i - max_j)
  , label = "Childbirth"
)

ggplot(mapping = aes(x = x, y = value)) +
  # lexis lines
  geom_line(
    aes(group = name)
    , data = lex
    , size = line_size
  ) +
  # vertical line
  geom_line(
    data = data.frame(x = shift + mu, value = c(0, mu))
    , linetype = "dashed"
    , size = line_size
  ) +
  # show child death with a cross
  geom_point(
    data = data.frame(x = max(r), value = max_j)
    , shape = 3
    , size = 2
    , stroke = 2
    # , colour = "red"
  )+
  # Child deaths + childbirth
  geom_label_repel(
    aes(label = label)
    , nudge_x = -14
    , nudge_y = 1
    , label.size = NA
    , arrow = arrow(length = unit(0.015, "npc"))
    , data = labs1 %>% bind_rows(labs3)
    , family = "Times"
  ) +
  # Mother cohort
  geom_text(
    aes(label = label)
    , angle = 45
    , data = labs2
    , family = "Times"
  ) +
  scale_x_continuous(breaks = br, labels = labs) +
  scale_y_continuous(breaks = br_y, labels = labs_y, position = "right") +
  coord_equal(xlim = c(min(r), max(r)+3), ylim = range(r), expand = F) +
  theme_bw() +
  theme(
    axis.title.x=element_blank()
    , axis.text.x = element_text(size = 14)
    , axis.title.y=element_blank()
    , axis.text.y = element_text(size = 14)
    # , axis.text.y=element_blank()
    # , axis.ticks.y=element_blank()
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , plot.margin=grid::unit(c(0,0,0,0), "mm")
    , text=element_text(family="Times")
  )

ggsave("figures/Fig1.pdf", width = 4.5, height = 4.5, units = "in")  
