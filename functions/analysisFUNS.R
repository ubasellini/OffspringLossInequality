## --------------------------------------------------------- ##
##
## Routines to reproduce the results of the paper:
## "When do parents bury a child? Quantifying uncertainty
## in the parental age at offspring loss"
## 
## Functions to reformat data and conduct analysis of age
## at offspring loss. 
##
## Code by Diego Alburez-Gutierrez and Ugofilippo Basellini (2022) 
## unless otherwise stated
##
## --------------------------------------------------------- ##

# 1. For estimating matrices ----------

# Takes df of period rates in
# 1 age and 1 year groups from UNWPP and returns 
# cohort estimates including all possible combination
approximate_cohort_from_period <- function(per_df, add_all_combinations = T){
  
  coh_df <- 
    per_df %>% 
    group_by(LocID) %>% 
    mutate(cohort = Year - age) %>% 
    arrange(LocID, cohort, age) %>% 
    ungroup
  
  # Add all cohort/age combinations
  if(add_all_combinations){
    coh <- sort( unique(coh_df$cohort) )
    a <- sort(unique(coh_df$age))
    con <- as.character( unique(coh_df$LocID) )
    
    len_con <- length(a) * length(coh)
    len_a <- length(con) * length(coh)
    
    coh2 <- sort(rep(coh, length(a)))
    
    comb <- data.frame(
      LocID = sort(rep(con, len_con))
      , cohort = rep(coh2, length(con))
      , age = rep(a, len_a)
    ) %>% 
      arrange(LocID, cohort, age) 
    
    # Merge with created data structure
    out <- 
      merge(
        coh_df
        , comb
        , all = T
      ) %>% 
      arrange(LocID, cohort, age) %>% 
      data.frame()
  } 
  
  return(coh_df)
  
}

# This fuction takes one country period LT data and returns
# the approximated cohort values
convert_period_lt_to_cohort_lt <- 
  function(lt_per, years = 1765:2100, ages = 1:100) {
    
    # 1. Create cohort life table
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # LTC ... Cohort life table
    LTP_df <- 
      lt_per %>% 
      select(Country = LocID, Year, Age = age, mx, qx, ax) %>% 
      mutate(
        Cohort = Year - Age
        , Country = as.character(Country)
      ) %>% 
      arrange(Country, Cohort, Age) %>% 
      filter(Cohort > min(years) - 1) %>% 
      # Assign to new cohort life tables:
      # lx (number of people left alive)
      # dx (annual number of deaths)
      # nLx (person-live yeares lived in interval)
      mutate(
        lx = ifelse(Age == 0, 1, NA)
        , dx = NA
        , nLx = NA
        , Tx = NA
        , ex = NA
      ) %>% 
      select(Country, Cohort, Year, Age, mx, qx, ax, lx, dx, nLx, Tx, ex)
    
    # Sequentially, the function takes about 1 minute per country (for all years)
    # to exeute currently, so estimate around 2-3 hours for all countries
    # Parallelised on 4 cores is much faster, taking about 8-10 minutes
    
    # Return a df with new columns added
    
    # Makes sure that maximum cohort is the one for which 
    # all age groups are available in the projected data
    year_limit_max <- max(years) - max(ages)
    age_limit_max <- max(ages)
    
    # 1. lx, dx, nLx columns (ages 0-99)
    
    for (year in years){
      # print(year)
      for (age in ages){
        # print(age)
        # Get lagged and normal indices
        current_position <- which(LTP_df$Cohort == year & LTP_df$Age == age)
        lag_position <- which(LTP_df$Cohort == year & LTP_df$Age == (age - 1))
        
        # Assign new lx
        # Since the vector 'ages' starts at 1, the lagged qx value will be that of age 0 always
        LTP_df$lx[current_position] <- LTP_df$lx[lag_position] * (1 - LTP_df$qx[lag_position])
        
        # Assign new dx
        LTP_df$dx[lag_position] <- LTP_df$lx[lag_position] * LTP_df$qx[lag_position]
        
        # Assign nLx column
        if(age < age_limit_max + 1 & year <= year_limit_max) {
          LTP_df$nLx[lag_position] <- 
            LTP_df$lx[current_position] + LTP_df$dx[lag_position] * LTP_df$ax[lag_position]
        }
      }
    }
    
    # 2. dx and nLx columns for open age interval
    
    open_age <-  which(LTP_df$Age == age_limit_max & LTP_df$Cohort <= year_limit_max)
    last_age <- which(LTP_df$Age== age_limit_max - 1  & LTP_df$Cohort <= year_limit_max)
    open_age2 <- which(LTP_df$Age == age_limit_max)
    
    LTP_df$dx[open_age] <- LTP_df$lx[last_age] * LTP_df$qx[last_age]
    LTP_df$nLx[open_age2] <- LTP_df$ax[open_age2] * LTP_df$dx[open_age2]
    
    # 3. Tx column (persons alive at any point in time above age x)
    
    for (year in years){
      if (year <= year_limit_max)
        LTP_df$Tx[LTP_df$Cohort == year] <- rev(cumsum(rev(LTP_df$nLx[LTP_df$Cohort == year])))
    }
    
    # 4. Compute ex
    
    LTP_df$ex<-LTP_df$Tx/LTP_df$lx
    
    # 5. Return only relevant rows
    
    LTP_df$Year <- NULL
    
    LTC_df <- LTP_df
    
    return(LTC_df)
  }


# Read different data sources containing population and period rate data and 
# transform them to return, where appropriate, cohort rates for single-age
# and single-years.
get_denmark_hmd_hfd_data_clean2 <- function(max.age, years_unwpp, cohorts_unwpp, expand_lt_to_future_stable = T, export = F){
  
  country_name <- "Denmark"
  lookup_codes <- c("SWE", "DNK")
  names(lookup_codes) <- c("Sweden", "Denmark")
  
  code <- lookup_codes[country_name]
  
  # Check that all the needed data are available
  datasets_needed <- c(
    # HMD period life tables
    "ltper_1x1" = paste0("data/input/", code, "/fltper_1x1.txt")
    # HMD period exposures
    , "Exposures_1x1" = paste0("data/input/", code, "/Exposures_1x1.txt")
    # HFD Period age-specific fertility rates
    , "asfrRR" = paste0("data/input/", code, "/", code, "asfrRR.txt")
    # HFC Period age-specific fertility rates
    , "ASFRstand_TOT" = paste0("data/input/", code, "/", code, "_ASFRstand_TOT.txt")
    # HMD Yearly population figures
    , "Population" =  paste0("data/input/", code, "/Population.txt")
    # UNWPP period 1x1 life tables (created by diego from abridged data)
    , "lt_per_1_1_F" = paste0("data/input/", code, "/lt_per_1_1_F_",tolower(code),".csv")
    # UNWPP period 1x1 asfr (created by diego from abridged data)
    , "fert_per_1_1" = paste0("data/input/", code, "/fert_per_1_1_",tolower(code),".csv")
    # UNWPP population by sinle-age and sex
    , "WPP2019_Population" = paste0("data/input/", code, "/WPP2019_PopulationBySingleAgeSex_2020-2100_",tolower(code),".csv")
  )
  
  miss <- !file.exists(datasets_needed)
  if(any(miss)) stop(paste("Missing data:", datasets_needed[miss]))
  else print(paste("Reading dataset: ", datasets_needed))
  
  # Get the datasets that do not change
  
  # Period life tables
  # 1. Reformat HMD data to lower upper age group to 100
  
  # Get data period life tables
  rate <- 
    read.table(datasets_needed["ltper_1x1"], header = T, skip = 2) %>% 
    select(Year, age = Age, value = mx) %>% 
    filter(Year %in% years_unwpp) %>% 
    pivot_wider(names_from = "Year", values_from = "value") %>% 
    select(-age) %>% 
    as.matrix()
  
  pop <- 
    read.table(datasets_needed["Exposures_1x1"], header = T, skip = 2) %>% 
    select(Year, age = Age, value = Total) %>% 
    filter(Year %in% years_unwpp) %>% 
    pivot_wider(names_from = "Year", values_from = "value") %>% 
    select(-age) %>% 
    as.matrix()
  
  rownames(rate) <- rownames(pop) <- 0:(dim(pop)[1]-1)
  
  data <- list(
    age = as.numeric(rownames(rate))
    , year = as.numeric(colnames(rate))
    , pop = pop
    , rate = rate
  ) 
  
  # Lower ope age group to 100+
  new_data <- lower_lt_open_age_int(data = data, max.age = max.age)
  
  # To long format
  nmx_df <- 
    new_data$rate %>% 
    data.frame() %>% 
    rownames_to_column(var = "age") %>% 
    pivot_longer(-age, names_to = "year") %>% 
    mutate(
      year = as.numeric(gsub("X", "", year))
      , age = ifelse(age == "100+", 100, age)
      , age = as.numeric(age)
    ) %>% 
    arrange(year, age) 
  
  # Build period life tables from nmx column
  # This uses Ugo's lifeable.mx function
  lt_1_1 <- 
    nmx_df %>% 
    group_by(year) %>% 
    group_map(~ lifetable.mx(x = .$age, mx = .$value, sex= "F", ax = NULL, radix = 1e5)) %>% 
    bind_rows() %>% 
    mutate(year = nmx_df$year) %>% 
    select(Year = year, age = x, nax = ax, lx, nmx = mx, nqx = qx, ndx = dx, nLx = Lx, Tx, ex) %>% 
    # Since this is already given in UNWPP
    filter(Year < 2020)
  
  # print("Getting ASFR from https://www.humanfertility.org/cgi-bin/getfile.plx?f=SWE\20201118\SWEasfrRR.txt")
  
  # Period age-specific fertility rates
  asfr_1_1 <-
    read.table(datasets_needed["asfrRR"], header = T, skip = 2) %>% 
    mutate(LocID = 1) %>% 
    select(LocID, Year, age = Age, asfr = ASFR) %>% 
    filter(age %in% 15:50) %>% 
    filter(Year %in% years_unwpp) %>% 
    filter(Year < 2020) %>% 
    mutate(age = as.numeric(age))    
  
  # Add HFC data for historical data
  print("Getting ASFR from Human Fertility Collection")
  
  # Codebook:
  # https://www.fertilitydata.org/docs/formats.pdf
  # Period age-specific fertility rates
  asfr_hfc_5_1 <-
    read.table(datasets_needed["ASFRstand_TOT"], header = T, sep = ",", stringsAsFactors = F) %>% 
    # 1916 is first year for which HFD data is avilable
    filter(Year1 < 1916) %>%
    mutate(Collection = trimws(Collection)) %>% 
    arrange(Year1, Age, Collection) %>% 
    filter(Collection == "RE") %>% 
    mutate(LocID = 1) %>% 
    select(LocID, Year = Year1, age = Age, asfr = ASFR) %>% 
    filter(age %in% 15:50) %>%
    mutate(
      asfr = trimws(asfr) 
      , asfr = as.numeric(ifelse(asfr == ".", 0, asfr))
      , age = as.numeric(age)
    )    
  
  # Use linear interpolation between calendar years to get single-years
  
  y_int <- min(asfr_hfc_5_1$Year):1915
  
  asfr_hfc_1_1 <-
    asfr_hfc_5_1 %>% 
    group_by(LocID, age) %>%
    group_map(
      ~ interpolate_COLUMN_calendar_years(., column = "asfr", years = y_int)
      , .keep = T
    ) %>% 
    bind_rows() %>% 
    ungroup() %>% 
    filter(!is.na(Year)) %>% 
    arrange(LocID, Year, age)
  
  asfr_1_1 <- bind_rows(asfr_hfc_1_1, asfr_1_1)
  
  # HMD Yearly population figures
  pop_1_1_temp <- 
    read.table(datasets_needed["Population"], header = T, skip = 2) %>% 
    mutate(LocID = 1) %>% 
    select(LocID, Year, age = Age, PopMale = Male, PopFemale = Female, PopTotal = Total) %>% 
    mutate(
      age = ifelse(age == "110+", 110, age)
      , age = as.numeric(age)
    ) %>% 
    filter(Year < 2020)
  
  # There are two values for Denmark 1921, keep first
  
  pop_1_1_temp <- 
    pop_1_1_temp %>% 
    filter(Year != "1921-") %>% 
    mutate(Year = ifelse(Year == "1921+", "1921", Year))
  
  # Sum all pop over 100 to 100
  over100 <- 
    pop_1_1_temp %>% 
    filter(age >= 100) %>% 
    group_by(LocID, Year) %>% 
    summarise(
      PopMale = sum(PopMale)
      , PopFemale = sum(PopFemale)
      , PopTotal = sum(PopTotal)
      , age = 100
    ) %>% 
    ungroup()
  
  pop_1_1 <- 
    pop_1_1_temp %>% 
    filter(age < 100) %>% 
    bind_rows(over100) %>% 
    mutate(Year = as.numeric(Year)) %>% 
    arrange(LocID, Year, age)
  
  print("Using mortality projetions from UNWPP after 2019.")
  
  # Now, convert period to cohort
  # Steps: 
  # 1. Get 1x1 UNWPP and HMD life tables
  # 2. Merge with historical Sweden data (make sure that upper age bounds fit)
  # 3. Approximate cohort from period mortality rates
  
  # 1. Get data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  lt_un <- 
    read.csv(datasets_needed["lt_per_1_1_F"], stringsAsFactors = F) %>% 
    data.frame() %>% 
    filter(country == tolower(country_name)) %>% 
    mutate(LocID = 1) %>% 
    select(LocID, Year = year, age, ax, lx, mx, qx, dx, Lx, Tx, ex) %>%
    # Keep only projections
    filter(Year >= 2020)
  
  
  # HMD Period life tables
  lt_hmd <- 
    lt_1_1 %>% 
    mutate(LocID = 1) %>% 
    select(LocID, Year, age, ax = nax, lx, mx = nmx, qx = nqx, dx = ndx, Lx = nLx, Tx, ex) %>% 
    # rescale so that radix matches UN radix
    mutate(lx = lx * 10, dx = dx* 10, Lx = Lx * 10, Tx = Tx * 10)
  
  # 2. merge historical and projected LTs
  lt_1_1_temp <- 
    bind_rows(lt_un, lt_hmd)  %>% 
    arrange(LocID, Year, age)
  
  if(expand_lt_to_future_stable){
    
    # Expand cohort life tables to the future assuming
    # demographic stability. 
    # This is needed because the function that estimates rates
    # needs to have cohort lifetable for cohorts c:(c+50), where
    # 'c' is the mother's birth cohort. The approximation of cohort 
    # LTs above only produces complete LT for cohorts 1950:2100 (i.e. with 
    # values gor ages 0:100). We have data for older cohorts but it is incomplete 
    # for all ages. For example, LT for cohort 1951 only goes up to age 99 and so on. 
    # One solution is to assume that the values for age 100 in the 1951 cohort LT
    # is equivalent to the last observed values, i.e., the values for age 100 in 
    # the 1950 cohort LT, and so on. 
    
    print(paste0("Expanding period LT beyond 2000 by assuming that values remains constant after period 2100 at the 2100 levels"))
    
    lt_per <- 
      lt_1_1_temp %>% 
      bind_rows(
        lt_1_1_temp %>% 
          group_by(LocID) %>% 
          filter(Year == max(.$Year)) %>% 
          group_map(~repeat_lt(., up_to_year = up_to_year), .keep = T) %>% 
          bind_rows() %>% 
          ungroup()     
      ) %>% 
      arrange(LocID, Year, age)
    
  } else {
    lt_per <- 
      lt_1_1_temp %>% 
      arrange(LocID, Year, age)
  }
  
  # 3. Approximate cohort from period mortality rates
  print("Converting period to cohort LT...")
  
  lt_coh <- 
    convert_period_lt_to_cohort_lt(lt_per, years = years_unwpp, ages = 1:max.age)  %>% 
    select(LocID = Country, age = Age, nax = ax, lx, nmx = mx, nqx = qx, ndx = dx, nLx, Tx, ex, cohort = Cohort)
  
  # B. Combine period LTs
  print("Combining period LT from HMD and UNWPP" )
  
  lt_1_1 <- 
    lt_1_1_temp %>% 
    select(LocID, Year, age, nax = ax, lx, nmx = mx, nqx = qx, ndx = dx, nLx = Lx, Tx, ex)
  
  # B. Merge HFD and UNWPP period fertility rates
  print("Combining period fertility from HFD and projections > 2020 from UNWPP")
  
  # These are UNWPP period fert projections
  asfr_wpp_proj <- 
    read.csv(datasets_needed["fert_per_1_1"], stringsAsFactors = F) %>% 
    # 752 is Sweden
    filter(country == tolower(country_name)) %>% 
    filter(year >= 2020) %>% 
    data.frame() %>% 
    mutate(
      LocID = 1
      , value = value /1000
    ) %>% 
    select(LocID, Year = year, age, asfr = value)
  
  asfr_1_1 <- 
    bind_rows(asfr_1_1, asfr_wpp_proj) %>% 
    arrange(LocID, Year, age)
  
  # C. Combine population
  print("Combining population totals from HMD and UNWPP" )
  
  # UNWPP population by sinle-age and sex
  pop_wpp <- 
    read.csv(datasets_needed["WPP2019_Population"], stringsAsFactors = F) %>% 
    data.frame() %>% 
    filter(Location == country_name, Variant == "Medium") %>% 
    mutate(
      LocID = 1
      , PopFemale = PopFemale * 1000
      , PopMale  = PopMale * 1000
      , PopTotal = PopTotal * 1000
    ) %>% 
    select(LocID, Year = Time, age = AgeGrp, PopMale, PopFemale, PopTotal) %>% 
    filter(Year >= 2020)
  
  pop_1_1 <- 
    bind_rows(pop_1_1, pop_wpp)  %>% 
    arrange(LocID, Year, age)
  
  # Filter before export
  lt_coh <- lt_coh  %>% 
    filter(cohort %in% cohorts_unwpp)
  
  lt_1_1 <- lt_1_1 %>% 
    filter(Year %in% years_unwpp)
  
  asfr_1_1 <- asfr_1_1 %>%
    filter(Year %in% years_unwpp)
  
  pop_1_1 <- pop_1_1 %>% 
    filter(Year %in% years_unwpp)
  
  if(export){
    print("Saving ASFR datasets to data/input/asfr_1_1_DNK.csv")
    write.csv(asfr_1_1, "data/input/asfr_1_1_DNK.csv", row.names = F)
  }
  
  out <- list(lt_coh = lt_coh, lt_1_1 = lt_1_1, asfr_1_1 =  asfr_1_1, pop_1_1 = pop_1_1)
  # asfr_coh = asfr_coh, 
  return(out)
}

# Optmised function to get cohort estiamtes from DemoKin
# for multuple countries and selected relatives only
# Note that this differs from 'standard' DemoKin 
# implementarion but should produce the same results
get_kin_trees_unwpp <- function(
  asfr_1_1, lt_1_1, pop_1_1, cohorts
  , age = 100, max_age = 100, years_unwpp, base = "Output/demokin/", pat = "tree"
  , parallel = F, cores = 2
){
  
  if(!parallel){
    # Iterate over all countries 

    for(con in unique(asfr_1_1$LocID)){
      for(cohort in cohorts){
        print(paste0("Working on country: ", con, " - ",cohort))
        print(Sys.time())
        
        out_file <- paste0(base, pat, "_", con,"_", cohort, ".rds")
        
        # Get rates to input in DemoKin function
        if(file.exists(out_file)){
          print(paste0(out_file, " already exists. I'll just start with the next one."))
        } else {
          ages_v <- 0:max_age
          ages <- length(ages_v)
          w <- max_age
          
          L <- 
            lt_1_1 %>% 
            filter(LocID == con) %>% 
            filter(age <= w) %>% 
            mutate(nLx = ifelse(age == w, Tx, nLx)) %>% 
            select(Year, age, nLx) %>% 
            spread(Year, nLx) %>% 
            select(-age)
          
          # The argument P in DemoKin::kins() is S(x)=L(x+1)/L(x)
          # The first row should be L(1)/L(0) 
          
          age_plus1 <- max_age + 1
          
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
            filter(LocID == con) %>% 
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
            filter(LocID == con) %>% 
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
    } # end country loop
  } else if(parallel) {
    
    print(paste0("Estimating kin in parallel using ", cores, " cores. Good luck!"))
    
    cl <- makeCluster(cores) #not to overload your computer
    registerDoParallel(cl)
    
    # Iterate over all countries 
    # con <- 752; cohort <- 1995
    # for(con in unique(asfr_1_1$LocID)){
    foreach(
      con = unique(asfr_1_1$LocID)
      , .packages = c("dplyr", "tidyr")
      , .export = c("kins_fast", "kins_non_stable_fast")
    ) %dopar% {
      for(cohort in cohorts){
        print(paste0("Working on country: ", con, " - ",cohort))
        print(Sys.time())
        
        out_file <- paste0(base, pat, "_", con,"_", cohort, ".rds")
        
        # Get rates
        if(file.exists(out_file)){
          print(paste0(out_file, " already exists. I'll just start with the next one."))
        } else {
          # The argument P in the matrix models is S(x)=L(x+1)/L(x)
          # The first row should be L(1)/L(0) 
          ages_v <- 0:max_age
          ages <- length(ages_v)
          w <- max(ages_v)
          
          L <- 
            lt_1_1 %>% 
            filter(LocID == con) %>% 
            filter(age <= w) %>% 
            mutate(nLx = ifelse(age == w, Tx, nLx)) %>% 
            select(Year, age, nLx) %>% 
            spread(Year, nLx) %>% 
            select(-age)
          
          age_plus1 <- max_age + 1
          
          P <- 
            rbind(
              L[c(-1, -age_plus1), ]/L[-c(max_age:age_plus1), ]
              , L[age_plus1, ]/(L[max_age, ] + L[age_plus1, ])
              , L[age_plus1, ]/(L[max_age, ] + L[age_plus1, ])
            )
          
          rownames(P) <- ages_v
          
          asfr_mat <- 
            asfr_1_1 %>% 
            filter(LocID == con) %>% 
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
            filter(LocID == con) %>% 
            select(Year, age, value = PopFemale) %>% 
            pivot_wider(names_from = Year, values_from = value) %>% 
            select(-age) %>% 
            data.frame()
          
          colnames(pop_mat) <- years_unwpp
          
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
    } # end country loop
    
    stopCluster(cl)
    
  }
  
  print("Kinship models done!")
  
}

# gET LX COLUMNS OF A LIFE TABLE GIVEN A RADIX
# AND VECTOR OF MORT PROBS
get_lx <- function(l0, nqx){ 
  lx <- cumprod(c(l0, 1 - nqx)) 
  lx[-length(lx)]
}


# Same but returns df
# Use with 1x1 df
get_lt_bulk_per <- function(data, age, add_country_year = F){
  
  lx <- data$lx
  
  n <- 1
  nage <- length(lx)
  nax <- 0.5 * n
  ndx <- -diff(lx)
  ndx <- c(ndx, lx[nage])
  
  lxpn <- lx[-1]
  nLx <- n * lxpn + ndx[-nage] * nax
  nLx[nage] <- lx[nage]/0.5  
  
  nmx <- ndx/nLx
  nqx <- 1 - lead(lx)/lx
  nqx[nage] <- 1
  
  Tx <- rev(cumsum(rev(nLx)))
  ex <- Tx/lx
  
  out <- data.frame(
    age = age
    , nax
    , lx
    , nmx
    , nqx
    , ndx
    , nLx
    , Tx
    , ex
  )  
  
  if(add_country_year){
    out$LocID <- unique(data$LocID)  
    out$Year <- unique(data$Year)
  } else {
    out$id <- sample(1:1e6, 1)
  }
  
  return(out)
  
}

# Given a set of parameters and objects, produce the rdata file
# that Ugo uses for the decomposition
get_matrices_wpp <- function(country_keep, return_list = F){
  
  print(paste("Starting to work on ", country_keep))
  
  base2 <- ifelse(base == "", ".", base)
  full.names <- base != ""
  
  files <- list.files(path = base2, pattern = pat, full.names = full.names)
  
  # A df to identify the DemoKin output files
  
  files_df <- 
    data.frame(files = files) %>% 
    mutate(
      f = gsub(paste0(base, pat, "_"), "", files)
      , f = gsub(".rds", "", f)
      , LocID = str_extract(f, "^[0-9]+")
      , country = lookup_c[LocID]
      , cohort = str_extract(f, "[0-9]{4}$")
    ) %>% 
    select(-f) %>% 
    filter(!is.na(country))
  
  files_n <- 
    files_df %>% 
    filter(LocID %in% lookup_c2[country_keep]) %>% 
    pull(files)
  
  print("Reformating DemoKin outputs into arrays...")
  
  # Get a list of arrays, one for each relatve type
  # The idea is to have the rows represnt ego's age
  # columns reresent ego's age, and the third dimension
  # be ego's birth cohort
  arys <- reformat_kin_trees(files_n, dataset = "kins", relatives = c("d", "m"))
  
  # 'd' stands for daughters in DemoKin
  H <- arys[["d"]]
  
  # Let's assign the cohorts as names to this array
  dimnames(H)[[3]] <-
    (1:dim(H)[3]) + min(cohorts_unwpp) - 1
  
  lt  <-
    lt_coh %>%
    # Keep country of interest
    filter(LocID %in% lookup_c2[country_keep]) %>%
    group_by(LocID, cohort) %>%
    mutate(
      # lx should be normalised to be a density function
      lx = lx/first(lx)
    ) %>%
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
  # 20220210
  # NEW APPROACH TO GET MATRIX OF pOPULATION
  # Estimate from cohort life tables given a birth cohort size
  
  l0_df <- 
    pop_1_1 %>% 
    filter(LocID %in% lookup_c2[country_keep]) %>%
    filter(age == 0) %>% 
    filter(Year %in% cohorts_demokin) %>% 
    select(year = Year, age, value = PopFemale)
  
  l0 <- l0_df$value
  names(l0) <- l0_df$year
  
  POP <- mapply(
    get_lx
    , l0
    , as.data.frame(Q)
    , SIMPLIFY = T
  )
  
  # dim(POP)
  
  # Export arrays 
  out_name <- countrycode::countrycode(country_keep, "country.name", "iso3c")
  out_file <- paste0("data/clean/", out_name,"_arrays_clean.rdata")
  print(paste0("Saving file to disk: ", out_file))
  
  # rcode lksdn437
  save(H, D, L, Q, POP, file = out_file)
  
  print("Done! existing function...")
  
  if(return_list){
    return(list(H, D, L, Q, POP))
  }
  
}

# Read different data sources containing population and period rate data and 
# transform them to return, where appropriate, cohort rates for single-age
# and single-years.
get_sweden_hmd_hfd_data_clean2 <- function(max.age, years_unwpp, cohorts_unwpp, expand_lt_to_future_stable = T, export = F){
  
  # Check that all the needed data are available
  datasets_needed <- c(
    # HMD period life tables
    "data/input/SWE/fltper_1x1.txt"
    # HMD period exposures
    , "data/input/SWE/Exposures_1x1.txt"
    # HFD Period age-specific fertility rates
    , "data/input/SWE/SWEasfrRR.txt"
    # HDC Period age-specific fertility rates
    , "data/input/SWE/SWE_ASFRstand_TOT.txt"
    # HMD Yearly population figures
    , "data/input/SWE/Population.txt"
    # UNWPP period 1x1 life tables (created by diego from abridged data)
    , "data/input/SWE/lt_per_1_1_F_swe.csv"
    # UNWPP period 1x1 asfr (created by diego from abridged data)
    , "data/input/SWE/fert_per_1_1_swe.csv"
    # UNWPP population by sinle-age and sex
    , "data/input/SWE/WPP2019_PopulationBySingleAgeSex_2020-2100_swe.csv"
  )
  
  miss <- !file.exists(datasets_needed)
  if(any(miss)) stop(paste("Missing data:", datasets_needed[miss]))
  else print(paste("Reading dataset: ", datasets_needed))
  
  # Get the datasets that do not change
  
  # Period life tables
  # 1. Reformat HMD data to lower upper age group to 100
  
  # Get data period life tables
  rate <- 
    read.table("data/input/SWE/fltper_1x1.txt", header = T, skip = 2) %>% 
    select(Year, age = Age, value = mx) %>% 
    filter(Year %in% years_unwpp) %>% 
    pivot_wider(names_from = "Year", values_from = "value") %>% 
    select(-age) %>% 
    as.matrix()
  
  pop <- 
    read.table("data/input/SWE/Exposures_1x1.txt", header = T, skip = 2) %>% 
    select(Year, age = Age, value = Total) %>% 
    filter(Year %in% years_unwpp) %>% 
    pivot_wider(names_from = "Year", values_from = "value") %>% 
    select(-age) %>% 
    as.matrix()
  
  rownames(rate) <- rownames(pop) <- 0:(dim(pop)[1]-1)
  
  data <- list(
    age = as.numeric(rownames(rate))
    , year = as.numeric(colnames(rate))
    , pop = pop
    , rate = rate
  ) 
  
  # Lower ope age group to 100+
  new_data <- lower_lt_open_age_int(data = data, max.age = max.age)
  
  # To long format
  nmx_df <- 
    new_data$rate %>% 
    data.frame() %>% 
    rownames_to_column(var = "age") %>% 
    pivot_longer(-age, names_to = "year") %>% 
    mutate(
      year = as.numeric(gsub("X", "", year))
      , age = ifelse(age == "100+", 100, age)
      , age = as.numeric(age)
    ) %>% 
    arrange(year, age) 
  
  # Build period life tables from nmx column
  # This uses Ugo's lifeable.mx function
  lt_1_1 <- 
    nmx_df %>% 
    group_by(year) %>% 
    group_map(~ lifetable.mx(x = .$age, mx = .$value, sex= "F", ax = NULL, radix = 1e5)) %>% 
    bind_rows() %>% 
    mutate(year = nmx_df$year) %>% 
    select(Year = year, age = x, nax = ax, lx, nmx = mx, nqx = qx, ndx = dx, nLx = Lx, Tx, ex) %>% 
    # Since this is already given in UNWPP
    filter(Year < 2020)
  
  # print("Getting ASFR from https://www.humanfertility.org/cgi-bin/getfile.plx?f=SWE\20201118\SWEasfrRR.txt")
  
  # Period age-specific fertility rates
  asfr_1_1 <-
    read.table("data/input/SWE/SWEasfrRR.txt", header = T, skip = 2) %>% 
    mutate(LocID = 1) %>% 
    select(LocID, Year, age = Age, asfr = ASFR) %>% 
    filter(age %in% 15:50) %>% 
    filter(Year %in% years_unwpp) %>% 
    mutate(age = as.numeric(age))    
  
  # Add HFC data for historical data
  print("Getting ASFR from Human Fertility Collection")
  
  # Codebook:
  # https://www.fertilitydata.org/docs/formats.pdf
  # Period age-specific fertility rates
  asfr_hfc_5_1 <-
    read.table("data/input/SWE/SWE_ASFRstand_TOT.txt", header = T, sep = ",", stringsAsFactors = F) %>% 
    filter(Year1 < 1891) %>% 
    mutate(Collection = trimws(Collection)) %>% 
    arrange(Year1, Age, Collection) %>% 
    filter(Collection == "STAT") %>% 
    mutate(LocID = 1) %>% 
    select(LocID, Year = Year1, age = Age, asfr = ASFR) %>% 
    filter(age %in% 15:50) %>%
    mutate(
      asfr = trimws(asfr) 
      , asfr = as.numeric(ifelse(asfr == ".", 0, asfr))
      , age = as.numeric(age)
    )    
  
  # Use linear interpolation between calendar years to get single-years
  
  y_int <- min(asfr_hfc_5_1$Year):1890
  
  asfr_hfc_1_1 <-
    asfr_hfc_5_1 %>% 
    group_by(LocID, age) %>%
    group_map(
      ~ interpolate_COLUMN_calendar_years(., column = "asfr", years = y_int)
      , .keep = T
    ) %>% 
    bind_rows() %>% 
    ungroup() %>% 
    filter(!is.na(Year)) %>% 
    arrange(LocID, Year, age)
  
  asfr_1_1 <- bind_rows(asfr_hfc_1_1, asfr_1_1)
  
  # HMD Yearly population figures
  pop_1_1_temp <- 
    read.table("data/input/SWE/Population.txt", header = T, skip = 2) %>% 
    mutate(LocID = 1) %>% 
    select(LocID, Year, age = Age, PopMale = Male, PopFemale = Female, PopTotal = Total) %>% 
    mutate(
      age = ifelse(age == "110+", 110, age)
      , age = as.numeric(age)
    ) %>% 
    filter(Year %in% years_unwpp)
  
  # Sum all pop over 100 to 100
  over100 <- 
    pop_1_1_temp %>% 
    filter(age >= 100) %>% 
    group_by(LocID, Year) %>% 
    summarise(
      PopMale = sum(PopMale)
      , PopFemale = sum(PopFemale)
      , PopTotal = sum(PopTotal)
      , age = 100
    ) %>% 
    ungroup()
  
  pop_1_1 <- 
    pop_1_1_temp %>% 
    filter(age < 100) %>% 
    bind_rows(over100) %>% 
    arrange(LocID, Year, age)
  
  
  print("Using mortality projetions from UNWPP after 2019.")
  
  # Now, convert period to cohort
  # Steps: 
  # 1. Get 1x1 UNWPP and HMD life tables
  # 2. Merge with historical Sweden data (make sure that upper age bounds fit)
  # 3. Approximate cohort from period mortality rates
  
  # 1. Get data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  lt_un <- 
    read.csv("data/input/SWE/lt_per_1_1_F_swe.csv", stringsAsFactors = F) %>% 
    data.frame() %>% 
    filter(country == "sweden") %>% 
    mutate(LocID = 1) %>% 
    # select(LocID, Year = year, age, nax = ax, lx, nmx = mx, nqx = qx, ndx = dx, nLx = Lx, Tx, ex) %>% 
    select(LocID, Year = year, age, ax, lx, mx, qx, dx, Lx, Tx, ex) %>%
    # Kepp only projections
    filter(Year >= 2020)
  
  
  # HMD Period life tables
  lt_hmd <- 
    lt_1_1 %>% 
    mutate(LocID = 1) %>% 
    select(LocID, Year, age, ax = nax, lx, mx = nmx, qx = nqx, dx = ndx, Lx = nLx, Tx, ex) %>% 
    # rescale so that radix matches UN radix
    mutate(lx = lx * 10, dx = dx* 10, Lx = Lx * 10, Tx = Tx * 10)
  
  # 2. merge historical and projected LTs
  lt_1_1_temp <- 
    bind_rows(lt_un, lt_hmd)  %>% 
    arrange(LocID, Year, age)
  
  if(expand_lt_to_future_stable){
    
    # Expand cohort life tables to the future assuming
    # demographic stability. 
    # This is needed because the function that estimates rates
    # needs to have cohort lifetable for cohorts c:(c+50), where
    # 'c' is the mother's birth cohort. The approximation of cohort 
    # LTs above only produces complete LT for cohorts 1950:2100 (i.e. with 
    # values gor ages 0:100). We have data for older cohorts but it is incomplete 
    # for all ages. For example, LT for cohort 1951 only goes up to age 99 and so on. 
    # One solution is to assume that the values for age 100 in the 1951 cohort LT
    # is equivalent to the last observed values, i.e., the values for age 100 in 
    # the 1950 cohort LT, and so on. 
    
    print(paste0("Expanding period LT beyond 2000 by assuming that values remains constant after period 2100 at the 2100 levels"))
    
    lt_per <- 
      lt_1_1_temp %>% 
      bind_rows(
        lt_1_1_temp %>% 
          group_by(LocID) %>% 
          filter(Year == max(.$Year)) %>% 
          group_map(~repeat_lt(., up_to_year = up_to_year), .keep = T) %>% 
          bind_rows() %>% 
          ungroup()     
      ) %>% 
      arrange(LocID, Year, age)
    

  } else {
    lt_per <- 
      lt_1_1_temp %>% 
      arrange(LocID, Year, age)
  }
  
  # 3. Approximate cohort from period mortality rates
  print("Converting period to cohort LT...")
  
  lt_coh <- 
    # convert_period_lt_to_cohort_lt(lt_per, years = 1751:2100, ages = 1:100)  %>%
    convert_period_lt_to_cohort_lt(lt_per, years = years_unwpp, ages = 1:max.age)  %>% 
    select(LocID = Country, age = Age, nax = ax, lx, nmx = mx, nqx = qx, ndx = dx, nLx, Tx, ex, cohort = Cohort)
  
  # B. Combine period LTs
  print("Combining period LT from HMD and UNWPP" )
  lt_1_1 <- 
    lt_1_1_temp %>% 
    select(LocID, Year, age, nax = ax, lx, nmx = mx, nqx = qx, ndx = dx, nLx = Lx, Tx, ex)
  
  # B. Merge HFD and UNWPP period fertility rates
  print("Combining period fertility from HFD and projections > 2020 from UNWPP")
  
  # These are UNWPP period fert projections
  asfr_wpp_proj <- 
    read.csv("data/input/SWE/fert_per_1_1_swe.csv", stringsAsFactors = F) %>% 
    # 752 is Sweden
    filter(country == "sweden") %>% 
    filter(year >= 2020) %>% 
    data.frame() %>% 
    mutate(
      LocID = 1
      , value = value /1000
    ) %>% 
    select(LocID, Year = year, age, asfr = value)
  
  asfr_1_1 <- 
    bind_rows(asfr_1_1, asfr_wpp_proj) %>% 
    arrange(LocID, Year, age)
  
  # C. Combine population
  print("Combining population totals from HMD and UNWPP" )
  
  # UNWPP population by sinle-age and sex
  pop_wpp <- 
    read.csv("data/input/SWE/WPP2019_PopulationBySingleAgeSex_2020-2100_swe.csv", stringsAsFactors = F) %>% 
    data.frame() %>% 
    filter(Location == "Sweden", Variant == "Medium") %>% 
    mutate(
      LocID = 1
      , PopFemale = PopFemale * 1000
      , PopMale  = PopMale * 1000
      , PopTotal = PopTotal * 1000
    ) %>% 
    select(LocID, Year = Time, age = AgeGrp, PopMale, PopFemale, PopTotal) %>% 
    filter(Year > 2020)
  
  pop_1_1 <- 
    bind_rows(pop_1_1, pop_wpp)  %>% 
    arrange(LocID, Year, age)
  
  # Filter before export
  lt_coh <- lt_coh  %>% 
    filter(cohort %in% cohorts_unwpp)
  
  lt_1_1 <- lt_1_1 %>% 
    filter(Year %in% years_unwpp)
  
  asfr_1_1 <- asfr_1_1 %>%
    filter(Year %in% years_unwpp)
  
  pop_1_1 <- pop_1_1 %>% 
    filter(Year %in% years_unwpp)
  
  if(export){
    print("Saving ASFR datasets to data/input/asfr_1_1_SWE.csv")
    write.csv(asfr_1_1, "data/input/asfr_1_1_SWE.csv", row.names = F)
  }
  
  out <- list(lt_coh = lt_coh, lt_1_1 = lt_1_1, asfr_1_1 =  asfr_1_1, pop_1_1 = pop_1_1)
  # asfr_coh = asfr_coh, 
  return(out)
}

# Read different data sources containing population and period rate data and 
# transform them to return, where appropriate, cohort rates for single-age
# and single-years.
get_unwpp_cohort_rates <- function(lookup_c, reestimate = F, expand_lt_to_future_stable, up_to_year = NA, cohorts_unwpp, cohort_lt, country_key, export = F){
  
  # if reestimate, downloads raw data from unwpp, ungroups it and does the
  # period approximation
  
  if(reestimate){
    
    # Download from UNWPP website directly if they don't exist
    
    undir <- "data/input/unwpp/"
    fert_name <- "WPP2019_Fertility_by_Age.csv"
    mort_name <- "WPP2019_Life_Table_Medium.csv"
    pop1_name <- "WPP2019_PopulationBySingleAgeSex_1950-2019.csv"
    pop2_name <- "WPP2019_PopulationBySingleAgeSex_2020-2100.csv"
    
    diri_mort <- paste0(undir, mort_name)
    diri_fert <- paste0(undir, fert_name)
    diri_pop1 <- paste0(undir, pop1_name)
    diri_pop2 <- paste0(undir, pop2_name)
    
    # Download takes a while
    if(!file.exists(diri_mort)) {
      print("Downloading data from website...")
      print("Will probably take a while")
      
      url_base <- "https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/"
      
      # FERT
      
      url <- paste0(url_base, fert_name)
      httr::GET(url, write_disk(diri_fert, overwrite = T))
      
      # MORT
      
      url <- paste0(url_base, mort_name)
      httr::GET(url, write_disk(diri_mort, overwrite = T))  
      
      # Pop1
      url <- paste0(url_base, pop1_name)
      httr::GET(url, write_disk(diri_pop1, overwrite = T))  
      
      # Pop2
      url <- paste0(url_base, pop2_name)
      httr::GET(url, write_disk(diri_pop2, overwrite = T))  
      
    }
    
    fert_un <- 
      fread(diri_fert, stringsAsFactors = F) %>% 
      filter(Variant == "Medium") %>% 
      filter(LocID %in% names(lookup_c)) %>% 
      data.frame() %>% 
      select(LocID, Year = MidPeriod, ASFR) 
    # mutate(Year = Year - 3)
    
    lx_un <- 
      fread(diri_mort, stringsAsFactors = F) %>% 
      data.frame() %>% 
      filter(Sex == "Female") %>% 
      filter(LocID %in% names(lookup_c)) %>% 
      select(LocID, Year = MidPeriod, x = AgeGrpStart, lx)
    
    pop_1_1 <- 
      # historical rates
      fread(diri_pop1, stringsAsFactors = F) %>% 
      data.frame() %>% 
      filter(LocID %in% names(lookup_c)) %>% 
      select(LocID, Year = Time, age = AgeGrp, starts_with("Pop")) %>% 
      # projected rates
      bind_rows(
        fread(diri_pop2, stringsAsFactors = F) %>% 
          data.frame() %>% 
          filter(LocID %in% names(lookup_c)) %>% 
          select(LocID, Year = Time, age = AgeGrp, starts_with("Pop"))    
      ) %>% 
      mutate(
        PopMale = PopMale * 1000
        , PopFemale = PopFemale * 1000
        , PopTotal = PopTotal * 1000
      )
    
    print("UNWPP data read!")
    
    # 1.2. Fertility 
    
    print("processing fertility data...")
    
    # Expand period rates to 5x1 
    
    # UNWPP data is in 5 y age groups
    # Convert to to 1-age
    
    asfr_5_1 <-
      fert_un %>% 
      mutate(ASFR = ASFR / 1000) %>% 
      group_by(LocID, Year) %>%
      group_map( ~ ungroup_fert2(., add_country_year = T), .keep = T) %>% 
      bind_rows() %>% 
      ungroup() %>% 
      filter(!is.na(age)) %>% 
      filter(age %in% reprod_ages) %>% 
      select(LocID, Year, everything())
    
    # Use linear interploation between calendar years to get sinlge-years
    
    asfr_1_1 <-
      asfr_5_1 %>% 
      group_by(LocID, age) %>%
      group_map(
        ~ interpolate_COLUMN_calendar_years(., column = "asfr", years = years_unwpp)
        , .keep = T
      ) %>% 
      bind_rows() %>% 
      ungroup() %>% 
      filter(!is.na(Year)) %>% 
      arrange(LocID, Year, age)
    
    # Convert Period to Cohort 
    
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
      filter(cohort %in% cohorts_unwpp)
    
    print("processing mortality data...")
    
    # 1.3. Mortality 
    
    # UNWPP data is in 5 y age groups
    # Convert to to 1-age
    
    lx_5_1 <-
      lx_un %>% 
      group_by(LocID, Year) %>%
      group_map( 
        ~ lx_fun2(
          data = .
          , age = 0:100
          # , correct_last_age = T
        )
        , .keep = T
      ) %>%
      bind_rows() %>% 
      ungroup() %>% 
      select(LocID, Year, everything())
    
    # Use linear interploation between calendar years to get sinlge-years
    
    lx_1_1 <-
      lx_5_1 %>% 
      group_by(LocID, age) %>%
      group_map(
        ~ interpolate_COLUMN_calendar_years(., column = "lx", years = years_unwpp)
        , .keep = T
      ) %>% 
      bind_rows() %>% 
      ungroup() %>% 
      arrange(LocID, Year, age)
    
    # Get life tables from lx column
    
    lt_1_1_temp <-
      lx_1_1 %>% 
      group_by(LocID, Year) %>%
      group_map( ~ get_lt_bulk_per(data = ., age = 0:100, add_country_year = T), .keep = T) %>%
      bind_rows() %>% 
      ungroup() %>% 
      select(LocID, Year, everything())
    
    if(expand_lt_to_future_stable){
      
      # Expand cohort life tables to the future assuming
      # demographic stability. 
      # This is needed because the function that estimates rates
      # needs to have cohort lifetable for cohorts c:(c+50), where
      # 'c' is the mother's birth cohort. The approximation of cohort 
      # LTs above only produces complete LT for cohorts 1950:2100 (i.e. with 
      # values gor ages 0:100). We have data for older cohorts but it is incomplete 
      # for all ages. For example, LT for cohort 1951 only goes up to age 99 and so on. 
      # One solution is to assume that the values for age 100 in the 1951 cohort LT
      # is equivalent to the last observed values, i.e., the values for age 100 in 
      # the 1950 cohort LT, and so on. 
      
      print(paste0("Expanding period LT beyond 2000 by assuming that values remains constant after period 2100 at the 2100 levels"))
      
      lt_1_1 <- 
        lt_1_1_temp %>% 
        bind_rows(
          lt_1_1_temp %>% 
            group_by(LocID) %>% 
            filter(Year == max(.$Year)) %>% 
            group_map(~repeat_lt(., up_to_year = up_to_year), .keep = T) %>% 
            bind_rows() %>% 
            ungroup()     
        ) %>% 
        arrange(LocID, Year, age)
      
    } else {
      # Keep only cohorts for which wew have 'complete' (approximated) 
      # cohorts values
      lt_1_1 <- 
        lt_1_1_temp %>% 
        arrange(LocID, Year, age)
    }
    
    # Convert Period to Cohort 
    
    lt_coh <- 
      lt_1_1 %>% 
      group_by(LocID) %>%
      arrange(age) %>% 
      group_map( 
        ~ approximate_cohort_from_period(., add_all_combinations = F)
        , .keep = T
      ) %>%
      bind_rows() %>% 
      ungroup() %>% 
      select(-Year) %>% 
      arrange(LocID, cohort, age) %>% 
      filter(cohort %in% cohort_lt)
    
  # Use old cohort rates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  } else if(!reestimate) {
    
    lookup_c_key <- country_key$LocID
    names(lookup_c_key) <- country_key$old_name
    
    print("Reading rates from child death demography paper...")
    
    # Fertility  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    print("Reading fertility data...")
    
    ASFRC <- read.csv(file = paste0("data/input/unwpp/","ASFRC.csv"), stringsAsFactors = F)
    
    asfr_coh <-
      ASFRC %>% 
      mutate(LocID = lookup_c_key[country]) %>% 
      filter(LocID %in% names(lookup_c)) %>%
      filter(Cohort %in% cohorts_unwpp) %>% 
      select(LocID, age = Age, asfr = ASFR, cohort = Cohort) %>% 
      mutate(asfr = asfr / 1000)
    
    asfr_1_1 <- 
      fread(file = paste0("data/input/unwpp/","fert_per_1_1.csv"), stringsAsFactors = F) %>% 
      data.frame() %>% 
      mutate(LocID = lookup_c_key[country]) %>% 
      filter(LocID %in% names(lookup_c)) %>% 
      select(LocID, Year = year, age, asfr = value) %>% 
      mutate(asfr = asfr / 1000)
    
    # Mortality  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    print("Reading mortality data...")
    
    # Female only
    LTCF <- 
      data.table::fread(file = paste0("data/input/unwpp/","LTCF.csv"), stringsAsFactors = F) %>% 
      data.frame()
    
    lt_coh_temp <- 
      LTCF %>% 
      mutate(LocID = lookup_c_key[Country]) %>% 
      filter(LocID %in% names(lookup_c)) %>%
      filter(Cohort %in% cohort_lt) %>% 
      select(LocID, age = Age, nax = ax, lx, nmx = mx, nqx = qx, ndx = dx, nLx, Tx, ex, cohort = Cohort) %>% 
      mutate(across(c(lx, ndx, nLx, Tx), ~ . * 1e5))
    
    lt_1_1 <- 
      fread(file = paste0("data/input/unwpp/","lt_per_1_1_F.csv"), stringsAsFactors = F) %>% 
      data.frame() %>% 
      mutate(LocID = lookup_c_key[country]) %>% 
      filter(LocID %in% names(lookup_c)) %>%
      select(LocID, Year = year, age, nax = ax, lx, nmx = mx, nqx = qx, ndx = dx, nLx = Lx, Tx, ex) %>% 
      # mutate(across(c(lx, ndx, nLx, Tx), ~ . / 1e6))
      mutate(across(c(lx, ndx, nLx, Tx), ~ . / 10))
    
    if(expand_lt_to_future_stable){
      
      # Expand cohort life tables to the future assuming
      # demographic stability. 
      # This is needed because the function that estimates rates
      # needs to have cohort lifetable for cohorts c:(c+50), where
      # 'c' is the mother's birth cohort. The approximation of cohort 
      # LTs above only produces complete LT for cohorts 1950:2100 (i.e. with 
      # values gor ages 0:100). We have data for older cohorts but it is incomplete 
      # for all ages. For example, LT for cohort 1951 only goes up to age 99 and so on. 
      # One solution is to assume that the values for age 100 in the 1951 cohort LT
      # is equivalent to the last observed values, i.e., the values for age 100 in 
      # the 1950 cohort LT, and so on. 
      
      print(paste0("Expanding cohort LT beyond 2000 by assuming that values remains constant after cohort 2000 at the 2000 levels"))
      
      lt_coh <- 
        lt_coh_temp %>% 
        group_by(LocID) %>% 
        group_map(
          ~ expand_cohort_lt_stable(., up_to_year = up_to_year)
          , .keep = T
        ) %>% 
        bind_rows() %>% 
        ungroup() %>% 
        arrange(LocID, cohort, age) %>% 
        data.frame()
      
    } else {
      lt_coh <- lt_coh_temp
    } # end expand_lt_to_future_stable
    
    print("Reading pop data...")
    
    # Pop data
    
    pop_1_1 <- 
      # historical rates
      fread("data/input/unwpp/WPP2019_PopulationBySingleAgeSex_1950-2019.csv", stringsAsFactors = F) %>% 
      data.frame() %>% 
      filter(LocID %in% names(lookup_c)) %>% 
      select(LocID, Year = Time, age = AgeGrp, starts_with("Pop")) %>% 
      # projected rates
      bind_rows(
        fread("data/input/unwpp/WPP2019_PopulationBySingleAgeSex_2020-2100.csv", stringsAsFactors = F) %>% 
          data.frame() %>% 
          filter(LocID %in% names(lookup_c)) %>% 
          select(LocID, Year = Time, age = AgeGrp, starts_with("Pop"))    
      ) %>% 
      mutate(
        PopMale = PopMale * 1000
        , PopFemale = PopFemale * 1000
        , PopTotal = PopTotal * 1000
      )
    
  } # end reestimate
  
  if(export){
    print("Saving ASFR datasets to data/input/asfr_1_1_xxx.csv")
    write.csv(asfr_1_1, "data/input/asfr_1_1_unwpp.csv", row.names = F)
  }
  
  out <- 
    list(
      lt_coh = lt_coh
      , asfr_coh = asfr_coh
      , asfr_1_1 = asfr_1_1
      , lt_1_1 = lt_1_1
      , pop_1_1 = pop_1_1
    )
  
  return(out)
  
}



# Linear interpolation of values to create single-year values
# from 5-y groupped values.
interpolate_COLUMN_calendar_years <- function(df, column = "asfr",  years = 1950:2100) {
  
  y_range_length <- length(years)
  
  val <- approx(df[, column], method = "linear", n = y_range_length)$y    
  
  out <- 
    data.frame(
      LocID = unique(df$LocID)
      , Year = years
      , age = unique(df$age)
    )
  
  out[ , column] <- val
  out
}

kins_fast <- function(ego_age = NULL, year = NULL, # where ego is, it will be at half year
                      P = NULL, asfr = NULL, N = NULL,
                      stable = FALSE,
                      age = 0:100,
                      birth_female = 1/2.04,
                      alive = "yes")
{
  # if stable or not
  if(stable){
    kins <- kins_stable(P = P[,as.character(year)],
                        asfr = asfr[,as.character(year)],
                        birth_female = birth_female) %>%
      filter(x <= ego_age)
  }else{
    kins <- kins_non_stable_fast(ego_age = ego_age, year = year,
                                 P = P, asfr = asfr, N = N,
                                 birth_female = birth_female)
  }
  
  # living results
  alive_yes <- kins %>% filter(alive=="yes")
  alive_yes$alive <- NULL
  kins_by_age_ego <- alive_yes %>% group_by(x) %>% select(-x_kin) %>% summarise_all(sum)
  kins_by_age_kin <- alive_yes[alive_yes$x == ego_age,]
  kins_mean_age   <- colSums(kins_by_age_kin[,3:ncol(alive_yes)]*0:100)/colSums(kins_by_age_kin[,3:ncol(alive_yes)])
  kins_var_age    <- colSums(kins_by_age_kin[,3:ncol(alive_yes)]*(0:100)^2)/colSums(kins_by_age_kin[,3:ncol(alive_yes)]) - kins_mean_age^2
  kins_total      <- colSums(kins_by_age_kin[,c(3:ncol(alive_yes))])
  out_yes <- list(kins = alive_yes,
                  kins_by_age_ego = kins_by_age_ego,
                  kins_by_age_kin = kins_by_age_kin,
                  kins_mean_age = kins_mean_age,
                  kins_total = kins_total)
  
  # death results
  alive_no <- kins %>% filter(alive=="no")
  alive_no$alive <- NULL
  freq_d_by_age_ego <- alive_no %>% group_by(x) %>% select(-x_kin) %>% summarise_all(sum)
  cum_d_by_age_ego  <- freq_d_by_age_ego %>% ungroup() %>%
    summarise_at(.vars = 2:ncol(freq_d_by_age_ego),cumsum) %>%
    mutate(x=freq_d_by_age_ego$x)
  cum_d_total       <- cum_d_by_age_ego %>% filter(x == ego_age) %>% select(-x)
  lost_mean_age     <- colSums(freq_d_by_age_ego[,2:ncol(freq_d_by_age_ego)]*freq_d_by_age_ego$x)/
    colSums(freq_d_by_age_ego[,2:ncol(freq_d_by_age_ego)])
  lost_var_age      <- colSums(freq_d_by_age_ego[,2:ncol(freq_d_by_age_ego)]*freq_d_by_age_ego$x^2)/lost_mean_age^2
  out_no <- list(kins = alive_no,
                 kins_death_by_age_ego = freq_d_by_age_ego,
                 kins_cum_death_by_age_ego = cum_d_by_age_ego,
                 kins_cum_death_total = cum_d_total,
                 lost_mean_age = lost_mean_age,
                 lost_var_age  = lost_var_age)
  
  # not return all
  if(alive=="yes"){
    kins=out_yes
  }else if(alive=="no"){
    kins=out_no
  } else{
    kins=list(kins_living=out_yes,kins_death=out_no)
  }
  
  return(kins)
}

kins_non_stable_fast <- 
  function (ego_age = NULL, year = NULL, P = NULL, asfr = NULL, 
            N = NULL, age = 0:100, birth_female = 1/2.04, Pb = FALSE) {
    
    stopifnot(!is.null(P) & !is.null(asfr) & !is.null(N))
    if (!any(as.integer(colnames(P)) == as.integer(colnames(asfr)))) 
      stop("Data should be from same years.")
    years_data <- as.integer(colnames(P))
    ages <- length(age)
    w <- last(age)
    ego_cohort <- year - ego_age
    zeros = matrix(0, nrow = ages, ncol = ages)
    if (Pb) {
      stopifnot(length(years_data) == ncol(Pb))
    }
    else {
      Pb = P[1, ]
    }
    U = f = list()
    for (t in 1:length(years_data)) {
      Ut = Mt = Dcum = matrix(0, nrow = ages, ncol = ages)
      Ut[row(Ut) - 1 == col(Ut)] <- P[-101, t]
      Ut[ages, ages] = P[101, t]
      diag(Mt) = 1 - P[, t]
      U[[as.character(years_data[t])]] <- rbind(cbind(Ut, 
                                                      zeros), cbind(Mt, Dcum))
      ft = matrix(0, nrow = ages * 2, ncol = ages * 2)
      ft[1, 1:ages] = asfr[, t] * birth_female * (1 + P[, 
                                                        t])/2 * Pb[1, t]
      f[[as.character(years_data[t])]] <- ft
    }
    W <- rbind(t(t(N * asfr)/colSums(N * asfr)), matrix(0, ages, 
                                                        length(years_data)))
    cat(paste0("Rates before ", min(years_data), " assumed as constant."))
    for (y in 1500:(min(years_data) - 1)) {
      U[[as.character(y)]] = U[[as.character(min(years_data))]]
      f[[as.character(y)]] = f[[as.character(min(years_data))]]
      W = cbind(W, W[, as.character(min(years_data))])
      colnames(W)[ncol(W)] = as.character(y)
    }
    if ((year - 1) > max(years_data)) {
      cat(paste0("Rates after ", max(years_data), " assumed as constant."))
      for (y in (max(years_data) + 1):(year - 1)) {
        U[[as.character(y)]] = U[[as.character(max(years_data))]]
        f[[as.character(y)]] = f[[as.character(max(years_data))]]
        W = cbind(W, W[, as.character(max(years_data))])
        colnames(W)[ncol(W)] = as.character(y)
      }
    }
    osM = nosM = oaM = coaM = yaM = cyaM = osM = gmMy = gmM = ggmMy = ggmM = matrix(0, 
                                                                                    ages * 2, ages)
    oaMy = coaMy = yaMy = cyaMy = osMy = nosMy = matrix(0, ages * 
                                                          2, ages)
    e = matrix(0, ages * 2, ages * 2)
    diag(e[1:ages, 1:ages]) = 1
    for (m_age in age) {
      m = W[, as.character(ego_cohort)]
      m_cohort = ego_cohort - m_age - 1
    }
    d =matrix(0, ages * 2, 1)
    m[(ages + 1):(2 * ages)] = 0
    kins = data.frame(x = 0, x_kin = rep(age, 2), alive = c(rep("yes", 
                                                                ages), rep("no", ages)),
                      m = m, d = d)
    
    for (x in 1:ego_age) {
      Ut = U[[as.character(ego_cohort + x - 1)]]
      ft = f[[as.character(ego_cohort + x - 1)]]
      
      m = Ut %*% m
      d = Ut %*% d + ft %*% e[, x]
      
      kins <- rbind(kins, data.frame(x = x, x_kin = rep(age, 
                                                        2), alive = c(rep("yes", ages), rep("no", ages)), 
                                     m = m, d = d))
    }
    return(kins)
  }


lx_fun2 <- function(data, age){
  lx <- splinefun(data$x, data$lx, method = "monoH.FC")(age)
  
  data.frame(
    LocID = unique(data$LocID)
    , Year = unique(data$Year)
    , age = age
    , lx = lx
  )
}



# based on 'set.upperage' function of 
# package 'demography'
# data is a list with x elements
# data$age is a vector of ages
# data$year is a vector of years of interest
# data$pop is a matrix of exposures (from "Exposures_1x1.txt") in which years
# are columns and rows are ages
# data$rate is a similar matrix of nmx values
# Note that data$pop and data$rate can also be lists, each of which 
# contains sex-specific matrices. 
lower_lt_open_age_int <- function(data, max.age){
  
  if (max(data$age) < max.age) 
    stop("max.age too large")
  else if (max(data$age) == max.age) 
    return(data)
  if (is.null(data$pop)) 
    stop("This procedure needs the population data")
  no.rate <- is.null(data$rate)
  multiple.series <- is.list(data$pop)
  if (multiple.series) {
    nn <- length(data$pop)
    fred <- data
    for (j in 1:nn) {
      if (!no.rate) 
        fred$rate <- data$rate[[j]]
      fred$pop <- data$pop[[j]]
      tmp <- set.upperage(fred, max.age)
      if (!no.rate) 
        data$rate[[j]] <- tmp$rate
      data$pop[[j]] <- tmp$pop
    }
    data$age <- tmp$age
  }
  else {
    idx <- data$age >= max.age
    if (sum(!idx) > 0) {
      age <- data$age[!idx]
      rnames <- rownames(data$pop)[!idx]
      if (max(age) < max.age) {
        age <- c(age, max.age)
        rnames <- c(rnames, NA)
      }
    }
    else {
      age <- max.age
      rnames <- ""
    }
    rnames[length(rnames)] <- paste(max.age, "+", sep = "")
    upper.pop <- data$pop[idx, ]
    pop <- apply(matrix(upper.pop, nrow = sum(idx)), 2, sum, 
                 na.rm = TRUE)
    data$pop <- rbind(matrix(data$pop[!idx, ], ncol = ncol(data$pop)), 
                      pop)
    rownames(data$pop) <- rnames
    colnames(data$pop) <- data$year
    if (!no.rate) {
      upper.rate <- data$rate[idx, ]
      actuals = apply(matrix(upper.rate * upper.pop, nrow = sum(idx)), 
                      2, sum, na.rm = TRUE)
      data$rate <- rbind(matrix(data$rate[!idx, ], ncol = ncol(data$rate)), 
                         actuals/pop)
      rownames(data$rate) <- rnames
      colnames(data$rate) <- data$year
    }
    data$age <- age
    if (is.null(data$obs.var)) {
      if (length(dim(data$obs.var)) == 2) 
        data$obs.var <- data$obs.var[data$age <= max.age, 
        ]
      else data$obs.var <- data$obs.var[data$age <= max.age]
    }
  }
  return(data)
}

# Ugo's lifetable function
## function for constructing a classic (& rather general) lifetable
## from mortality rates
lifetable.mx <- function(x, mx, sex="M", ax=NULL, radix = 100000){
  m <- length(x)
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*radix
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, mx, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}

# R function to convert a list of dataframes into a 3-D array
# Obtained from 
# https://www.researchgate.net/publication/277670582_R_function_to_convert_a_list_of_dataframes_into_a_3-D_array
list2ary <- function(input.list){  #input a list of lists
  rows.cols <- dim(input.list[[1]])
  sheets <- length(input.list)
  output.ary <- array(unlist(input.list), dim = c(rows.cols, sheets))
  colnames(output.ary) <- (1:dim(output.ary)[2]) - 1
  row.names(output.ary) <- (1:dim(output.ary)[1]) - 1
  
  return(output.ary)    # output as a 3-D array
}


## Arguments of 'QOSplit' function ## 
## Fx - vector of aggregated ASFRs to be split into single ages ##
## L - vector of ages  (lower limits. e.g 10-15 is 10, 20-24 is 20, and so on) ##
## AgeInt - vector of age intervals ##
QOSplit<-function(Fx,L,AgeInt){
  
  ## (0) Closing left (-99) and right (99) open intervals, if any ##
  n<-length(AgeInt)
  if (AgeInt[1]==-99) {L[1]<-L[1]-AgeInt[2]+1;AgeInt[1]<-AgeInt[2] } 
  if (AgeInt[n]==99) {AgeInt[n]<-AgeInt[n-1]}
  
  h<-AgeInt			# number of 1-years age groups in the 5-years age groups
  N<-sum(h)			# numberof columns  in G
  
  ## (1) Getting G matrix
  G<-matrix(0,nrow=n,ncol=N)
  G[1,1:h[1]]<-1/h[1]			
  j1<-h[1]+1				
  for (i in 2:n){ 
    j2<-j1+h[i]-1	
    G[i,j1:j2]<-1/h[i]			
    j1<-j2+1  } 
  
  # (2) Quadratic minimization
  
  # Matrix for the 2-nd order derivatives
  F<-matrix(0,nrow=N-2,ncol=N)
  for (i in 1:(N-2)){  F[i,i]<-1; F[i,i+1]<--2; F[i,i+2]<-1;}
  G<-G[,-N]; G<-G[,-1];
  F<-F[,-N]; F<-F[,-1];
  Dmat=t(F)%*%F			# matrix for Quadratic form	
  dvec<-rep(0,N-2)
  Amat<-t(rbind(G,diag(N-2)))	# constrains matrix
  meq<-n					# number of equations
  
  bvec<-c(Fx,dvec)		# vector of birth
  aux<-solve.QP(Dmat, dvec, Amat, bvec, meq=meq)
  aux<-c(0,aux$solution,0)	#  the first and the last values are assumed to have no fertility!
  aux[aux<0]<-0   ## extremely small negative values can be neglected and treated as '0's 
  
  solutionN<-round(aux,8)
  Age<-L[1]+c(1:sum(AgeInt))-1
  out<-as.data.frame(cbind(Age,solutionN))
  names(out)<-c("Age","ASFR")
  return(out)  
}

# Get arrays for each relatve type
# The idea is to have the rows represnt ego's age
# columns reresent ego's age, and the third dimension
# be ego's birth cohort
reformat_kin_trees <- function(files, dataset = "kins", relatives = c("d", "m")){
  
  # a. Get data and transform to long format

  df <- lapply(
    files
    , function(f, dataset, relatives){
      readRDS(f)[[dataset]] %>% 
        pivot_longer(-c(x, x_kin), names_to = "kin", values_to = "age_dist") %>% 
        select(kin, age_ego = x, age_kin = x_kin, age_dist) %>% 
        filter(kin %in% relatives) %>% 
        mutate(
          cohort = gsub(".rds", "", f)
          , cohort = str_extract(cohort, "[0-9]{4}$")
        )
    }
    , dataset = dataset
    , relatives = relatives
  ) %>% 
    bind_rows() 
  
  # b. Reformat as 3-dimensional array
  # There should be 1 array per data type
  
  arrays <- 
    lapply(relatives, function(r, df){
      df2 <- 
        df %>% 
        filter(kin == r) 
      
      df_l <- split(
        df2 %>% select(starts_with("age"))   
        , list(df2$cohort)
      )
      
      # Convert to matrix
      df_wide <- 
        lapply(df_l, function(d, ...) {
          out <- 
            d %>%
            pivot_wider(names_from = age_ego, values_from = age_dist) %>% 
            select(-age_kin) 
          # now, this is key: since we want rows to always represent mothers ages
          # and columns to represent children, we have to transpose the matrix
          # when it refers to daughters. The reason for this is that if we don't
          # do it, rows will refer to offspring. To understand why, consider that
          # Demokin works from the pesrpective of an ego. So, when we consider 
          # the mothers of ego, the column 'age_ego' represents the ages of the 
          # 'children' of the mothesr and 'age_kin' are the ages of the mothers.
          # The opposite is true when we consider the children of Ego. In this 
          # case, age-ego are the ages of 'mothers', and age_kin are the ages
          # of the offspring. Therefore, we need to rotate the matrix so that 
          # children are in the columns. However, this does not affect the 
          # values at all, so that the matrix still refers to ego's birth
          # cohort. 
          if(r == "d") out <- data.frame(t(out))
          out
        }, r)
      
      ar <- list2ary(input.list = df_wide)    
      
    }, df)
  
  names(arrays) <- relatives
  
  return(arrays)
}

# Get arrays for each relatve type
# The idea is to have the rows represnt ego's age
# columns reresent ego's age, and the third dimension
# be ego's birth cohort
reformat_kin_trees2 <- function(files, dataset = "kins", relatives = c("d", "m")){
  
  # First, a df with all combinations of cohorts and countries
  
  # Which files to load
  data <- 
    lt_coh %>% 
    distinct(LocID, cohort) %>% 
    mutate(country = lookup_c[as.character(LocID)])
  
  country_keep <- unique(data$country)
  
  base2 <- ifelse(base == "", ".", base)
  full.names <- base != ""
  
  files <- list.files(path = base2, pattern = pat, full.names = full.names)
  
  # A df to identify the DemoKin output files
  
  files_df <- 
    data.frame(files = files) %>% 
    mutate(
      f = gsub(paste0(base, pat, "_"), "", files)
      , f = gsub(".rds", "", f)
      , LocID = str_extract(f, "^[0-9]+")
      , country = lookup_c[LocID]
      , cohort = str_extract(f, "[0-9]{4}$")
    ) %>% 
    select(-f) %>% 
    filter(!is.na(country))
  
  files_n <- 
    files_df %>% 
    filter(LocID %in% lookup_c2[country_keep]) %>% 
    pull(files)
  
  
  # a. Get data and transform to long format
  

  df <- lapply(
    files
    , function(f, dataset, relatives){
      readRDS(f)[[dataset]] %>% 
        pivot_longer(-c(x, x_kin), names_to = "kin", values_to = "age_dist") %>% 
        select(kin, age_ego = x, age_kin = x_kin, age_dist) %>% 
        filter(kin %in% relatives) %>% 
        mutate(
          cohort = gsub(".rds", "", f)
          , cohort = str_extract(cohort, "[0-9]{4}$")
        )
    }
    , dataset = dataset
    , relatives = relatives
  ) %>% 
    bind_rows() 
  
  # b. Reformat as 3-dimensional array
  # There should be 1 array per data type
  
  arrays <- 
    lapply(relatives, function(r, df){
      df2 <- 
        df %>% 
        filter(kin == r) 
      
      df_l <- split(
        df2 %>% select(starts_with("age"))   
        , list(df2$cohort)
      )
      
      # Convert to matrix
      df_wide <- 
        lapply(df_l, function(d, ...) {
          out <- 
            d %>%
            pivot_wider(names_from = age_ego, values_from = age_dist) %>% 
            select(-age_kin) 
          # now, this is key: since we want rows to always represent mothers ages
          # and columns to represent children, we have to transpose the matrix
          # when it refers to daughters. The reason for this is that if we don't
          # do it, rows will refer to offspring. To understand why, consider that
          # Demokin works from the pesrpective of an ego. So, when we consider 
          # the mothers of ego, the column 'age_ego' represents the ages of the 
          # 'children' of the mothesr and 'age_kin' are the ages of the mothers.
          # The opposite is true when we consider the children of Ego. In this 
          # case, age-ego are the ages of 'mothers', and age_kin are the ages
          # of the offspring. Therefore, we need to rotate the matrix so that 
          # children are in the columns. However, this does not affect the 
          # values at all, so that the matrix still refers to ego's birth
          # cohort. 
          if(r == "d") out <- data.frame(t(out))
          out
        }, r)
      
      ar <- list2ary(input.list = df_wide)    
      
    }, df)
  
  names(arrays) <- relatives
  
  return(arrays)
}

# Just replicates a given LT a number of years
repeat_lt <- function(lt, up_to_year = 2150){
  extra_years <- up_to_year - max(lt$Year)
  lt_expanded <-
    lt[ rep(row.names(lt), extra_years), ]
  base_year <- unique(lt$Year)
  new_years <- sort(rep((base_year+1):(base_year+extra_years), length(lt$age)))
  lt_expanded$Year <- new_years
  lt_expanded
}

ungroup_fert2 <- function(df, add_country_year = F) {
  tryCatch(
    expr = {
      out <- 
        c(0, df$ASFR) %>% 
        QOSplit(L = seq(10,45,5), AgeInt=rep(5,8)) %>% 
        rename(age = Age, asfr = ASFR) %>% 
        mutate(age = factor(age, levels = 0:100)) %>% 
        complete(age, fill = list(asfr = 0)) %>% 
        mutate(age = as.numeric(as.character(age))) 
      
      if(add_country_year){
        out$LocID <- unique(df$LocID)  
        out$Year <- unique(df$Year)
      } else {
        out$id <- sample(1:1e6, 1)
      }
      
      out
    }, error = function(e) data.frame(age = NA, asfr = NA, id = NA)
  )
  
}


# 2. Analysis and results------------

# What proportion of a woman\s birth cohort offspring
# died in the mother's lifetime?
get_share_offspring_died <- function(cou, asfr_coh, last_cohort = 2000){
  
  print(cou)
  load(paste0("data/clean/",cou,"_arrays_clean.rdata"))
  
  # Get ASFR
  asfr <- 
    asfr_coh %>% 
    filter(country == cou)
  
  # What % of all the offspring born to mothers in each cohorts
  # died before their mothers?
  
  ## inputs
  ages_mother <- 0:100
  reprod_ages <- 15:50
  ages_child <- 0:(max(ages_mother) - min(reprod_ages))
  
  # Identify correct cohorts
  cohorts <-  min(as.numeric(colnames(D))):last_cohort
  n_row <- length(ages_mother)
  n_col <- length(ages_child)
  n_coh <- length(cohorts)
  
  died_df <- 
    lapply(cohorts, share_offspring_died, H, Q, POP, asfr, cohorts, n_col, n_row, reprod_ages, ages_mother) %>% 
    bind_rows() %>% 
    pivot_longer(-cohort) %>% 
    mutate(country = cou) %>% 
    rename(measure = name)
  
  died_df
}

# What proportion of a woman\s birth cohort offspring
# died in the mother's lifetime?
# Worker function!
share_offspring_died <- function(my.coh, H, Q, POP, asfr, cohorts, n_col, n_row, reprod_ages, ages_mother){
  whi.coh <- which(cohorts==my.coh)
  myH <- H[,1:n_col,whi.coh]
  mypop <- POP[,whi.coh]
  myQ <- matrix(NA,n_row,n_col)
  ## adjusting D 
  for(j in 1:n_col){
    for(i in 1:n_row){
      wrong <- 
        (i - 1) < min(reprod_ages) | 
        (i - 1) > max(ages_mother) | 
        (j - 1) < max(c((i - 1) - max(reprod_ages), 0)) |
        (j - 1) > ((i - 1) - min(reprod_ages))
      if(wrong){
        myQ[i,j] <- NA
      }else{
        myQ[i,j] <- Q[j, whi.coh + i - j]
      }
    } # end rows
  } # end cols
  
  # ALL OFFSPRING DEATHS
  
  m <- ncol(myH)
  MYPOP <- mypop %*% t(rep(1,m))
  ## multiplying matrices
  NU <- myH*MYPOP*myQ
  nu <- rowSums(NU,na.rm = T)
  
  deaths <- sum(nu)
  
  # ALL OFFSPRING BIRTHS
  # 20220311 This is wrong!
  # The previous approaches multiplied excepted children surviving by populatino 
  # of mothers. We should actually multiply population of mothers by fertility 
  # rates to get estimate of how many offspring will bo born to a cohort of mothers
  # CB <- myH*MYPOP
  # cb <- rowSums(CB,na.rm = T)
  # births <- sum(cb)
  
  fert <- 
    asfr %>% 
    filter(cohort == my.coh) %>% 
    pull(asfr)
  
  fert <- c(rep(0, 15), fert, rep(0, 50))
  births <- sum(mypop * fert)
  
  prop_died <- deaths/births
  
  return(
    data.frame(
      cohort = my.coh
      , deaths = deaths
      , births = births
      , prop_died = prop_died
    )
  )  
}


get_maol_iqr <- function(H, D, L, cohorts_kin = 1836:1950, by_cohs = 1){
  # 0. Functions
  
  # These are functions used below
  make_df_from_array <- function(a, cohorts_kin, by_cohs = 1){
    
    l <- 
      apply(a, 3, function(mat){
        mat %>% 
          data.frame() %>% 
          rownames_to_column("age_kin") %>% 
          pivot_longer(-age_kin,names_to = "age_ego", values_to = "value") %>% 
          mutate(
            age_ego = as.numeric(gsub("[A-Z]+", "", age_ego)) - 1
            , age_kin = as.numeric(age_kin) - 1
          )
      }) 
    
    cohorts_to_add <- cohorts_kin[seq(1, length(cohorts_kin), by = by_cohs)]
    
    dfs <-
      lapply(seq_along(l), function(n, ...){
        df <- l[[n]]
        df$cohort <- cohorts_to_add[n]
        df
      }, cohorts_to_add) %>%
      bind_rows()
    
    return(dfs)  
  }
  
  ## Ugo's function to compute median age at death using splines
  cont_med.fun <- function(x,y){
    fit <- spline(x=x, y=y, n=10000) 
    xi <- fit$x 
    yi <- fit$y 
    if (any(round(yi,3)==0.5)){
      q2 <- xi[which(round(yi,3)==0.5)]
    }else if (any(round(yi,3)==0.501)){
      q2 <- xi[which(round(yi,3)==0.501)]
    }else if (any(round(yi,3)==0.499)){
      q2 <- xi[which(round(yi,3)==0.499)]
    }
    med <- q2[1]
    return(med)
  }
  
  ## Ugo's function to compute IQR median age at death using splines
  cont_iqr.fun <- function(x,y, max_age = 100){
    
    appr <- min(y) > 0.251
    
    fit <- spline(x=x, y=y, n=10000)
    xi <- fit$x
    yi <- fit$y
    if (any(round(yi,3)==0.25)){
      q1 <- xi[which(round(yi,3)==0.25)]
    }else if (any(round(yi,3)==0.251)){
      q1 <- xi[which(round(yi,3)==0.251)]
    }else if (any(round(yi,3)==0.249)){
      q1 <- xi[which(round(yi,3)==0.249)]
    } else if(appr){
      q1 <- max_age
    }
    if (any(round(yi,3)==0.75)){
      q3 <- xi[which(round(yi,3)==0.75)]
    }else if (any(round(yi,3)==0.751)){
      q3 <- xi[which(round(yi,3)==0.751)]
    }else if (any(round(yi,3)==0.749)){
      q3 <- xi[which(round(yi,3)==0.749)]
    }
    iqr <- abs(q1[1] - q3[1])
    
    iqr <- ifelse(appr, paste0(iqr, "_flagged"), iqr)
    iqr <- as.character(iqr)
    return(iqr)
  }
  
  
  # Pre-defined parameters for Sweden
  
  ages_mother <- 0:100
  reprod_ages <- 15:50
  ages_child <- 0:(max(ages_mother) - min(reprod_ages))
  # cohorts_kin <- 1836:1950
  # by_cohs <- 1
  
  # Create objects for loop
  n_row <- length(ages_mother)
  n_col <- length(ages_child)
  
  # A hack so that we can run it for a single cohort
  
  cohs <- seq(1, length(cohorts_kin), by = by_cohs)  
  
  h <- d <- l <- num_arr <- wrong <- NULL
  hdl <- numeric(0)
  num_arr <- array(NA, dim = c(n_row, n_col, 1))
  
  # Loop over cohorts
  for(coh in cohs){
    # Loop over columns
    for(j in 1:n_col){
      # Loop over rows
      for(i in 1:n_row){
        
        # See if we should just skip a combination of maternal and child ages
        # From the paper: 
        # We restrict the range of possible maternal and offspring ages $(i,j)$ to avoid impossible cases such as mothers aged 65 experiencing the death of an infant. 
        # For mothers, we consider ages $\alpha \leq i \leq 100$, where $[\alpha,\beta] = [15,49]$ are the limits of female reproductive life.
        # Offspring ages are restricted to $max(i-\beta,0) \leq j \leq i - \alpha$, where $max$ is a function that returns the maximum of two values. 
        
        wrong <- 
          (i - 1) < min(reprod_ages) | 
          (i - 1) > max(ages_mother) | 
          (j - 1) < max(c((i - 1) - max(reprod_ages), 0)) |
          (j - 1) > ((i - 1) - min(reprod_ages))

        
        if(wrong) {
          hdl <- c(hdl, NA)
        } else {
          # Numerator ~~~~~~~~~~~~~~~
          # h_{i,j,c} * d_{j,c+i-j} * l_{i,c}
          
          h <- H[i, j, coh]
          d <- D[j, 1, coh + i - j]
          l <- L[i, 1, coh]
          
          hdl <- c(hdl, h*d*l)
          
        } # end else
      } # end rows
    } # end cols    
    
    # convert back to matrix for given cohort
    num <- matrix(hdl, ncol = n_col, byrow = F) 
    # den <- matrix(bp, ncol = n_col, byrow = F) 
    
    #   # Add to array
    num_arr <-
      array(
        c(num_arr, num)
        , dim = c(n_row, n_col, dim(num_arr)[3] + 1)
      )
    
    num <-  NULL
    hdl <- numeric(0)
    
  } # end cohorts loop
  
  # Remove extra first dimension in array
  num_arr <- num_arr[, , -1]
  
  # If only one cohort, it will be a matrix, but should be array
  if(length(dim(num_arr)) == 2){
    num_arr <- array(num_arr, dim = c(n_row, n_col, 1))
  }
  
  ### 3. Get Median/IQR age at child loss 
  
  # Here, we use Ugo's continuous approach to get the median and IQR. Namely using a spline on the lx function and derive the quartiles with greater precision.
  # There are a very small number of cases where mortality is so low that more than 25% of the population are still alive at age 100+. 
  # For these cases, we assume that Q1 is 100. These cases are flagged as being approximated (this underestimates the IQR). 
  # The idea is to flag out these cases in a potential plot. 
  
  # convert the arrays to a df in long format:
  
  num_df <- make_df_from_array(num_arr, cohorts_kin, by_cohs = by_cohs)
  
  # Normalise numeratort to sum to 1
  nu <- 
    num_df %>% 
    select(cohort, age_kin, age_ego, value) %>%
    # Group by mother's age
    group_by(cohort, age_kin) %>% 
    summarise(value = sum(value, na.rm = T)) %>% 
    ungroup() %>% 
    # Normalise numerator to add to 1 within each country/year
    group_by(cohort) %>%
    mutate(numerator = value/sum(value)) %>% 
    ungroup() %>% 
    rename(age = age_kin) 
  
  # Get IQR and median using Ugo's continuous approach
  median_iqr_ol <-
    nu %>% 
    select(cohort, age, var = numerator) %>% 
    arrange(cohort, age) %>% 
    # filter(country == "Mali", cohort == 1950) %>% 
    group_by(cohort) %>% 
    mutate(cum = cumsum(var)) %>% 
    summarise(
      median = cont_med.fun(x = age, y = cum) 
      , iqr = cont_iqr.fun(x = age, y = cum, max_age = max(ages_mother)) 
    ) %>%
    ungroup() 
  
  return(median_iqr_ol)
  
}

get_maol_all <- function(my.coh, H, deaths_mat , surv_mat, what = "median"){
  
  whi.coh <- which(cohorts==my.coh)
  surv_vec <- surv_mat[,whi.coh]
  myH <- H[,1:n_col,whi.coh]
  DEATHS <- matrix(NA,n_row,n_col)
  ## adjusting the deaths matrix (Q or D)
  for(j in 1:n_col){
    for(i in 1:n_row){
      wrong <- 
        (i - 1) < min(reprod_ages) | 
        (i - 1) > max(ages_mother) | 
        (j - 1) < max(c((i - 1) - max(reprod_ages), 0)) |
        (j - 1) > ((i - 1) - min(reprod_ages))

      if(wrong){
        DEATHS[i,j] <- NA
      }else{
        DEATHS[i,j] <- deaths_mat[j, whi.coh + i - j]
      }
    } # end rows
  } # end cols
  
  out <- worker_MAOL(myH = myH, surv_vec = surv_vec, DEATHS = DEATHS, ages_mother = ages_mother, what = what)
  return(out)
}

# Inspired by get_maol_all, except that it lists all the arguments (get_maol_all just
# inhertis them from the upper function environment)
get_maol_all2 <- function(my.coh, H, deaths_mat , surv_mat, what = "median", cohorts, n_col, n_row, reprod_ages, ages_mother){
  
  whi.coh <- which(cohorts==my.coh)
  surv_vec <- surv_mat[,whi.coh]
  myH <- H[,1:n_col,whi.coh]
  DEATHS <- matrix(NA,n_row,n_col)
  ## adjusting the deaths matrix (Q or D)
  for(j in 1:n_col){
    for(i in 1:n_row){
      wrong <- 
        (i - 1) < min(reprod_ages) | 
        (i - 1) > max(ages_mother) | 
        (j - 1) < max(c((i - 1) - max(reprod_ages), 0)) |
        (j - 1) > ((i - 1) - min(reprod_ages))

      if(wrong){
        DEATHS[i,j] <- NA
      }else{
        DEATHS[i,j] <- deaths_mat[j, whi.coh + i - j]
      }
    } # end rows
  } # end cols
  
  out <- worker_MAOL(myH = myH, surv_vec = surv_vec, DEATHS = DEATHS, ages_mother = ages_mother, what = what)
  return(out)
}

# Use this to get estimates for many countries using lapply
get_maol_all_batch <- function(iso){
  
  # Load data
  # iso <- countrycode(country_name, "country.name", "iso3c")
  print(iso)
  
  f <- paste0("data/clean/",iso,"_arrays_clean.rdata")
  load(f)
  
  # 0. Parameters 
  
  ages_mother <- 0:100
  reprod_ages <- 15:50
  ages_child <- 0:(max(ages_mother) - min(reprod_ages))
  
  if(iso == "SWE"){
    cohorts_kin <-  1765:2000
  } else if (iso == "DNK"){
    cohorts_kin <- 1878:2000  
  } else {
    # UNWPP
    cohorts_kin <- 1950:2000
  }
  
  cohorts <- cohorts_kin
  n_row <- length(ages_mother)
  n_col <- length(ages_child)
  n_coh <- length(cohorts)
  n <- length(reprod_ages)
  
  # 1. All-age MAOL and IQR 
  
  # nu = h_{i,j}*pop_i*q_j
  
  maol_all <- 
    unlist(
      lapply(
        cohorts_kin, 
        get_maol_all2
        , H = H
        , deaths_mat = Q
        , surv_mat = POP
        , what = "mean"
        , cohorts = cohorts, n_col = n_col, n_row = n_row
        , reprod_ages = reprod_ages, ages_mother = ages_mother
      )
    )
  
  
  sd_all <- 
    unlist(
      lapply(
        cohorts_kin, get_maol_all2
        , H = H, deaths_mat = Q, surv_mat = POP, what = "sd"
        , cohorts = cohorts, n_col = n_col, n_row = n_row
        , reprod_ages = reprod_ages, ages_mother = ages_mother
      )
    )
  
  out <- 
    data.frame(
      cohort = cohorts_kin
      , country = iso
      , M = maol_all
      , SD = sd_all
    ) %>% 
    pivot_longer(c(M, SD), names_to = "measure")
  
  out  
}

get_maol_by_age <- function(my.coh, H, deaths_mat , surv_mat, what = "median", rescale_nu = T){
  
  whi.coh <- which(cohorts==my.coh)
  surv_vec <- surv_mat[,whi.coh]
  myH <- H[,1:n_col,whi.coh]
  DEATHS <- matrix(NA,n_row,n_col)
  ## adjusting the deaths matrix (Q or D)
  for(j in 1:n_col){
    for(i in 1:n_row){
      wrong <- 
        (i - 1) < min(reprod_ages) | 
        (i - 1) > max(ages_mother) | 
        (j - 1) < max(c((i - 1) - max(reprod_ages), 0)) |
        (j - 1) > ((i - 1) - min(reprod_ages))

      if(wrong){
        DEATHS[i,j] <- NA
      }else{
        DEATHS[i,j] <- deaths_mat[j, whi.coh + i - j]
      }
    } # end rows
  } # end cols
  
  out <- worker_MAOL_1(myH = myH, surv_vec = surv_vec, DEATHS = DEATHS, ages_mother = ages_mother, what = what, rescale_nu = rescale_nu)
  out
}

# Transform matrix Z to a lonf data frame
# on a country by country basis
get_curves <- function(cou, last_cohort = 2000){
  print(cou)
  load(paste0("results/decomp/",cou,"_SSE.Rdata"))
  
  # Identify correct cohorts
  first_cohort <- last_cohort - ncol(Z) + 1
  lookup_coh <- as.character(first_cohort:last_cohort)
  names(lookup_coh) <- paste0("X", seq_along(first_cohort:last_cohort))
  
  counts <- 
    data.frame(Z) %>% 
    mutate(age = 15:100) %>% 
    pivot_longer(-age, names_to = "cohort", values_to = "value") %>% 
    mutate(
      cohort = lookup_coh[cohort]
      , country = cou
      , measure = "count"
    )
  
  # Summary data
  mu_hat <- 
    data.frame(fitSSE$MU.hat) %>% 
    mutate(age = 15:100) %>% 
    pivot_longer(-age, names_to = "cohort", values_to = "value") %>% 
    mutate(
      cohort = lookup_coh[cohort]
      , country = cou
      , measure = "MU.hat"
    )
  
  gamma1_hat <- 
    data.frame(fitSSE$GAMMA1.hat) %>% 
    mutate(age = fitSSE$x1) %>%
    pivot_longer(-age, names_to = "cohort", values_to = "value") %>% 
    mutate(
      cohort = lookup_coh[cohort]
      , country = cou
      , measure = "GAMMA1.hat"
    )
  
  gamma2_hat <- 
    data.frame(fitSSE$GAMMA2.hat) %>% 
    mutate(age = fitSSE$x2) %>%
    pivot_longer(-age, names_to = "cohort", values_to = "value") %>% 
    mutate(
      cohort = lookup_coh[cohort]
      , country = cou
      , measure = "GAMMA2.hat"
    )
  
  
  # Components
  comp <- 
    data.frame(
      age = NA
      , cohort = lookup_coh
      , country = cou
      , w1 = SSEsummary$w1
      , M1 = SSEsummary$M1
      , M2 = SSEsummary$M2 
      , SD1  = SSEsummary$SD1
      , SD2 = SSEsummary$SD2
    ) %>% 
    pivot_longer(-c(age, cohort, country), names_to = "measure", values_to = "value")
  
  out <- bind_rows(counts, mu_hat, gamma1_hat, gamma2_hat, comp)  
  
  return(out)
  
}


get_cd_by_age_batch <- function(iso){

  # Load data
  print(iso)
  
  f <- paste0("data/clean/",iso,"_arrays_clean.rdata")
  load(f)
  
  # 0. Parameters 
  
  ages_mother <- 0:100
  reprod_ages <- 15:50
  ages_child <- 0:(max(ages_mother) - min(reprod_ages))
  
  if(iso == "SWE"){
    cohorts_kin <-  1765:2000
  } else if (iso == "DNK"){
    cohorts_kin <- 1878:2000  
  } else {
    # UNWPP
    cohorts_kin <- 1950:2000
  }
  
  cohorts <- cohorts_kin
  n_row <- length(ages_mother)
  n_col <- length(ages_child)
  n_coh <- length(cohorts)
  n <- length(reprod_ages)
  
  # 1. All-age MAOL and IQR 
  
  # nu = h_{i,j}*pop_i*q_j
  
  out <- 
    lapply(
      cohorts_kin
      , get_cd_by_age_batch2
      , H = H
      , deaths_mat = Q
      , surv_mat = POP
      , cohorts = cohorts, n_col = n_col, n_row = n_row
      , reprod_ages = reprod_ages, ages_mother = ages_mother
    ) %>% 
    bind_rows() %>% 
    mutate(country = iso)
  
  out  
}

get_cd_by_age_batch_counterfactual <- function(iso){
  
  # Load data
  print(iso)
  
  f <- paste0("data/clean/",iso,"_arrays_clean.rdata")
  load(f)
  
  # 0. Parameters 
  
  ages_mother <- 0:100
  reprod_ages <- 15:50
  ages_child <- 0:(max(ages_mother) - min(reprod_ages))
  
  if(iso == "SWE"){
    cohorts_kin <-  1765:2000
  } else if (iso == "DNK"){
    cohorts_kin <- 1878:2000  
  } else {
    # UNWPP
    cohorts_kin <- 1950:2000
  }
  
  cohorts <- cohorts_kin
  n_row <- length(ages_mother)
  n_col <- length(ages_child)
  n_coh <- length(cohorts)
  n <- length(reprod_ages)
  
  # 1. All-age MAOL and IQR 
  
  # nu = h_{i,j}*pop_i*q_j
  
  out <- 
    lapply(
      cohorts_kin
      , get_cd_by_age_batch2_counterfactual
      , H = H
      , deaths_mat = Q
      , surv_mat = POP
      , cohorts = cohorts, n_col = n_col, n_row = n_row
      , reprod_ages = reprod_ages, ages_mother = ages_mother
    ) %>% 
    bind_rows() %>% 
    mutate(country = iso)
  
  out  
}

get_cd_by_age_batch2 <- function(my.coh, H, deaths_mat , surv_mat, cohorts, n_col, n_row, reprod_ages, ages_mother){

  whi.coh <- which(cohorts==my.coh)
  surv_vec <- surv_mat[,whi.coh]
  myH <- H[,1:n_col,whi.coh]
  DEATHS <- matrix(NA,n_row,n_col)
  ## adjusting the deaths matrix (Q or D)
  for(j in 1:n_col){
    for(i in 1:n_row){
      wrong <- 
        (i - 1) < min(reprod_ages) | 
        (i - 1) > max(ages_mother) | 
        (j - 1) < max(c((i - 1) - max(reprod_ages), 0)) |
        (j - 1) > ((i - 1) - min(reprod_ages))

      if(wrong){
        DEATHS[i,j] <- NA
      }else{
        DEATHS[i,j] <- deaths_mat[j, whi.coh + i - j]
      }
    } # end rows
  } # end cols
  
  out <- 
    worker_cd_by_age(myH = myH, surv_vec = surv_vec, DEATHS = DEATHS, ages_mother = ages_mother) %>% 
    mutate(cohort = my.coh)
  
  return(out)
}

get_cd_by_age_batch2_counterfactual <- function(my.coh, H, deaths_mat , surv_mat, cohorts, n_col, n_row, reprod_ages, ages_mother){
  
  whi.coh <- which(cohorts==my.coh)
  surv_vec <- surv_mat[,whi.coh]
  myH <- H[,1:n_col,whi.coh]
  DEATHS <- matrix(NA,n_row,n_col)
  ## adjusting the deaths matrix (Q or D)
  for(j in 1:n_col){
    for(i in 1:n_row){
      wrong <- 
        (i - 1) < min(reprod_ages) | 
        (i - 1) > max(ages_mother) | 
        (j - 1) < max(c((i - 1) - max(reprod_ages), 0)) |
        (j - 1) > ((i - 1) - min(reprod_ages))
      
      if(wrong){
        DEATHS[i,j] <- NA
      }else{
        DEATHS[i,j] <- deaths_mat[1, whi.coh + i - j]
      }
    } # end rows
  } # end cols
  
  out <- 
    worker_cd_by_age(myH = myH, surv_vec = surv_vec, DEATHS = DEATHS, ages_mother = ages_mother) %>% 
    mutate(cohort = my.coh)
  
  return(out)
}

worker_cd_by_age <- function(myH,surv_vec,DEATHS,ages_mother){
  ## making L a matrix
  m <- ncol(myH)
  SURV <- surv_vec %*% t(rep(1,m))
  ## multiplying matrices
  NU <- myH*SURV*DEATHS
  
  
  out <- 
    NU %>% 
    data.frame() %>% 
    rownames_to_column("mother_age") %>% 
    pivot_longer(-mother_age, names_to = "child_age") %>% 
    mutate(
      child_age = as.numeric(gsub("X", "", child_age))
      , mother_age = as.numeric(mother_age)
    ) 
  
  return(out)
}

# the arugments death_vec and a survival_vec depends on the 
# estimations method (hdl or hpopq)
worker_MAOL <- function(myH,surv_vec,DEATHS,ages_mother, what = "median"){
  ## making L a matrix
  m <- ncol(myH)
  SURV <- surv_vec %*% t(rep(1,m))
  ## multiplying matrices
  NU <- myH*SURV*DEATHS
  ## summing rows and scale to 1
  nu <- rowSums(NU,na.rm = T)
  nus <- nu/sum(nu)
  ## cumulative sum
  nu_cum <- cumsum(nus) 
  ## computing the median
  fit <- spline(x=ages_mother, y=nu_cum, n=10000) 
  xi <- fit$x 
  yi <- fit$y 
  
  if("median" %in% what){
    if (any(round(yi,3)==0.5)){
      q2 <- xi[which(round(yi,3)==0.5)]
    }else if (any(round(yi,3)==0.501)){
      q2 <- xi[which(round(yi,3)==0.501)]
    }else if (any(round(yi,3)==0.499)){
      q2 <- xi[which(round(yi,3)==0.499)]
    }
    median <- q2[1]  
  }
  if("iqr" %in% what){
    if (any(round(yi,3)==0.25)){
      q1 <- xi[which(round(yi,3)==0.25)]
    }else if (any(round(yi,3)==0.251)){
      q1 <- xi[which(round(yi,3)==0.251)]
    }else if (any(round(yi,3)==0.249)){
      q1 <- xi[which(round(yi,3)==0.249)]
    } else if(appr){
      q1 <- max_age
    }
    if (any(round(yi,3)==0.75)){
      q3 <- xi[which(round(yi,3)==0.75)]
    }else if (any(round(yi,3)==0.751)){
      q3 <- xi[which(round(yi,3)==0.751)]
    }else if (any(round(yi,3)==0.749)){
      q3 <- xi[which(round(yi,3)==0.749)]
    }
    iqr <- abs(q1[1] - q3[1])
  }
  if("mean" %in% what){
    ## mean and sd
    ax <- 0.5
    n <- length(ages_mother)
    x1 <- ages_mother
    
    mean <- sum((x1+ax)*nu)/sum(nu)
  }
  if("sd" %in% what){
    ## mean and sd
    ax <- 0.5
    n <- length(ages_mother)
    x1 <- ages_mother
    
    mean <- sum((x1+ax)*nu)/sum(nu)
    v <- sum(nu * ((x1+ax) - mean)^2) / sum(nu)
    sd <- sqrt(v)
  }
  
  out <- get(what)
  return(out)
}

worker_MAOL_1 <- function(myH,surv_vec,DEATHS,ages_mother, what = "median", rescale_nu = T){
  ## making L a matrix
  m <- ncol(myH)
  SURV <- surv_vec %*% t(rep(1,m))
  ## multiplying matrices
  NU <- myH*SURV*DEATHS
  
  # Do for each child age
  out_l <- 
    lapply(as.data.frame(NU), function(col, ...){
      nu <- col[!is.na(col)]
      target_ages <- which(!is.na(col)) - 1
      ## summing rows and scale to 1
      if(rescale_nu) nu <- nu/sum(nu)
      ## cumulative sum
      nu_cum <- cumsum(nu) 
      
      ## computing the median
      Pinv <- approxfun(nu_cum, target_ages, yleft=NA, yright=NA)
      # Pinv <- approxfun(nu_cum, 1:length(nu), yleft=NA, yright=NA)
      if("median" %in% what) out <- Pinv(0.5)
      if("iqr" %in% what) out <- abs(Pinv(0.75) - Pinv(0.25))
      out
    }, rescale_nu)
  
  unlist(out_l)
  
}
