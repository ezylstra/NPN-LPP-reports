# ABQ BioPark Botanic Garden
# Erin Zylstra
# erinz@usanpn.org
# 30 January 2026

require(rnpn)
require(tidyr)
require(dplyr)
require(lubridate)
require(stringr)
require(ggplot2)
require(ggforce) # Nice facet options for ggplot
require(terra)
require(lme4)
require(ggeffects) # Predictions & plotting for mixed-effect models
library(emmeans)
# require(lmerTest) # Don't want to load this (but will call in script)

source("geom_abs_text.R")

# require(geosphere) # for optional site clustering
# require(fpc) # for optional site clustering
# require(RColorBrewer)

#- Specify data of interest ---------------------------------------------------#
# Name and nickname of LPP
lpp <- "ABQ BioPark Botanic Garden"
lpp_short <- "abq"

# Specify years of interest
all_yrs <- 2009:2025

#- Load and format status/intensity data --------------------------------------#

# Check if the status/intensity data have already been downloaded and saved. 
# If not, do so now
si_csv_name <- paste0("npn-data/si-", lpp_short, "-thru", max(all_yrs), ".csv")

if (!file.exists(si_csv_name)) {
  # Name of person requesting NPN data
  requestor <- "erinz"
  
  # Get list of LPP network ID(s)
  network_ids <- npn_groups()
  network_id <- network_ids %>%
    filter(str_detect(network_name, lpp)) %>%
    filter(str_detect(network_name, "deprecated", negate = TRUE))
  if (nrow(network_id) == 1) {
    message("Good news: One network_id found for LPP")
    lpp_id <- network_id$network_id
  } else if (nrow(network_id) == 0) {
    stop("No matches for LPP name", call. = FALSE)
  } else {
    warning("Multiple network_ids found for LPP", call. = FALSE)
    lpp_id <- c(network_id$network_id)
  }
  
  # Download status/intensity data
  si <- npn_download_status_data(
    request_source = requestor,
    network_ids = lpp_id,
    years = all_yrs,
    climate_data = FALSE,
    additional_fields = c("site_name", "observedby_person_id")
  )

  # Remove some fields that we don't need (to make file smaller)
  si <- si %>%
    select(-c(update_datetime, species))
  
  # Write to file
  write.csv(si, si_csv_name, row.names = FALSE)
} 

# Load csv
si <- read.csv(si_csv_name)

# Format date and observer
si <- si %>%
  mutate(obsdate = ymd(observation_date)) %>%
  mutate(yr = year(observation_date)) %>%
  select(-observation_date) %>%
  rename(person_id = observedby_person_id)

# See if there's a decent amount of information each year
count(si, yr) # Will remove 2013 (only 5 observations)
si <- si %>% filter(yr > 2013)

# Remove any duplicate records (all entries the same except observation_id)
si <- si %>% distinct(across(-observation_id))

# Append phenophase "groups" (csv created in phenophases-intensities.R)
pheno_list <- read.csv("phenophases.csv")
si <- left_join(si, select(pheno_list, phenophase_id, class_id, pheno_group),
                by = "phenophase_id")

# Append intensity single values (csv created in phenophase-intensities.R)
intensity_list <- read.csv("intensities.csv")
# Going to keep intensity values numeric, but want to differentiate between the
# different types (number, percent, qualitative). Will convert percentages to 
# proportions and qualitative values to huge numbers.
intensity_list <- intensity_list %>%
  rename(value_orig = value) %>%
  mutate(value = case_when(
    type == "percent" ~ value_orig / 100,
    type == "qualitative" ~ value_orig * 100000,
    .default = value_orig
  ))
intensity_list <- intensity_list %>%
  select(category_id, value_name, type, value) %>%
  rename(intensity_category_id = category_id, 
         intensity_value = value_name, 
         intensity_type = type,
         intensity = value)
si <- left_join(si, intensity_list, 
                by = c("intensity_category_id", "intensity_value"))

#- Extract site, species information ------------------------------------------#
# Removing to keep working objects smaller. Can merge site, species information 
# back in at any point

sites <- si %>%
  group_by(site_id, site_name, latitude, longitude, 
           elevation_in_meters, state) %>%
  summarize(n_plant_spp = length(unique(species_id[kingdom == "Plantae"])),
            n_plants = length(unique(individual_id[kingdom == "Plantae"])),
            n_animal_spp = length(unique(species_id[kingdom == "Animalia"])),
            n_observers = length(unique(person_id)),
            n_yrs = length(unique(yr)),
            yr_first = min(yr),
            yr_last = max(yr),
            .groups = "keep") %>%
  data.frame()

species <- si %>%
  group_by(kingdom, species_id, genus, common_name) %>%
  summarize(n_sites = length(unique(site_id)),
            n_individuals = length(unique(individual_id)),
            n_observers = length(unique(person_id)),
            n_yrs = length(unique(yr)),
            yr_first = min(yr),
            yr_last = max(yr),
            .groups = "keep") %>%
  data.frame()

# Add some information into species dataframe
ids <- species$species_id
spp_info <- npn_species() %>%
  data.frame() %>%
  filter(species_id %in% ids) %>%
  select(species_id, species, functional_type, family_common_name, family_name)
species <- species %>%
  left_join(spp_info, by = "species_id")

# Remove columns we no longer need
si <- si %>%
  select(-c(site_name, latitude, longitude, elevation_in_meters, state,
            species_id, genus, abundance_value))

#- Clean up phenophase status data --------------------------------------------#
# Any phenophase_status == -1 or 0 when intensity value present?
count(si, phenophase_status, is.na(intensity_value))
# Apparently, this does happen (at least for -1)

# Change unknown phenophase_status to 1 if an intensity value reported
# and then remove any remaining records with phenophase status unknown
si <- si %>%
  mutate(phenophase_status = case_when(
    !is.na(intensity_value) ~ 1,
    .default = phenophase_status
  )) %>%
  filter(phenophase_status %in% c(0, 1))

#- Deal with multiple observations by the same person -------------------------#
# Same plant or animal (same site), same phenophase, same date - same observer

# For now, will keep record with more advanced phenophase or higher 
# intensity value. Will do this by sorting observations in descending
# order and keeping only the first
inddateobsp <- si %>%
  group_by(common_name, individual_id, obsdate, person_id, phenophase_id) %>%
  summarize(n_obs = n(),
            .groups = "keep") %>%
  data.frame()
inddateobsp$obsnum <- 1:nrow(inddateobsp)

si <- si %>%
  arrange(person_id, individual_id, obsdate, phenophase_id, 
          desc(phenophase_status), desc(intensity)) %>%
  left_join(select(inddateobsp, -c(n_obs, common_name)), 
            by = c("person_id", "individual_id", 
                   "obsdate", "phenophase_id")) %>%
  # Create "dups" column, where dups > 1 indicates that the observation can be
  # removed since there's another observation that same day with more advanced
  # phenology or higher intensity/abundance.
  mutate(dups = sequence(rle(as.character(obsnum))$lengths))

# Remove extra observations and unnecessary columns
si <- si %>%
  filter(dups == 1) %>%
  select(-c(obsnum, dups)) %>%
  arrange(common_name, obsdate, person_id, phenophase_id)

#- Accounting for multiple observations by different people -------------------#
# Same plant or animal, same phenophase, same date - different observers

# For plants, multiple observations of the same individual and phenophase on the 
# same date can tell us about inter-observer variation (assuming phenophase 
# status and intensity doesn't change).

# To summarize sampling effort (e.g., when and where observations were made), 
# probably best to keep all the observations in the dataset.

# To analyze phenophase onsets, probably best
# to subset or aggregate observations so there is a maximum of one record for 
# each phenophase associated with a plant each day.

# Function to calculate sum, but if all values in vector = NA, then return NA
na_max <- function(x) {
  y <- ifelse(sum(is.na(x)) == length(x), NA, max(x, na.rm = TRUE))
  return(y)
}

#- Creating datasets for plant observations -----------------------------------#

# For an individual plant, we have some choices about how to characterize
# the phenological state and intensity (if recorded) when multiple observations
# were made: 
# 1) retain one observation of the plant per day.
# 2) retain all data, use status = 1 observations, and for intensity, average 
#    over values

si_sub_plants <- si %>%
  filter(kingdom == "Plantae") %>%
  arrange(common_name, individual_id, obsdate, phenophase_id, 
          desc(phenophase_status), desc(intensity)) 

  # For each species, summarize number of phenophases, pheno classes, pheno 
  # groups, and observations per pheno class
    spp_pc <- si_sub_plants %>%
      group_by(common_name) %>%
      summarize(n_phenophases = length(unique(phenophase_id)),
                n_phenoclasses = length(unique(class_id)),
                n_phenogroups = length(unique(pheno_group)),
                initial_leaf = sum(class_id == 1),
                young_leaf = sum(class_id == 2),
                leaf = sum(class_id == 3),
                color_leaf = sum(class_id == 4),
                fall_leaf = sum(class_id == 5),
                flower = sum(class_id == 6),
                open_flower = sum(class_id == 7),
                pollen = sum(class_id == 8),
                fruit = sum(class_id == 10),
                ripe = sum(class_id == 12),
                drop_fruit = sum(class_id == 13)) %>%
      data.frame()
    spp_pc

# For now, will go with Option 1 above (selecting one observation instead of 
# averaging information among observations)
# Sequential process for selecting a single observation:
# 1) select the complete observation (ideally, status and intensity values for
#    multiple phenophases).
# 2) for observations with an equal ranking of information, select the one 
#    from the more consistent observer of that plant in that year.
    
# To do this we'll need to:
# 1) Put data in wide form with each row having all information from an
#    observer of a plant in that day.
# 2) Fix inconsistencies within that observer-plant-date combo
# Then we can start to see how often we have conflicting information from 
# multiple observers

# Putting data in wide form:
# First, will remove any full flowering observations for lilac (if they're in 
# there). Then can use phenophase class_id instead of phenophase_id because 
# there's a maximum of one phenophase per class for each species.

# Create a table with info/names for each plant phenophase class    
pheno_list <- pheno_list %>%
  mutate(class_short = case_when(
    class_id == 1 ~ "initial_growth",
    class_id == 2 ~ "young_leaf",    
    class_id == 3 ~ "leaf",
    class_id == 4 ~ "color_leaf",
    class_id == 5 ~ "fall_leaf",
    class_id == 6 ~ "flower",
    class_id == 7 ~ "open_flower",
    class_id == 8 ~ "pollen",
    class_id == 10 ~ "fruit",
    class_id == 12 ~ "ripe_fruit",
    class_id == 13 ~ "drop_fruit",
    .default = NA
  )) %>%
  mutate(group_short = case_when(
    pheno_group == "Leaves" ~ "leaf",
    pheno_group == "Leaf senescence" ~ "leaf_end",
    pheno_group == "Flowers" ~ "flower",
    pheno_group == "Open flowers" ~ "flower_open",
    pheno_group == "Fruits" ~ "fruit",
    pheno_group == "Ripe fruits" ~ "fruit_ripe"
  ))

pl_pheno_classes <- pheno_list %>%
  filter(class_id %in% 1:13) %>%
  select(class_id, class_short, group_short) %>%
  distinct() %>%
  mutate(group_code = case_when(
    group_short == "leaf" ~ "lf",
    group_short == "leaf_end" ~ "lfe",
    group_short == "flower" ~ "fl",
    group_short == "flower_open" ~ "flo",
    group_short == "fruit" ~ "fr",
    group_short == "fruit_ripe" ~ "frr"
  )) %>%
  filter(!is.na(class_short)) %>%
  arrange(class_id) %>%
  mutate(class_id2 = str_pad(class_id, width = 2, pad = 0))

# Summary table for phenogroups:
pheno_list <- pheno_list %>%
  mutate(pheno_group = replace(pheno_group, 
                               pheno_group == "Leaf senescence", 
                               "Colored leaves"))
pgs_pl <- pheno_list %>%
  filter(class_id %in% pl_pheno_classes$class_id) %>%
  group_by(pheno_group, class_id, class_name) %>%
  summarize(n_phenophases = n(), .groups = "keep") %>%
  data.frame() %>%
  arrange(class_id)

# Put in wide form
plant_obs <- si_sub_plants %>%
  filter(phenophase_id != 206) %>%
  left_join(select(pl_pheno_classes, class_id, class_id2), 
            by = "class_id") %>%
  select(-c(phenophase_id, phenophase_description, pheno_group, intensity_value, 
            intensity_category_id, intensity_type, class_id)) %>%
  rename(status = phenophase_status) %>%
  pivot_wider(names_from = class_id2,
              values_from = c(status, intensity),
              names_sort = TRUE) %>%
  data.frame()
# Add in columns for any pheno classes that aren't already in there
all_cols <- c(paste0("status_", pl_pheno_classes$class_id2),
              paste0("intensity_", pl_pheno_classes$class_id2))
missing_cols <- setdiff(all_cols, colnames(plant_obs))
if (length(missing_cols) > 0) {
  add <- data.frame(matrix(ncol = length(missing_cols), nrow = 1, NA))
  colnames(add) <- missing_cols
  add <- add %>%
    mutate(across(starts_with("status"), as.numeric)) %>%
    mutate(across(starts_with("intensity"), as.integer))
  plant_obs <- cbind.data.frame(plant_obs, add)
}
# Put columns in order
first_status_col <- which(grepl("status", colnames(plant_obs)))[1]
plant_obs <- plant_obs %>%
  select(colnames(plant_obs)[1:(first_status_col - 1)], all_of(all_cols))

# Resolve any inconsistencies with statuses of phenophase classes
  # If young leaves = 1, leaves = 1 (ok if young = 1, leaves = NA)
  # If color leaves = 1, leaves = 1
  # If open flower = 1, flower = 1
  # If pollen = 1, flower = 1
  # If ripe fruit = 1, fruit = 1
# For now, leaving instances where more specific category = 1 and general
# category (leaves, flower, fruit) = NA

# To make this a little easier to follow, will temporarily rename columns
pl_pheno_classes <- pl_pheno_classes %>%
  mutate(stat_cols_c = paste0("status_", class_short),
         int_cols_c = paste0("intensity_", class_short),
         stat_cols_g = paste0("status_", group_code),
         int_cols_g = paste0("intensity_", group_code))
plant_obs <- plant_obs %>%
  rename_with(~ pl_pheno_classes$stat_cols_c,
              all_of(paste0("status_", pl_pheno_classes$class_id2))) %>%
  rename_with(~ pl_pheno_classes$int_cols_c,
              all_of(paste0("intensity_", pl_pheno_classes$class_id2)))  

# First, see how often these issues come up?
  count(plant_obs, status_young_leaf, status_leaf)
  count(plant_obs, status_color_leaf, status_leaf)
  count(plant_obs, status_open_flower, status_flower)
  count(plant_obs, status_pollen, status_flower)
  count(plant_obs, status_ripe_fruit, status_fruit)

# Replace values where needed
plant_obs <- plant_obs %>%
  mutate(status_leaf = case_when(
    (status_young_leaf == 1 & status_leaf == 0) ~ 1,
    (status_color_leaf == 1 & status_leaf == 0) ~ 1,
    .default = status_leaf
  )) %>%
  mutate(status_flower = case_when(
    (status_open_flower == 1 & status_flower == 0) ~ 1,
    (status_pollen == 1 & status_flower == 0) ~ 1,
    .default = status_flower
  )) %>%
  mutate(status_fruit = case_when(
    (status_ripe_fruit == 1 & status_fruit == 0) ~ 1, 
    .default = status_fruit
  ))

# Summarize data available for each observer-plant-date combination
plant_obs <- plant_obs %>%
  mutate(n_status = rowSums(!is.na(pick(starts_with("status"))))) %>%
  mutate(n_yes = rowSums(pick(starts_with("status")), na.rm =TRUE)) %>%
  mutate(n_intensities = rowSums(!is.na(pick(starts_with("intensity")))))

# For each plant and year, create an observer rank based on the number of dates
# each observer monitored the plant (higher rank [lower number] goes to 
# the most consistent observer)
  
  # Find all plant-year combinations where at least once, multiple people made 
  # observations on the same day
  mult_obs_date <- plant_obs %>%
    group_by(individual_id, obsdate, yr) %>%
    summarize(n_observations = n(), .groups = "keep") %>%
    data.frame() %>%
    mutate(obsnum = row_number())
  mult_obs_yr <- mult_obs_date %>%  
    group_by(individual_id, yr) %>%
    summarize(mult_obs = 1 * any(n_observations > 1), .groups = "keep") %>%
    data.frame() %>%
    filter(mult_obs == 1)
  # Rank observers by number of days monitored, breaking ties with number of 
  # intensities recorded and number of statuses recorded
  for (i in 1:nrow(mult_obs_yr)) {
    plant_obs_sub <- plant_obs %>%
      filter(individual_id == mult_obs_yr$individual_id[i]) %>%
      filter(yr == mult_obs_yr$yr[i])
    observer_count <- plant_obs_sub %>%
      group_by(person_id) %>%
      summarize(n_obs = n(),
                n_status = sum(n_status),
                n_intensities = sum(n_intensities)) %>%
      data.frame() %>%
      arrange(desc(n_obs), desc(n_intensities), desc(n_status)) %>%
      mutate(observer_rank = row_number()) %>%
      mutate(individual_id = mult_obs_yr$individual_id[i], 
             yr = mult_obs_yr$yr[i]) 
    if (i == 1) {
      observer_rank <- observer_count
    } else {
      observer_rank <- rbind(observer_rank, observer_count)
    }
  }
  
# Attach number of observations and observer rank to dataframe
plant_obs <- plant_obs %>%
  left_join(mult_obs_date, by = c("individual_id", "obsdate", "yr")) %>%
  left_join(select(observer_rank, person_id, individual_id, yr, observer_rank),
            by = c("person_id", "individual_id", "yr"))   
# Check that there's always an observer rank if there are multiple 
# observations of that plant on that date (should not be any NA values):
summary(plant_obs$observer_rank[plant_obs$n_observations > 1])
# Check that if there's no observer rank, there's only one observation of that
# plant on that date (should all be 1):
summary(plant_obs$n_observations[is.na(plant_obs$observer_rank)])

# Finally, remove all but one observation of a plant per day. Selecting based
# on: 
# 1) Number of statuses reported
# 2) Number of "yes" statuses with an intensity value reported
# 3) Observer rank (preferentially selecting data from the person who observed
#    the plant more often than others)

plant_obs2 <- plant_obs %>%
  # Replance NAs in observer rank with 0 for easier sorting
  mutate(observer_rank = replace_na(observer_rank, 0)) %>%
  # Sort dataframe by above criteria
  arrange(obsnum, desc(n_status), desc(n_intensities), observer_rank) %>%
  # Remove all but one observation and remove unnecessary columns
  mutate(dups = sequence(rle(as.character(obsnum))$lengths)) %>%
  filter(dups == 1) %>%
  select(-c(obsnum, observer_rank, dups))

# Summary of objects we'll use moving forward:
  # si = dataframe where each row has data associated with unique combination of
    # site, species, individual_id, date, phenophase, and observer. Best for 
    # summarizing effort. May contain multiple observations of an individual and
    # phenophase in a day.
  # plant_obs = dataframe where each row has data for a unique combination of
    # site, species, individual_id, date, and observer. Best for summarizing 
    # effort. May contain multiple observations of a plant and phenophase in a 
    # day.
  # plant_obs2 = dataframe where each row has data for a unique combination of
    # site, species, individual_id, and date. Best for summarizing plant
    # phenology. No more than one observation of a plant each day.

#- Start creating data overviews / effort summaries ---------------------------#

#- Table: Site summaries ------------------------------------------------------#

# Site name, no. plant species, no. individual plants, no. animal species, 
# no. observers, year range, no. observations (plant/animal-date-observer)
sitest <- sites %>%
  select(site_id, site_name, n_plant_spp, n_plants, n_animal_spp, n_observers, 
         n_yrs, yr_first, yr_last) %>%
  rename(plant_spp = n_plant_spp,
         plants = n_plants,
         animal_spp = n_animal_spp,
         observers = n_observers,
         years = n_yrs)
# Find number of observations of plants
site_p <- plant_obs %>%
  group_by(site_id) %>%
  summarize(n_plant_obs = n()) %>%
  data.frame()
# Combine and clean up
sitest <- sitest %>%
  left_join(site_p, by = "site_id") %>%
  rename(obs_plant = n_plant_obs) %>%
  arrange(desc(years), desc(obs_plant))				 

sitest_sub <- sitest %>%
  select(-c(site_id, yr_first, yr_last)) %>%
  mutate(site_name = str_trim(site_name, side = "both"))

# write.table(sitest_sub, "clipboard", sep = "\t", row.names = FALSE)

#- Table: Yearly summaries ----------------------------------------------------#

# No. of sites, No. of observers, no. plant observations by year
yearst <- si %>%
  group_by(yr) %>%
  summarize(sites = length(unique(site_id)),
            observers = length(unique(person_id))) %>% 
  data.frame()
# Find number of observations of plants
year_p <- plant_obs %>%
  group_by(yr) %>%
  summarize(obs_plant = n()) %>%
  data.frame()
# Combine and clean up
yearst <- yearst %>%
  left_join(year_p, by = "yr")

# write.table(yearst, "clipboard", sep = "\t", row.names = FALSE)

#- Table: Species summaries ---------------------------------------------------#
# Recreating species table after removing years with few observations

yrs_subset <- yearst$yr[yearst$obs_plant > 25]

species2 <- si %>%
  filter(yr %in% yrs_subset) %>%
  group_by(kingdom, common_name) %>%
  summarize(n_sites = length(unique(site_id)),
            n_individuals = length(unique(individual_id)),
            n_observers = length(unique(person_id)),
            n_yrs = length(unique(yr)),
            yr_first = min(yr),
            yr_last = max(yr),
            .groups = "keep") %>%
  data.frame() %>%
  left_join(select(species, common_name, functional_type), by = "common_name")

# Use the plant_obs2 dataframes to calculate frequency (number of days
# between observations within a given year). Note that we're not differentiating
# between phenophases here. 
freq_pl <- plant_obs2 %>%
  select(site_id, common_name, individual_id, day_of_year, yr) %>%
  filter(yr %in% yrs_subset) %>%
  arrange(common_name, individual_id, yr, day_of_year)
freq_pl$interval <- NA
for (i in 2:nrow(freq_pl)) {
  if (freq_pl$individual_id[i] == freq_pl$individual_id[i - 1] &
      freq_pl$yr[i] == freq_pl$yr[i - 1]) {
    freq_pl$interval[i] <- freq_pl$day_of_year[i] - freq_pl$day_of_year[i - 1]
  } 
}
freqs_pl <- freq_pl %>%
  group_by(common_name, individual_id, yr) %>%
  summarize(nobs = n(), 
            interval = median(interval, na.rm = TRUE),
            .groups = "keep") %>%
  data.frame() %>%
  group_by(common_name) %>%
  summarize(nobs_mn = round(mean(nobs)),
            interval_mn = round(mean(interval, na.rm = TRUE), 1)) %>%
  data.frame()

freqs <- freqs_pl %>%
  rename(nobs_mn_per_yrind = nobs_mn,
         interval_mn_per_yrind = interval_mn)
species2 <- species2 %>%
  left_join(freqs, by = "common_name") %>%
  mutate(functional_type = factor(functional_type, 
                                  levels = c("Cactus",
                                             "Deciduous broadleaf",
                                             "Drought deciduous broadleaf",
                                             "Evergreen broadleaf",
                                             "Forb"))) %>%
  arrange(desc(kingdom), functional_type, common_name, .locale = "en")

# write.table(species2, "clipboard", sep = "\t", row.names = FALSE)

# Summary stats for plants species
  # Total number of plants
  (nplants <- sum(species2$n_individuals[species2$kingdom == "Plantae"]))
  # Mean number of plants per species
  nplants / sum(species2$kingdom == "Plantae")
  # Mean number of observations per plant per year (mean of species means)
  (obsperplant <- mean(species2$nobs_mn_per_yrind[species2$kingdom == "Plantae"]))
  # Mean observation interval (mean of species means)
  (obsinterval <- mean(species2$interval_mn_per_yrind[species2$kingdom == "Plantae"]))
  
# Date range within year
  count(plant_obs2, day_of_year)

#- Map(s) ---------------------------------------------------------------------#

# Go to site-maps.R

#- Aggregate status information within phenophase group -----------------------#

# Create table summarizing phenophase class groups
count(si, class_id, pheno_group, phenophase_description) %>%
  left_join(distinct(pheno_list, class_id, class_name))
  
# Aggregate status information within pheno group (doesn't make sense to 
# aggregate intensity values, since they may not always be the same type)

# Create new data frame without intensity data and other data summaries
plg_status <- plant_obs2 %>% 
  select(-c(contains("intens"), "n_status", "n_yes", "n_observations"))

# Aggregate information across classes within each pheno group
for (group in unique(pl_pheno_classes$group_code)) {
  cols <- pl_pheno_classes$stat_cols_c[pl_pheno_classes$group_code == group]
  plg_status[,paste0("sum_", group)] <- apply(as.matrix(plg_status[,cols]), 
                                              1, na_max)
}
# Remove columns associated with classes and rename columns with status data
plg_status <- plg_status %>% select(-contains("status_"))
colnames(plg_status) <- str_replace(colnames(plg_status), "sum", "status")

#- Look at weekly proportions -------------------------------------------------#

df <- plg_status

# Add week to the data so we can calculate weekly proportions
# We'll also create wk_doy columns to assign each week with a day of the year:
# wk_doy1 = start of each week (eg, date for week 1 would be Jan 1)
# wk_doy4 = middle of each week (eg, date for week 1 would be Jan 4)
df <- df %>%
  mutate(wk = week(obsdate)) %>%
  # Remove observations in week 53
  filter(wk < 53) %>%
  # Create wk_doy columns
  mutate(wk_date1 = parse_date_time(paste(2024, wk, 1, sep = "/"), "Y/W/w"),
         wk_date1 =  as.Date(wk_date1),
         wk_doy1 = yday(wk_date1),
         wk_date4 = parse_date_time(paste(2024, wk, 4, sep = "/"), "Y/W/w"),
         wk_date4 =  as.Date(wk_date4),
         wk_doy4 = yday(wk_date4))

# Keep just one observation of each plant and phenophase, each week. Sort so the 
# most advanced phenophase gets kept (if more than one value in a week)
df_lf <- df %>%
  arrange(common_name, individual_id, yr, wk, desc(status_lf)) %>%
  distinct(individual_id, yr, wk, .keep_all = TRUE) %>%
  filter(!is.na(status_lf)) %>%
  mutate(phenogroup = "lf") %>%
  rename(status = status_lf) %>%
  select(-contains("status_"))
df_fl <- df %>%
  arrange(common_name, individual_id, yr, wk, desc(status_fl)) %>%
  distinct(individual_id, yr, wk, .keep_all = TRUE) %>%
  filter(!is.na(status_fl)) %>%
  mutate(phenogroup = "fl") %>%
  rename(status = status_fl) %>%
  select(-contains("status_"))
df_flo <- df %>%
  arrange(common_name, individual_id, yr, wk, desc(status_flo)) %>%
  distinct(individual_id, yr, wk, .keep_all = TRUE) %>%
  filter(!is.na(status_flo)) %>%
  mutate(phenogroup = "flo") %>%
  rename(status = status_flo) %>%
  select(-contains("status_"))
df_fr <- df %>%
  arrange(common_name, individual_id, yr, wk, desc(status_fr)) %>%
  distinct(individual_id, yr, wk, .keep_all = TRUE) %>%
  filter(!is.na(status_fr)) %>%
  mutate(phenogroup = "fr") %>%
  rename(status = status_fr) %>%
  select(-contains("status_"))
df_frr <- df %>%
  arrange(common_name, individual_id, yr, wk, desc(status_frr)) %>%
  distinct(individual_id, yr, wk, .keep_all = TRUE) %>%
  filter(!is.na(status_frr)) %>%
  mutate(phenogroup = "frr") %>%
  rename(status = status_frr) %>%
  select(-contains("status_"))
df_lfe <- df %>%
  arrange(common_name, individual_id, yr, wk, desc(status_lfe)) %>%
  distinct(individual_id, yr, wk, .keep_all = TRUE) %>%
  filter(!is.na(status_lfe)) %>%
  mutate(phenogroup = "lfe") %>%
  rename(status = status_lfe) %>%
  select(-contains("status_"))
df <- bind_rows(df_lf, df_fl, df_flo, df_fr, df_frr, df_lfe)

# Exclude species monitored in only 1-2 years (2024-2025; winterfat, 
# gilia beardtongue, and indianhemp)
df <- df %>%
  filter(!common_name %in% c("winterfat",
                             "gilia beardtongue",
                             "indianhemp"))

# Calculate weekly proportions of yeses for each species, phenophase
wkprops <- df %>%
  group_by(common_name, phenogroup, wk, wk_date1, wk_doy1, 
           wk_date4, wk_doy4) %>%
  summarize(n_obs = n(),
            n_yes = sum(status),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(prop = n_yes / n_obs)
# 35 spp-phpg-wk combinations have only 4 observations (the 2 milkweeds),
# but otherwise all have 5 or more, so don't need to exclude more.

# Visualize weekly proportions
# # Leaves
# ggplot(filter(wkprops, phenogroup == "lf"),
#        aes(x = wk, y = prop)) +
#   geom_line() +
#   geom_point(aes(size = n_obs)) +
#   facet_wrap(~common_name, ncol = 2) +
#   geom_vline(xintercept = 40, color = "steelblue3") +
#   labs(title = "Leaves")
# # Leaf end
# ggplot(filter(wkprops, phenogroup == "lfe"),
#        aes(x = wk, y = prop)) +
#   geom_line() +
#   geom_point(aes(size = n_obs)) +
#   facet_wrap(~common_name, ncol = 2) +
#   geom_vline(xintercept = 40, color = "steelblue3") +
#   labs(title = "Leaf end")
# # Flower
# ggplot(filter(wkprops, phenogroup == "fl"),
#        aes(x = wk, y = prop)) +
#   geom_line() +
#   geom_point(aes(size = n_obs)) +
#   facet_wrap(~common_name, ncol = 2) +
#   geom_vline(xintercept = 40, color = "steelblue3") +
#   labs(title = "Flowers")
# # Flower open
# ggplot(filter(wkprops, phenogroup == "flo"),
#        aes(x = wk, y = prop)) +
#   geom_line() +
#   geom_point(aes(size = n_obs)) +
#   facet_wrap(~common_name, ncol = 2) +
#   geom_vline(xintercept = 40, color = "steelblue3") +
#   labs(title = "Flowers open")
# # Fruit
# ggplot(filter(wkprops, phenogroup == "fr"),
#        aes(x = wk, y = prop)) +
#   geom_line() +
#   geom_point(aes(size = n_obs)) +
#   facet_wrap(~common_name, ncol = 2) +
#   geom_vline(xintercept = 40, color = "steelblue3") +
#   labs(title = "Fruit")
# # Fruit ripe
# ggplot(filter(wkprops, phenogroup == "frr"),
#        aes(x = wk, y = prop)) +
#   geom_line() +
#   geom_point(aes(size = n_obs)) +
#   facet_wrap(~common_name, ncol = 2) +
#   geom_vline(xintercept = 40, color = "steelblue3") +
#   labs(title = "Fruit ripe")

# Leaves: calendar looks fine (some with leaves all yr: saltbush, rabbitbrush)
# Leaf end: calendar looks fine (rabbitbrush low values all year)
# Flowers: calendar looks better than water EXCEPT for Siberian elm (and maybe 
  # stretchberry), which start at very end of year
# Flowers open: calendar looks fine

# Fruit: 
  # no seasonal signal for horsetail milkweed (0s), stretchberry (~0.25)
  # calendar best for eastern cottonwood, golden currant, Siberian elm
  # 3 species have inverse patterns (high most of year with dip mid-year), so
    # calendar year would work: fourwing saltbush, screwbean mesquite, tree cholla. 
  # summer year best for rabbitbush, with 0's only ~ wks 25-38
  # showy milkweed has 0's in weeks 10-20....
# Fruit ripe similarly mixed
  # horsetail always 0 so REMOVE

# Create plots with  smooths. then it will be easier to figure out how we 
# want to deal with figures displaying onsets and estimating trends when 
# not all are good fits with a calendar year.

# Create plots with weekly proportions and smooths curves between -------------#

# Not going to worry about year for now - lumping all years together
# Could use a GAM to model weekly proportions and have the ends of the year
# match up (and get predictions with CIs), but this is something we'll probably 
# have to come back. For now will plot raw data (with smoothed spline)

# To create figures with ticks between month labels on the x axis, need to do 
# a little extra work. 
# Dates where we want month labels (15th of month)
x_lab <- as.Date(paste0("2024-", 1:12, "-15"))
# Dates where we want ticks (1st of month)
x_tick <- as.Date(c(paste0("2024-", 1:12, "-01"), "2024-12-31"))
# Will specify axis breaks & ticks at 1st and 15th of month. Make labels on
# the 1st black and change color of tick marks on the 15th to NA.

x_lab_doy <- yday(x_lab)
x_tick_doy <- yday(x_tick)
n_x_tick <- length(x_tick)
month_labels <- month.abb[month(x_lab)]

# Create a column in onsets df with nice names of phenophase groups, for plots
wkprops <- wkprops %>%
  mutate(group_labels = case_when(
    phenogroup == "lf" ~ "Leaves",
    phenogroup == "lfe" ~ "Colored leaves",
    phenogroup == "fl" ~ "Flowers",
    phenogroup == "flo" ~ "Open flowers",
    phenogroup == "fr" ~ "Fruit",
    phenogroup == "frr" ~ "Ripe fruit"
  )) %>%
  mutate(group_labels = factor(group_labels,
                               levels = c("Colored leaves",
                                          "Ripe fruit",
                                          "Fruit",
                                          "Open flowers",
                                          "Flowers",
                                          "Leaves"))) %>%
  mutate(phenogroup_f = factor(group_labels, 
                               levels = rev(levels(group_labels))))

# Assign colors for each phenophase group
color_vec <- c("#b2df8a",   # Leaves
               "#fdbf6f",   # Flowers
               "#ff7f00",   # Open flowers
               "#cab2d6",   # Fruits
               "#6a3d9a",   # Ripe fruits
               "#8c510a")   # Colored leaves
# Name vector so colors are consistent across figures (in case not all levels 
# are present in all figures)
names(color_vec) <- levels(wkprops$phenogroup_f)

# Make green for leaves a little darker in color
rgb2col = function(rgbmat){
  ProcessColumn = function(col){
    rgb(rgbmat[1, col], 
        rgbmat[2, col], 
        rgbmat[3, col], 
        maxColorValue = 255)
  }
  sapply(1:ncol(rgbmat), ProcessColumn)
}
darker_gr <- rgb2col(round(col2rgb(color_vec["Leaves"]) * 0.8))
color_vec["Leaves"] <- darker_gr

wkprops <- wkprops %>% 
  left_join(select(species, common_name, functional_type), 
            by = "common_name") %>%
  mutate(func_group = functional_type) %>%
  mutate(common_name_full = paste0(common_name, " (", 
                                   str_to_lower(functional_type), ")")) %>%
  arrange(func_group, functional_type, common_name)

# Look at horsetail milkweed fruit data
hmf <- wkprops %>%
  filter(common_name == "horsetail milkweed") %>%
  filter(phenogroup %in% c("fr", "frr"))
hmf %>%
  group_by(phenogroup) %>%
  summarize(n_wks = n_distinct(wk),
            prop_g0 = sum(prop > 0),
            max_prop = max(prop)) %>%
  data.frame()
hmf %>%
  filter(prop > 0)

# Remove fruit, fruit ripe data for horsetail milkweed (almost all 0s)
wkprops <- wkprops %>%
  filter(!(common_name == "horsetail milkweed" & phenogroup %in% c("fr", "frr")))

spp_list <- unique(wkprops$common_name_full)

# 4 species in each multi-panel plot
n_plots <- ceiling(length(spp_list) / 4)

pt_alpha <- 0.3

for (i in 1:n_plots) {
  spps <- spp_list[(i * 4 - 3):(i * 4)]
  n_spps <- sum(!is.na(spps))
  props4 <- filter(wkprops, common_name_full %in% spps)
  
  # Get smoothed spline for each species, phenophase
  for (k in 1:n_spps) {
    props <- filter(props4, common_name_full == spps[k])
    pg_list <- unique(props$phenogroup)
    
    for (j in pg_list) {
      
      props_pg <- props %>% filter(phenogroup == j)
      pg_preds <- data.frame(spline(x = props_pg$wk_doy4, 
                                    y = props_pg$prop,
                                    n = 1000)) %>%
        mutate(y01 = if_else(y < 0, 0, y)) %>%
        mutate(y01 = if_else(y01 > 1, 1, y01)) %>%
        mutate(common_name_full = spps[k]) %>%
        mutate(phenogroup_f = props_pg$phenogroup_f[1])
      
      if (j == pg_list[1]) {
        pg_preds_all <- pg_preds
      } else {
        pg_preds_all <- rbind(pg_preds_all, pg_preds)
      }
    }
    if (k == 1) {
      pg_preds4 <- pg_preds_all
    } else {
      pg_preds4 <- rbind(pg_preds4, pg_preds_all)
    }
  }
  
  assign(paste0("phenophase_periods_p", i),
         ggplot() +
           geom_point(data = props4,
                      aes(x = wk_doy4, y = prop,
                          size = n_obs, color = phenogroup_f), alpha = pt_alpha) +
           scale_size_continuous(range = c(0.5, 4), limits = c(0, 40),
                                 breaks = seq(10, 40, by = 10)) +
           geom_line(data = pg_preds4, 
                     aes(x = x, y = y01, color = phenogroup_f)) +
           scale_color_manual(values = color_vec) +
           scale_x_continuous(limits = c(min(x_tick_doy), max(x_tick_doy)),
                              # expand = c(0, 0),
                              breaks = c(x_lab_doy, x_tick_doy),
                              labels = c(month_labels, rep("", n_x_tick))) +
           scale_y_continuous(expand = c(0.04, 0.04), 
                              breaks = seq(0, 1, by = 0.2)) +
           facet_wrap(~ factor(common_name_full, levels = spps), ncol = 1) +
           labs(y = "Proportion of observations in phenophase", 
                color = "", size = "No. observations") +
           theme_bw() +
           theme(axis.ticks.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                       rep("black", n_x_tick))),
                 panel.grid.major = element_line(color = "gray95"),
                 panel.grid.minor = element_blank(), 
                 legend.position = "bottom",
                 legend.box = "vertical", 
                 legend.margin = margin(),
                 axis.text.x = element_text(size = 9), 
                 axis.text.y = element_text(size = 9),
                 legend.title = element_text(size = 9),
                 axis.title.x = element_blank()) +
           guides(size = guide_legend(order = 1),
                  color = guide_legend(nrow = 1, byrow = TRUE, order = 2))
  )
  
  # ggsave(paste0("output/phenophase-periods-spp-", lpp_short, "-", i, ".png"),
  #        get(paste0("phenophase_periods_p", i)),
  #        width = 6.5,
  #        height = ifelse(n_spps == 4, 7.5, 4.25),
  #        units = "in",
  #        dpi = 600)
}

# Create functions to calculate day of wateryr and summeryr -------------------#

wateryr_calc = function(x, start.month = 10){
  x = as.Date(x)
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1)
  as.integer(x - start.date + 1)
}
summeryr_calc = function(x, start.month = 7){
  x = as.Date(x)
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1)
  as.integer(x - start.date + 1)
}

# Create boxplots for onset dates ---------------------------------------------#

# Function to help add sample sizes to plot (95% of the way across date axis)
n_fun <- function(x) {
  return(data.frame(y = 0.99 * datelims[2],
                    label = length(x)))
}

# Create a unique ID for each plant-year to make various matches easier
df2 <- df %>%
  mutate(plantyr = paste0(individual_id, "_", yr)) %>%
  rename(doy = day_of_year)

# Create simplified dataframe, with water/summer year options
df2 <- df %>%
  rename(doy = day_of_year) %>%
  select(person_id, site_id, common_name, individual_id, phenogroup, status, 
         obsdate, yr, doy) %>%
  mutate(obsmonth = month(obsdate)) %>%
  mutate(wateryr = ifelse(obsmonth %in% 10:12,
                          paste0(yr, "-", yr + 1),
                          paste0(yr -1, "-", yr))) %>%
  mutate(summeryr = ifelse(obsmonth %in% 6:12,
                           paste0(yr, "-", yr + 1),
                           paste0(yr - 1, "-", yr))) %>%
  mutate(dowy = wateryr_calc(obsdate),
         dosy = summeryr_calc(obsdate))

# Steps:
# Isolate status data for each species-plant-phenophase
# Calculate prior no
# Keep only yes observations with prior no <= 14
# Group observations by year-period (may not always be calendar)
# Retain date of earliest yes

# Will do this by species to make it easier (though code will be longer)

# tree cholla ----------------------------#
  # using Jan-Dec and skipping fruit phps (no signal)
  trch <- df2 %>%
    filter(common_name == "tree cholla") %>%
    filter(phenogroup %in% c("fl", "flo")) %>%
    arrange(individual_id, phenogroup, obsdate) %>%
    mutate(prior_no = NA)
  # Calculate prior no
  for (i in 2:nrow(trch)) {
    trch$prior_no[i] <- ifelse(
      trch$individual_id[i] == trch$individual_id[i-1] & 
        trch$phenogroup[i] == trch$phenogroup[i-1],
      as.numeric(trch$obsdate[i] - trch$obsdate[i-1]),
      NA)
  }
  # Exclude observations with prior no > 14
  trch <- trch %>%
    filter(status == 1) %>%
    filter(!is.na(prior_no) & prior_no < 15)
  # Group by year-period
  trch_gr <- trch %>%
    group_by(common_name, individual_id, phenogroup, yr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = min(doy),
              .groups = "keep") %>%
    mutate(group_levels = case_when(
      phenogroup == "lf" ~ "Leaves",
      phenogroup == "fl" ~ "Flowers",
      phenogroup == "flo" ~ "Open flowers",
      phenogroup == "fr" ~ "Fruits",
      phenogroup == "frr" ~ "Ripe fruits",
      phenogroup == "lfe" ~ "Colored leaves"
      )) %>%
    mutate(group_levels = factor(group_levels,
                                 levels = c("Colored leaves",
                                            "Ripe fruits",
                                            "Fruits",
                                            "Open flowers",
                                            "Flowers",
                                            "Leaves"))) %>%
    data.frame()

  # Assign colors for each phenophase group
  color_vec <- c("#b2df8a",   # Leaves
                 "#fdbf6f",   # Flowers
                 "#ff7f00",   # Open flowers
                 "#cab2d6",   # Fruits
                 "#6a3d9a",   # Ripe fruits
                 "#8c510a")   # Colored leaves
  # Name vector so colors are consistent across figures (in case not all levels 
  # are present in all figures)
  names(color_vec) <- rev(levels(trch_gr$group_levels))
  
  # Create nice labels for dates on X axes in plots
  # mindoy <- min(onsets_sub$firstyes)
  # maxdoy <- max(onsets_sub$firstyes)
  doys <- seq(1, 365)
  plotdates <- as.Date(paste(2023, doys, sep = "-"), "%Y-%j")
  date1ind <- which(day(plotdates) == 1)
  doybreaks <- doys[date1ind]
  datebreaks <- plotdates[date1ind] %>% format("%m/%d")
  datelims <- c(min(doys), max(doys) + 10)

  boxplot_trch <- ggplot(data = trch_gr,
                         aes(x = group_levels, 
                             y = first_yes_doy, 
                             fill = group_levels)) +
    geom_boxplot() +
    stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, 
                 family = "sans", size = 3) +
    scale_y_continuous(limits = datelims, breaks = doybreaks, 
                       labels = datebreaks) +
    scale_fill_manual(values = color_vec) +
    facet_grid(~common_name) +
    labs(x = "", y = "First day observed") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9),
          strip.text = element_text(size = 9))
  
  # ggsave("output/onsets-plot-abq-trch.png",
  #        boxplot_trch,
  #        width = 6.5,
  #        height = 1.3,
  #        units = "in",
  #        dpi = 600)

# eastern cottonwood --------------------------------------#
  # using Jan-Dec
  eaco <- df2 %>%
    filter(common_name == "eastern cottonwood") %>%
    arrange(individual_id, phenogroup, obsdate) %>%
    mutate(prior_no = NA)
  # Calculate prior no
  for (i in 2:nrow(eaco)) {
    eaco$prior_no[i] <- ifelse(
      eaco$individual_id[i] == eaco$individual_id[i-1] & 
        eaco$phenogroup[i] == eaco$phenogroup[i-1],
      as.numeric(eaco$obsdate[i] - eaco$obsdate[i-1]),
      NA)
  }
  # Exclude observations with prior no > 14
  eaco <- eaco %>%
    filter(status == 1) %>%
    filter(!is.na(prior_no) & prior_no < 15)
  # Group by year-period
  eaco_gr <- eaco %>%
    group_by(common_name, individual_id, phenogroup, yr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = min(doy),
              .groups = "keep") %>%
    mutate(group_levels = case_when(
      phenogroup == "lf" ~ "Leaves",
      phenogroup == "fl" ~ "Flowers",
      phenogroup == "flo" ~ "Open flowers",
      phenogroup == "fr" ~ "Fruits",
      phenogroup == "frr" ~ "Ripe fruits",
      phenogroup == "lfe" ~ "Colored leaves"
    )) %>%
    mutate(group_levels = factor(group_levels,
                                 levels = c("Colored leaves",
                                            "Ripe fruits",
                                            "Fruits",
                                            "Open flowers",
                                            "Flowers",
                                            "Leaves"))) %>%
    data.frame()

  # Create nice labels for dates on X axes in plots
  # mindoy <- min(onsets_sub$firstyes)
  # maxdoy <- max(onsets_sub$firstyes)
  doys <- seq(1, 365)
  plotdates <- as.Date(paste(2023, doys, sep = "-"), "%Y-%j")
  date1ind <- which(day(plotdates) == 1)
  doybreaks <- doys[date1ind]
  datebreaks <- plotdates[date1ind] %>% format("%m/%d")
  datelims <- c(min(doys), max(doys) + 10)
  
  boxplot_eaco <- ggplot(data = eaco_gr,
                         aes(x = group_levels, 
                             y = first_yes_doy, 
                             fill = group_levels)) +
    geom_boxplot() +
    stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, 
                 family = "sans", size = 3) +
    scale_y_continuous(limits = datelims, breaks = doybreaks, 
                       labels = datebreaks) +
    scale_fill_manual(values = color_vec) +
    facet_grid(~common_name) +
    labs(x = "", y = "First day observed") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9),
          strip.text = element_text(size = 9))
  
  # ggsave("output/onsets-plot-abq-eaco.png",
  #        boxplot_eaco,
  #        width = 6.5,
  #        height = 2.5,
  #        units = "in",
  #        dpi = 600)

# golden currant -----#
  # using Jan-Dec
  gocu <- df2 %>%
    filter(common_name == "golden currant") %>%
    arrange(individual_id, phenogroup, obsdate) %>%
    mutate(prior_no = NA)
  # Calculate prior no
  for (i in 2:nrow(gocu)) {
    gocu$prior_no[i] <- ifelse(
      gocu$individual_id[i] == gocu$individual_id[i-1] & 
        gocu$phenogroup[i] == gocu$phenogroup[i-1],
      as.numeric(gocu$obsdate[i] - gocu$obsdate[i-1]),
      NA)
  }
  # Exclude observations with prior no > 14
  gocu <- gocu %>%
    filter(status == 1) %>%
    filter(!is.na(prior_no) & prior_no < 15)
  # Group by year-period
  gocu_gr <- gocu %>%
    group_by(common_name, individual_id, phenogroup, yr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = min(doy),
              .groups = "keep") %>%
    mutate(group_levels = case_when(
      phenogroup == "lf" ~ "Leaves",
      phenogroup == "fl" ~ "Flowers",
      phenogroup == "flo" ~ "Open flowers",
      phenogroup == "fr" ~ "Fruits",
      phenogroup == "frr" ~ "Ripe fruits",
      phenogroup == "lfe" ~ "Colored leaves"
    )) %>%
    mutate(group_levels = factor(group_levels,
                                 levels = c("Colored leaves",
                                            "Ripe fruits",
                                            "Fruits",
                                            "Open flowers",
                                            "Flowers",
                                            "Leaves"))) %>%
    data.frame()
  
  # Create nice labels for dates on X axes in plots
  # mindoy <- min(onsets_sub$firstyes)
  # maxdoy <- max(onsets_sub$firstyes)
  doys <- seq(1, 365)
  plotdates <- as.Date(paste(2023, doys, sep = "-"), "%Y-%j")
  date1ind <- which(day(plotdates) == 1)
  doybreaks <- doys[date1ind]
  datebreaks <- plotdates[date1ind] %>% format("%m/%d")
  datelims <- c(min(doys), max(doys) + 10)
  
  boxplot_gocu <- ggplot(data = gocu_gr,
                         aes(x = group_levels, 
                             y = first_yes_doy, 
                             fill = group_levels)) +
    geom_boxplot() +
    stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, 
                 family = "sans", size = 3) +
    scale_y_continuous(limits = datelims, breaks = doybreaks, 
                       labels = datebreaks) +
    scale_fill_manual(values = color_vec) +
    facet_grid(~common_name) +
    labs(x = "", y = "First day observed") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9),
          strip.text = element_text(size = 9))
  
  # ggsave("output/onsets-plot-abq-gocu.png",
  #        boxplot_gocu,
  #        width = 6.5,
  #        height = 2.5,
  #        units = "in",
  #        dpi = 600)

# Siberian elm -----#
  # Plot Oct - next Dec, for colored leaves, start in year 2
  siel <- df2 %>%
    filter(common_name == "Siberian elm") %>%
    arrange(individual_id, phenogroup, obsdate) %>%
    mutate(prior_no = NA)
  # Calculate prior no
  for (i in 2:nrow(siel)) {
    siel$prior_no[i] <- ifelse(
      siel$individual_id[i] == siel$individual_id[i-1] & 
        siel$phenogroup[i] == siel$phenogroup[i-1],
      as.numeric(siel$obsdate[i] - siel$obsdate[i-1]),
      NA)
  }
  # Exclude observations with prior no > 14
  siel <- siel %>%
    filter(status == 1) %>%
    filter(!is.na(prior_no) & prior_no < 15)
  # Group by year-period (calendar for non-flowers)
  siel_gr1 <- siel %>%
    filter(!phenogroup %in% c("fl", 'flo')) %>%
    group_by(common_name, individual_id, phenogroup, yr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = min(doy),
              .groups = "keep") %>%
    data.frame()
  # Group by year-period (water for flower)
  siel_gr2 <- siel %>%
    filter(phenogroup %in% c("fl", 'flo')) %>%
    group_by(common_name, individual_id, phenogroup, wateryr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = doy[dowy == min(dowy)],
              first_yes_dowy = min(dowy), 
              .groups = "keep") %>%
    data.frame()
  # Adding day of year (with negatives, so we can plot this with other 
  # phenophases that are calculated on a calendar year)
  siel_gr2 <- siel_gr2 %>%
    mutate(first_yes_doy_plot = ifelse(month(first_yes) > 9, 
                                       first_yes_doy - 365, 
                                       first_yes_doy))
  # So x-axis will range from Oct in year one (~ negative 62) to 365
  siel_gr2 %>% arrange(phenogroup)
  siel_gr <- siel_gr2 %>%
    select(-c(first_yes_doy, first_yes_dowy)) %>%
    rename(yr = wateryr,
           first_yes_doy = first_yes_doy_plot) %>%
    rbind(siel_gr1) %>%
    mutate(group_levels = case_when(
      phenogroup == "lf" ~ "Leaves",
      phenogroup == "fl" ~ "Flowers",
      phenogroup == "flo" ~ "Open flowers",
      phenogroup == "fr" ~ "Fruits",
      phenogroup == "frr" ~ "Ripe fruits",
      phenogroup == "lfe" ~ "Colored leaves"
    )) %>%
    mutate(group_levels = factor(group_levels,
                                 levels = c("Colored leaves",
                                            "Ripe fruits",
                                            "Fruits",
                                            "Open flowers",
                                            "Flowers",
                                            "Leaves"))) %>%
    data.frame()
  
  # Create nice labels for dates on X axes in plots
  doys_25 <- seq(yday(as.Date("2023-10-01")),365)
  doys_26 <- seq(1, 365)
  plotdates_25 <- as.Date(paste(2025, doys_25, sep = "-"), "%Y-%j")
  plotdates_26 <- as.Date(paste(2026, doys_26, sep = "-"), "%Y-%j")
  plotdates <- c(plotdates_25, plotdates_26)
  date1ind <- which(day(plotdates) == 1)[c(1,3,5,7,9,11,13,15)]
  doybreaks <- doys[date1ind]
  datebreaks <- plotdates[date1ind] %>% format("%m/%d")
  datelims <- c(yday(as.Date("2023-10-01")) - 365, 365 + 10)
  
  boxplot_siel <- ggplot(data = siel_gr,
                         aes(x = group_levels, 
                             y = first_yes_doy, 
                             fill = group_levels)) +
    geom_boxplot() +
    geom_hline(yintercept = 1, linetype = "dotted") +
    stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, 
                 family = "sans", size = 3) +
    scale_y_continuous(limits = datelims, breaks = doybreaks, 
                       labels = datebreaks) +
    scale_fill_manual(values = color_vec) +
    facet_grid(~common_name) +
    labs(x = "", y = "First day observed") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9),
          strip.text = element_text(size = 9))
  
  # ggsave("output/onsets-plot-abq-siel.png",
  #        boxplot_siel,
  #        width = 6.5,
  #        height = 2.5,
  #        units = "in",
  #        dpi = 600)

# rubber rabitbrush -----#
  # using Jan-Dec for fruits
  # using Jul-Jun for flowers
  # not modeling leaves (no signal)
  rura <- df2 %>%
    filter(common_name == "rubber rabbitbrush") %>%
    filter(!phenogroup %in% c("lf", "lfe")) %>%
    arrange(individual_id, phenogroup, obsdate) %>%
    mutate(prior_no = NA)
  # Calculate prior no
  for (i in 2:nrow(rura)) {
    rura$prior_no[i] <- ifelse(
      rura$individual_id[i] == rura$individual_id[i-1] & 
        rura$phenogroup[i] == rura$phenogroup[i-1],
      as.numeric(rura$obsdate[i] - rura$obsdate[i-1]),
      NA)
  }
  # Exclude observations with prior no > 14
  rura <- rura %>%
    filter(status == 1) %>%
    filter(!is.na(prior_no) & prior_no < 15)
  # Group by year-period for flowers
  rura_gr1 <- rura %>%
    filter(phenogroup %in% c("fl", "flo")) %>%
    group_by(common_name, individual_id, phenogroup, yr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = min(doy),
              .groups = "keep") %>%
    data.frame()
  # Group by year-period for fruits
  rura_gr2 <- rura %>%
    filter(phenogroup %in% c("fr", "frr")) %>%
    group_by(common_name, individual_id, phenogroup, summeryr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = doy[dosy == min(dosy)],
              first_yes_dosy = min(dosy), 
              .groups = "keep") %>%
    data.frame()
  rura_gr2 %>% arrange(phenogroup)
  # Can use doy for all, but fruit limited to Jul-Dec (first_yes_doy > 180)
  rura_gr <- rura_gr2 %>%
    filter(first_yes_doy > 180) %>%
    select(-first_yes_dosy) %>%
    rename(yr = summeryr) %>%
    rbind(rura_gr1) %>%
    mutate(group_levels = case_when(
      phenogroup == "lf" ~ "Leaves",
      phenogroup == "fl" ~ "Flowers",
      phenogroup == "flo" ~ "Open flowers",
      phenogroup == "fr" ~ "Fruits",
      phenogroup == "frr" ~ "Ripe fruits",
      phenogroup == "lfe" ~ "Colored leaves"
    )) %>%
    mutate(group_levels = factor(group_levels,
                                 levels = c("Colored leaves",
                                            "Ripe fruits",
                                            "Fruits",
                                            "Open flowers",
                                            "Flowers",
                                            "Leaves"))) %>%
    data.frame()
  
  # Create nice labels for dates on X axes in plots
  # mindoy <- min(onsets_sub$firstyes)
  # maxdoy <- max(onsets_sub$firstyes)
  doys <- seq(1, 365)
  plotdates <- as.Date(paste(2023, doys, sep = "-"), "%Y-%j")
  date1ind <- which(day(plotdates) == 1)
  doybreaks <- doys[date1ind]
  datebreaks <- plotdates[date1ind] %>% format("%m/%d")
  datelims <- c(min(doys), max(doys) + 10)
  
  boxplot_rura <- ggplot(data = rura_gr,
                         aes(x = group_levels, 
                             y = first_yes_doy, 
                             fill = group_levels)) +
    geom_boxplot() +
    stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, 
                 family = "sans", size = 3) +
    scale_y_continuous(limits = datelims, breaks = doybreaks, 
                       labels = datebreaks) +
    scale_fill_manual(values = color_vec) +
    facet_grid(~common_name) +
    labs(x = "", y = "First day observed") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9),
          strip.text = element_text(size = 9))
  
  # ggsave("output/onsets-plot-abq-rura.png",
  #        boxplot_rura,
  #        width = 6.5,
  #        height = 1.9,
  #        units = "in",
  #        dpi = 600)

# screwbean mesquite -----#
  # using Jan-Dec for most
  # using Jul-Jun for ripe fruit
  # not modeling fruit (no signal)
  scme <- df2 %>%
    filter(common_name == "screwbean mesquite") %>%
    filter(phenogroup != "fr") %>%
    arrange(individual_id, phenogroup, obsdate) %>%
    mutate(prior_no = NA)
  # Calculate prior no
  for (i in 2:nrow(scme)) {
    scme$prior_no[i] <- ifelse(
      scme$individual_id[i] == scme$individual_id[i-1] & 
        scme$phenogroup[i] == scme$phenogroup[i-1],
      as.numeric(scme$obsdate[i] - scme$obsdate[i-1]),
      NA)
  }
  # Exclude observations with prior no > 14
  scme <- scme %>%
    filter(status == 1) %>%
    filter(!is.na(prior_no) & prior_no < 15)
  # Group by year-period for everything but ripe fruit
  scme_gr1 <- scme %>%
    filter(phenogroup != "frr") %>%
    group_by(common_name, individual_id, phenogroup, yr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = min(doy),
              .groups = "keep") %>%
    data.frame()
  # Group by year-period for ripe fruit
  scme_gr2 <- scme %>%
    filter(phenogroup  == "frr") %>%
    group_by(common_name, individual_id, phenogroup, summeryr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = doy[dosy == min(dosy)],
              first_yes_dosy = min(dosy), 
              .groups = "keep") %>%
    data.frame()
  # Can use doy for all
  scme_gr <- scme_gr2 %>%
    select(-first_yes_dosy) %>%
    rename(yr = summeryr) %>%
    rbind(scme_gr1) %>%
    mutate(group_levels = case_when(
      phenogroup == "lf" ~ "Leaves",
      phenogroup == "fl" ~ "Flowers",
      phenogroup == "flo" ~ "Open flowers",
      phenogroup == "fr" ~ "Fruits",
      phenogroup == "frr" ~ "Ripe fruits",
      phenogroup == "lfe" ~ "Colored leaves"
    )) %>%
    mutate(group_levels = factor(group_levels,
                                 levels = c("Colored leaves",
                                            "Ripe fruits",
                                            "Fruits",
                                            "Open flowers",
                                            "Flowers",
                                            "Leaves"))) %>%
    data.frame()
  
  # Create nice labels for dates on X axes in plots
  # mindoy <- min(onsets_sub$firstyes)
  # maxdoy <- max(onsets_sub$firstyes)
  doys <- seq(1, 365)
  plotdates <- as.Date(paste(2023, doys, sep = "-"), "%Y-%j")
  date1ind <- which(day(plotdates) == 1)
  doybreaks <- doys[date1ind]
  datebreaks <- plotdates[date1ind] %>% format("%m/%d")
  datelims <- c(min(doys), max(doys) + 10)
  
  boxplot_scme <- ggplot(data = scme_gr,
                         aes(x = group_levels, 
                             y = first_yes_doy, 
                             fill = group_levels)) +
    geom_boxplot() +
    stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, 
                 family = "sans", size = 3) +
    scale_y_continuous(limits = datelims, breaks = doybreaks, 
                       labels = datebreaks) +
    scale_fill_manual(values = color_vec) +
    facet_grid(~common_name) +
    labs(x = "", y = "First day observed") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9),
          strip.text = element_text(size = 9))
  
  # ggsave("output/onsets-plot-abq-scme.png",
  #        boxplot_scme,
  #        width = 6.5,
  #        height = 2.2,
  #        units = "in",
  #        dpi = 600)

# stretchberry -----#
  # using Jan-Dec, but not fruit (no signal)
  stbe <- df2 %>%
    filter(common_name == "stretchberry") %>%
    filter(!phenogroup %in% c("fr", "frr")) %>%
    arrange(individual_id, phenogroup, obsdate) %>%
    mutate(prior_no = NA)
  # Calculate prior no
  for (i in 2:nrow(stbe)) {
    stbe$prior_no[i] <- ifelse(
      stbe$individual_id[i] == stbe$individual_id[i-1] & 
        stbe$phenogroup[i] == stbe$phenogroup[i-1],
      as.numeric(stbe$obsdate[i] - stbe$obsdate[i-1]),
      NA)
  }
  # Exclude observations with prior no > 14
  stbe <- stbe %>%
    filter(status == 1) %>%
    filter(!is.na(prior_no) & prior_no < 15)
  # Group by year-period
  stbe_gr <- stbe %>%
    group_by(common_name, individual_id, phenogroup, yr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = min(doy),
              .groups = "keep") %>%
    mutate(group_levels = case_when(
      phenogroup == "lf" ~ "Leaves",
      phenogroup == "fl" ~ "Flowers",
      phenogroup == "flo" ~ "Open flowers",
      phenogroup == "fr" ~ "Fruits",
      phenogroup == "frr" ~ "Ripe fruits",
      phenogroup == "lfe" ~ "Colored leaves"
    )) %>%
    mutate(group_levels = factor(group_levels,
                                 levels = c("Colored leaves",
                                            "Ripe fruits",
                                            "Fruits",
                                            "Open flowers",
                                            "Flowers",
                                            "Leaves"))) %>%
    data.frame()
  
  # Create nice labels for dates on X axes in plots
  doys <- seq(1, 365)
  plotdates <- as.Date(paste(2023, doys, sep = "-"), "%Y-%j")
  date1ind <- which(day(plotdates) == 1)
  doybreaks <- doys[date1ind]
  datebreaks <- plotdates[date1ind] %>% format("%m/%d")
  datelims <- c(min(doys), max(doys) + 10)
  
  boxplot_stbe <- ggplot(data = stbe_gr,
                         aes(x = group_levels, 
                             y = first_yes_doy, 
                             fill = group_levels)) +
    geom_boxplot() +
    stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, 
                 family = "sans", size = 3) +
    scale_y_continuous(limits = datelims, breaks = doybreaks, 
                       labels = datebreaks) +
    scale_fill_manual(values = color_vec) +
    facet_grid(~common_name) +
    labs(x = "", y = "First day observed") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9),
          strip.text = element_text(size = 9))
  
  # ggsave("output/onsets-plot-abq-stbe.png",
  #        boxplot_stbe,
  #        width = 6.5,
  #        height = 1.9,
  #        units = "in",
  #        dpi = 600)

# fourwing saltbush -----#
  # using Jan-Dec for flowers
  # using Jul-Jun for ripe fruit
  # not modeling leaves or fruit (no signal)
  fwsa <- df2 %>%
    filter(common_name == "fourwing saltbush") %>%
    filter(phenogroup %in% c("fl", "flo", "frr")) %>%
    arrange(individual_id, phenogroup, obsdate) %>%
    mutate(prior_no = NA)
  # Calculate prior no
  for (i in 2:nrow(fwsa)) {
    fwsa$prior_no[i] <- ifelse(
      fwsa$individual_id[i] == fwsa$individual_id[i-1] & 
        fwsa$phenogroup[i] == fwsa$phenogroup[i-1],
      as.numeric(fwsa$obsdate[i] - fwsa$obsdate[i-1]),
      NA)
  }
  # Exclude observations with prior no > 14
  fwsa <- fwsa %>%
    filter(status == 1) %>%
    filter(!is.na(prior_no) & prior_no < 15)
  # Group by year-period for everything but ripe fruit
  fwsa_gr1 <- fwsa %>%
    filter(phenogroup != "frr") %>%
    group_by(common_name, individual_id, phenogroup, yr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = min(doy),
              .groups = "keep") %>%
    data.frame()
  # Group by year-period for ripe fruit
  fwsa_gr2 <- fwsa %>%
    filter(phenogroup  == "frr") %>%
    group_by(common_name, individual_id, phenogroup, summeryr) %>%
    summarize(first_yes = min(obsdate),
              first_yes_doy = doy[dosy == min(dosy)],
              first_yes_dosy = min(dosy), 
              .groups = "keep") %>%
    data.frame()
  # Can use doy for all, but ripe fruit limited to Jul-Dec (first_yes_doy > 180)
  fwsa_gr <- fwsa_gr2 %>%
    filter(first_yes_doy > 180) %>%
    select(-first_yes_dosy) %>%
    rename(yr = summeryr) %>%
    rbind(fwsa_gr1) %>%
    mutate(group_levels = case_when(
      phenogroup == "lf" ~ "Leaves",
      phenogroup == "fl" ~ "Flowers",
      phenogroup == "flo" ~ "Open flowers",
      phenogroup == "fr" ~ "Fruits",
      phenogroup == "frr" ~ "Ripe fruits",
      phenogroup == "lfe" ~ "Colored leaves"
    )) %>%
    mutate(group_levels = factor(group_levels,
                                 levels = c("Colored leaves",
                                            "Ripe fruits",
                                            "Fruits",
                                            "Open flowers",
                                            "Flowers",
                                            "Leaves"))) %>%
    data.frame()
  
  # Create nice labels for dates on X axes in plots
  doys <- seq(1, 365)
  plotdates <- as.Date(paste(2023, doys, sep = "-"), "%Y-%j")
  date1ind <- which(day(plotdates) == 1)
  doybreaks <- doys[date1ind]
  datebreaks <- plotdates[date1ind] %>% format("%m/%d")
  datelims <- c(min(doys), max(doys) + 10)
  
  boxplot_fwsa <- ggplot(data = fwsa_gr,
                         aes(x = group_levels, 
                             y = first_yes_doy, 
                             fill = group_levels)) +
    geom_boxplot() +
    stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, 
                 family = "sans", size = 3) +
    scale_y_continuous(limits = datelims, breaks = doybreaks, 
                       labels = datebreaks) +
    scale_fill_manual(values = color_vec) +
    facet_grid(~common_name) +
    labs(x = "", y = "First day observed") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9),
          strip.text = element_text(size = 9))
  
  # ggsave("output/onsets-plot-abq-fwsa.png",
  #        boxplot_fwsa,
  #        width = 6.5,
  #        height = 1.6,
  #        units = "in",
  #        dpi = 600)

# horsetail milkweed -----#
# Not enough data on most phenophases

# showy milkweed -----#
# Not enough data (and a lot of noise) for most phenophases

# Can we look at TRENDS in phenophase onsets and/or climate drivers? ----------#
# For now, will just do this for species where 5-6 plants were monitored
# across 10-12 years; phenophases that had distinct onsets and periods when
# the plants were out of phase.
  # eastern cottonwood (all phps)
  # fourwing saltbush (flowers, open flowoers, ripe fruits [2nd half of year])
  # stretchberry (leaves, flowers, open flowers, colored leaves)

php_match <- eaco_gr %>%
  distinct(phenogroup, group_levels) %>%
  mutate(group_levels = factor(group_levels,
                               levels = rev(levels(group_levels))))
  
# Eastern cottonwood
  # Going to restrict first yeses for all php but colored leaves to the 1st half
  # of year. Restrict first yeses for colored leaves to second half of year.
  trenddf <- eaco_gr %>%
    filter((phenogroup != "lfe" & first_yes_doy < 183) |
             (phenogroup == "lfe" & first_yes_doy > 182)) %>%
    mutate(yr0 = yr - min(yr)) %>%
    mutate(season = ifelse(phenogroup == "lfe", "Late season", "Early season"))
  
  # Model with separate trend for each phenophase, random slopes by individual
  m <- lmer(first_yes_doy ~ yr0 * phenogroup + (1|individual_id), data = trenddf)
  summary(m)
  # Get trends for each phenophase
  slopes <- emtrends(m, ~ phenogroup, var = "yr0") %>% 
    data.frame() %>%
    mutate(slope = sprintf("%.2f", round(yr0.trend, 2)),
           lower = sprintf("%.2f", round(lower.CL, 2)),
           upper = sprintf("%.2f", round(upper.CL, 2)))
  ann_text <- slopes %>%
    select(phenogroup, slope, upper, lower) %>%
    left_join(php_match, by = "phenogroup") %>%
    mutate(label = paste0("slope = ", slope, " (", lower, ", ", upper, ")"))
  
  preds <- data.frame(predict_response(m, terms = c("yr0", "phenogroup"))) %>%
    mutate(common_name = "eastern cottonwood") %>%
    rename(yr0 = x, phenogroup = group) %>%
    mutate(yr = yr0 + min(trenddf$yr)) %>%
    left_join(php_match, by = "phenogroup") %>%
    mutate(season = ifelse(phenogroup == "lfe", "Late season", "Early season"))
  
  # Remove years with no observations
  phpdat <- trenddf %>%
    group_by(phenogroup) %>%
    summarize(minyr = min(yr),
              maxyr = max(yr)) %>%
    data.frame() 
  for (php in phpdat$phenogroup) {
    preds$predicted[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.low[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.high[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA
    preds$std.error[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA 
    preds$predicted[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.low[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.high[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
    preds$std.error[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
  }
  
  text_size <- 8
  ci_alpha <- 0.17
  
  trends_eaco <- ggplot(data = preds, aes(x = yr, y = predicted)) +
    geom_line(aes(col = group_levels)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group_levels), 
                alpha = ci_alpha) +
    geom_point(data = trenddf,
               aes(x = yr, y = first_yes_doy, color = factor(individual_id)),
               position = position_dodge(width = 0.3),
               size = 0.8, alpha = 0.5) +
    facet_wrap(~ group_levels, scales = "free_y") +
    scale_fill_manual(values = color_vec) +
    scale_color_manual(values = color_vec) +
    geom_abs_text(data = ann_text,
                  aes(xpos = 0.95, ypos = 0.97, label = label, hjust = 1, vjust = 1,
                  family = "sans", x = NULL, y = NULL), size = text_size/.pt) +
    labs(x = "Year",
         y = "First observation day of year", 
         color = "Phenophase",
         fill = "Phenophase") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = "none",
          axis.text.x = element_text(size = text_size), 
          axis.text.y = element_text(size = text_size),
          strip.text = element_text(size = text_size),
          axis.title.y = element_text(size = text_size + 1))
    # ggsave("output/trends-plot-abq-eaco.png",
    #        trends_eaco,
    #        width = 6.5,
    #        height = 6.5,
    #        units = "in",
    #        dpi = 600)
  
# Fourwing saltbush
  # For ripe fruits, we only considered first yeses that occurred in the second 
  # half of the year (already did this filtering)
  trenddf <- fwsa_gr %>%
    mutate(yr = ifelse(nchar(yr) > 4, str_sub(yr, 1, 4), yr)) %>%
    mutate(yr = as.numeric(yr)) %>%
    mutate(yr0 = yr - min(yr))
  
  # Model with separate trend for each phenophase, random slopes by individual
  m <- lmer(first_yes_doy ~ yr0 * phenogroup + (1|individual_id), data = trenddf)
  summary(m)
  # Get trends for each phenophase
  slopes <- emtrends(m, ~ phenogroup, var = "yr0") %>% 
    data.frame() %>%
    mutate(slope = sprintf("%.2f", round(yr0.trend, 2)),
           lower = sprintf("%.2f", round(lower.CL, 2)),
           upper = sprintf("%.2f", round(upper.CL, 2)))
  ann_text <- slopes %>%
    select(phenogroup, slope, upper, lower) %>%
    left_join(php_match, by = "phenogroup") %>%
    mutate(label = paste0("slope = ", slope, " (", lower, ", ", upper, ")"))
  
  preds <- data.frame(predict_response(m, terms = c("yr0", "phenogroup"))) %>%
    mutate(common_name = "eastern cottonwood") %>%
    rename(yr0 = x, phenogroup = group) %>%
    mutate(yr = yr0 + min(trenddf$yr)) %>%
    left_join(php_match, by = "phenogroup") %>%
    mutate(season = ifelse(phenogroup == "lfe", "Late season", "Early season"))
  
  # Remove years with no observations
  phpdat <- trenddf %>%
    group_by(phenogroup) %>%
    summarize(minyr = min(yr),
              maxyr = max(yr)) %>%
    data.frame() 
  for (php in phpdat$phenogroup) {
    preds$predicted[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.low[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.high[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA
    preds$std.error[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA 
    preds$predicted[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.low[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.high[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
    preds$std.error[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
  }
  
  text_size <- 8
  ci_alpha <- 0.17
  
  trends_fwsa <- ggplot(data = preds, aes(x = yr, y = predicted)) +
    geom_line(aes(col = group_levels)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group_levels), 
                alpha = ci_alpha) +
    geom_point(data = trenddf,
               aes(x = yr, y = first_yes_doy, color = factor(individual_id)),
               position = position_dodge(width = 0.3),
               size = 0.8, alpha = 0.5) +
    facet_wrap(~ group_levels, scales = "free_y") +
    scale_fill_manual(values = color_vec) +
    scale_color_manual(values = color_vec) +
    geom_abs_text(data = ann_text,
                  aes(xpos = 0.95, ypos = 0.97, label = label, hjust = 1, vjust = 1,
                      family = "sans", x = NULL, y = NULL), size = text_size/.pt) +
    labs(x = "Year",
         y = "First observation day of year", 
         color = "Phenophase",
         fill = "Phenophase") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = "none",
          axis.text.x = element_text(size = text_size), 
          axis.text.y = element_text(size = text_size),
          strip.text = element_text(size = text_size),
          axis.title.y = element_text(size = text_size + 1))
  # ggsave("output/trends-plot-abq-fwsa.png",
  #        trends_fwsa,
  #        width = 6.5,
  #        height = 3.25,
  #        units = "in",
  #        dpi = 600)
  
# Stretchberry
  # Going to restrict first yeses for leaves to the 1st half of year. 
  # Restrict first yeses for colored leaves to second half of year.
  trenddf <- stbe_gr %>%
    filter(!(phenogroup == "lf" & first_yes_doy > 182)) %>%
    filter(!(phenogroup == "lfe" & first_yes_doy < 183)) %>%
    mutate(yr0 = yr - min(yr))

  # Model with separate trend for each phenophase, random slopes by individual
  m <- lmer(first_yes_doy ~ yr0 * phenogroup + (1|individual_id), data = trenddf)
  summary(m)
  # Get trends for each phenophase
  slopes <- emtrends(m, ~ phenogroup, var = "yr0") %>% 
    data.frame() %>%
    mutate(slope = sprintf("%.2f", round(yr0.trend, 2)),
           lower = sprintf("%.2f", round(lower.CL, 2)),
           upper = sprintf("%.2f", round(upper.CL, 2)))
  ann_text <- slopes %>%
    select(phenogroup, slope, upper, lower) %>%
    left_join(php_match, by = "phenogroup") %>%
    mutate(label = paste0("slope = ", slope, " (", lower, ", ", upper, ")"))
  
  preds <- data.frame(predict_response(m, terms = c("yr0", "phenogroup"))) %>%
    mutate(common_name = "eastern cottonwood") %>%
    rename(yr0 = x, phenogroup = group) %>%
    mutate(yr = yr0 + min(trenddf$yr)) %>%
    left_join(php_match, by = "phenogroup") %>%
    mutate(season = ifelse(phenogroup == "lfe", "Late season", "Early season"))
  
  # Remove years with no observations
  phpdat <- trenddf %>%
    group_by(phenogroup) %>%
    summarize(minyr = min(yr),
              maxyr = max(yr)) %>%
    data.frame() 
  for (php in phpdat$phenogroup) {
    preds$predicted[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.low[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.high[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA
    preds$std.error[preds$phenogroup == php & preds$yr < phpdat$minyr[phpdat$phenogroup == php]] <- NA 
    preds$predicted[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.low[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
    preds$conf.high[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
    preds$std.error[preds$phenogroup == php & preds$yr > phpdat$maxyr[phpdat$phenogroup == php]] <- NA 
  }
  
  text_size <- 8
  ci_alpha <- 0.17
  
  trends_stbe <- ggplot(data = preds, aes(x = yr, y = predicted)) +
    geom_line(aes(col = group_levels)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group_levels), 
                alpha = ci_alpha) +
    geom_point(data = trenddf,
               aes(x = yr, y = first_yes_doy, color = factor(individual_id)),
               position = position_dodge(width = 0.3),
               size = 0.8, alpha = 0.5) +
    facet_wrap(~ group_levels, scales = "free_y") +
    scale_fill_manual(values = color_vec) +
    scale_color_manual(values = color_vec) +
    scale_x_continuous(breaks = seq(2016, 2024, by = 4)) +
    geom_abs_text(data = ann_text,
                  aes(xpos = 0.95, ypos = 0.97, label = label, hjust = 1, vjust = 1,
                      family = "sans", x = NULL, y = NULL), size = text_size/.pt) +
    labs(x = "Year",
         y = "First observation day of year", 
         color = "Phenophase",
         fill = "Phenophase") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = "none",
          axis.text.x = element_text(size = text_size), 
          axis.text.y = element_text(size = text_size),
          strip.text = element_text(size = text_size),
          axis.title.y = element_text(size = text_size + 1))
  # ggsave("output/trends-plot-abq-stbe.png",
  #        trends_stbe,
  #        width = 6.5,
  #        height = 6.5,
  #        units = "in",
  #        dpi = 600)
  
# Climate effects??
  
  