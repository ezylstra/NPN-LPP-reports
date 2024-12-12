# Example LPP report
# Erin Zylstra
# erinz@usanpn.org
# 2024-11-21

require(rnpn)
require(tidyr)
require(dplyr)
require(lubridate)
require(stringr)
require(ggplot2)
require(terra)
require(lme4)
# require(lmerTest) # Not sure I want this...
# require(geosphere) # for optional site clustering
# require(fpc) # for optional site clustering

rm(list = ls())

#- Specify data of interest ---------------------------------------------------#
# Name and nickname of LPP
lpp <- "Appalachian Mountain Club"
lpp_short <- "amc"
# lpp <- "Earth"
# lpp_short <- 

# Specify years of interest
all_yrs <- 2007:2023

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
    filter(str_detect(network_name, lpp))
  if (nrow(network_id) == 1) {
    message("Good news: One network_id found for LPP")
    lpp_id <- network_id$network_id
  } else if (nrow(network_id) == 0) {
    stop("No matches for LPP name", call. = FALSE)
  } else {
    warning("Multiple network_ids found for LPP", call. = FALSE)
    lpp_id <- c(network_id$network_id)
  }
  
  # Download status/intensity data (one row for each observation of a plant or 
  # animal species and phenophase)
  si <- npn_download_status_data(
    request_source = requestor,
    network_ids = lpp_id,
    years = all_yrs,
    climate_data = FALSE,
    additional_fields = c("site_name", "observedby_person_id")
  )
  #TODO: figure out whether it's ever worth including climate data
  #TODO: figure out where we get info about search methods, time for animals
  # Was going to download observation_time too, but it looks like it wasn't
  # reported frequently.
  
  # Remove some fields that we don't need (to make file smaller)
  si <- si %>%
    select(-c(update_datetime, species))
  
  # Write to file
  write.csv(si, si_csv_name, row.names = FALSE)
} 

# Load csv
si <- read.csv(si_csv_name, na.strings = c(NA, "-9999"))

# Format date and observer
si <- si %>%
  mutate(obsdate = ymd(observation_date)) %>%
  mutate(yr = year(observation_date)) %>%
  select(-observation_date) %>%
  rename(person_id = observedby_person_id)

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

#- Is it necessary to do any checks of site locations? ------------------------#
# Given that LPPs probably have a limited number of sites that were set up by 
# a LPL, going to skip this step for now

#- Extract site, species information ------------------------------------------#
# Removing to keep working objects smaller. Can merge site, species information 
# back in at any point (Note: for animals, individual_id is unique for each 
# combination of species and site)

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
            species_id, genus))

#- Clean up phenophase status data --------------------------------------------#
# Any phenophase_status == -1 or 0 when intensity/abundance value present?
count(si, phenophase_status, is.na(intensity_value), is.na(abundance_value))
# Apparently, this does happen (at least for -1)
#TODO: QA/QC fix. Can this be prevented when entering data?

# Change unknown phenophase_status to 1 if an intensity/abundance value reported
# and then remove any remaining records with phenophase status unknown
si <- si %>%
  mutate(phenophase_status = case_when(
    !is.na(intensity_value) ~ 1,
    !is.na(abundance_value) ~ 1,
    .default = phenophase_status
  )) %>%
  filter(phenophase_status %in% c(0, 1))

#- Deal with multiple observations by the same person -------------------------#
# Same plant or animal (same site), same phenophase, same date - same observer

# For now, will keep record with more advanced phenophase or higher 
# intensity/abundance value. Will do this by sorting observations in descending
# order and keeping only the first
inddateobsp <- si %>%
  group_by(common_name, individual_id, obsdate, person_id, phenophase_id) %>%
  summarize(n_obs = n(),
            .groups = "keep") %>%
  data.frame()
inddateobsp$obsnum <- 1:nrow(inddateobsp)

si <- si %>%
  arrange(person_id, individual_id, obsdate, phenophase_id, 
          desc(phenophase_status), desc(intensity), desc(abundance_value)) %>%
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

# For animals, multiple observations of the same species at a site on the same
# date provides a confounded measure of inter-observer variation and variation 
# in abundance or activity of the species throughout the day.

# To summarize sampling effort (e.g., when and where observations were made), 
# probably best to keep all the observations in the dataset.

# To analyze phenophase onsets or animal activity/occurrence, probably best
# to subset or aggregate observations so there is a maximum of one record for 
# each phenophase associated with a plant/animal at a given site each day.

# Function to calculate sum, but if all values in vector = NA, then return NA
na_max <- function(x) {
  y <- ifelse(sum(is.na(x)) == length(x), NA, max(x, na.rm = TRUE))
  return(y)
}

#- Creating datasets for animal observations ----------------------------------#

# First a dataset that includes all animal observations (can have multiple 
# observations of a species at a site in one day [by different observers]). For
# now will summarize by pheno group (live obs, dead obs, trapped/baited), but 
# will likely revisit this.
animal_obs <- si %>%
  filter(kingdom == "Animalia") %>%
  group_by(person_id, site_id, common_name, individual_id, obsdate,
           day_of_year, yr, pheno_group) %>%
  summarize(status = na_max(phenophase_status),
            abundance = na_max(abundance_value),
            intensity = na_max(intensity),
            .groups = "keep") %>%
  arrange(common_name, site_id, obsdate, person_id) %>%
  data.frame()

# Put into wide form, so each row is a summary of an observation 
# ie, each row is unique combination of species, site, date, and observer
animal_obs <- animal_obs %>%
  mutate(pheno_group2 = case_when(
    pheno_group == "Live observation" ~ "live",
    pheno_group == "Dead observation" ~ "dead",
    .default = "trap"
  )) %>%
  select(-pheno_group) %>%
  pivot_wider(names_from = pheno_group2,
              values_from =c(status, abundance, intensity),
              names_sort = TRUE) %>%
  data.frame()

# Then create a dataset that has no more than one observation for each species, 
# site, date, and phenophase. (Not simplifying to phenogroup yet since I'm not
# sure what we're doing with this yet). For an animal species, multiple 
# observations aren't nearly as problematic as they might be for plants because 
# a status = 0 could easily be a detection issue. If we assume there are no 
# false positives, then we can simply use the maximum status value and maximum 
# intensity or abundance value across observations
inddatep <- si %>%
  group_by(common_name, individual_id, obsdate, phenophase_id) %>%
  summarize(n_obs = n(),
            .groups = "keep") %>%
  data.frame()
inddatep$obsnum <- 1:nrow(inddatep)

animal_obs2 <- si %>%
  filter(kingdom == "Animalia") %>%
  arrange(common_name, individual_id, obsdate, phenophase_id, 
          desc(phenophase_status), desc(abundance_value), desc(intensity)) %>%
  left_join(select(inddatep, -c(n_obs, common_name)),
            by = c("individual_id", "obsdate", "phenophase_id")) %>%
  mutate(dups = sequence(rle(as.character(obsnum))$lengths)) %>%
  filter(dups == 1) %>%
  select(-c(obsnum, dups, person_id)) %>%
  arrange(common_name, individual_id, obsdate, phenophase_id)

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

  # A little data exploration...
  # For each species, there is only 1 phenophase per pheno class.
    # si_sub_plants %>%
    #   group_by(common_name, class_id) %>%
    #   summarize(n_phenophases = length(unique(phenophase_id)),
    #             .groups = "keep") %>%
    #   data.frame()
  # For each species, may have more than 1-3 phenophases per pheno group
    # si_sub_plants %>%
    #   group_by(common_name, pheno_group) %>%
    #   summarize(n_phenophases = length(unique(phenophase_id)),
    #             .groups = "keep") %>%
    #   data.frame()
  
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
                end_flower = sum(class_id == 9),
                fruit = sum(class_id == 10),
                unripe = sum(class_id == 11),
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

  # A little more data exploration... 
  # Is there just one intensity_category_id for each plant phenophase?
    # int_cat_sp <- si_sub_plants %>%
    #   group_by(common_name, phenophase_id, phenophase_description) %>%
    #   summarize(int_cats = length(unique(intensity_category_id)),
    #             .groups = "keep") %>%
    #   data.frame()
    # count(filter(si_sub_plants, common_name == "yellow birch", phenophase_id == "483"),
    #       intensity_category_id) # 40, 73, or NA
  # No, but I think that's because the counts include intensity categories of NA
  # and because the intensity category may have shifted over time:
    # count(filter(si_sub_plants, common_name == "yellow birch", phenophase_id == "483"),
    #       intensity_category_id, yr) # 40 used in 2014-2015; 73 used in 2016-2023
    # test <- si_sub_plants %>%
    #   group_by(individual_id, phenophase_id, yr) %>%
    #   summarize(int_cat = length(unique(intensity_category_id[!is.na(intensity_category_id)])),
    #             .groups = "keep") %>%
    #   data.frame()
    # count(test, int_cat) 
  # In a single year, only one intensity category (if any) for that plant, phenophase

# All this means that we put data in wide form, with 2 columns for each 
# plant phenophase: status, intensity

# Putting data in wide form:
# Can use phenophase class_id instead of phenophase_id because there's a maximum
# of one phenophase per class for each species.

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
    class_id == 9 ~ "end_flower",
    class_id == 10 ~ "fruit",
    class_id == 11 ~ "unripe_fruit",
    class_id == 12 ~ "ripe_fruit",
    class_id == 13 ~ "drop_fruit",
    .default = NA
  )) %>%
  mutate(group_short = case_when(
    pheno_group == "Dead observation" ~ "dead",
    pheno_group == "Live observation" ~ "live",
    pheno_group == "Trapped/baited" ~ "trap",
    pheno_group == "Leaves" ~ "leaf",
    pheno_group == "Leaf senescence" ~ "leaf_end",
    pheno_group == "Flowers" ~ "flower",
    pheno_group == "Open flowers" ~ "flower_open",
    pheno_group == "Flower senescence" ~ "flower_end",
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
    group_short == "flower_end" ~ "fle",
    group_short == "fruit" ~ "fr",
    group_short == "fruit_ripe" ~ "frr"
  )) %>%
  arrange(class_id) %>%
  mutate(class_id2 = str_pad(class_id, width = 2, pad = 0))

# Put in wide form
plant_obs <- si_sub_plants %>%
  left_join(select(pl_pheno_classes, class_id, class_id2), 
            by = "class_id") %>%
  select(-c(phenophase_id, phenophase_description, pheno_group, intensity_value, 
            intensity_category_id, intensity_type, abundance_value, class_id)) %>%
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
  # If unripe fruit = 1, fruit = 1
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
  # count(plant_obs, status_young_leaf, status_leaf)
  # count(plant_obs, status_color_leaf, status_leaf)
  # count(plant_obs, status_open_flower, status_flower)
  # count(plant_obs, status_pollen, status_flower)
  # count(plant_obs, status_unripe_fruit, status_fruit)
  # count(plant_obs, status_ripe_fruit, status_fruit)

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
    (status_unripe_fruit == 1 & status_fruit == 0) ~ 1,
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
  # animal_obs = dataframe where each row has data for a unique combination of 
    # individual id (site * species), date, and observer. Best for 
    # summarizing effort (by pheno group). May contain multiple observations of 
    # a species at a site in a day.
  # animal_obs2 = dataframe where each row has data for a unique combination 
    # of individual id (site * species), date, and phenophase. Best for
    # summarizing animal activity/detections. No more than one observation of
    # a species and phenophase each day.
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
# Find number of observations of animals
site_a <- animal_obs %>%
  group_by(site_id) %>%
  summarize(n_animal_obs = n()) %>%
  data.frame()
# Find number of observations of plants
site_p <- plant_obs %>%
  group_by(site_id) %>%
  summarize(n_plant_obs = n()) %>%
  data.frame()
# Combine and clean up
sitest <- sitest %>%
  left_join(site_a, by = "site_id") %>%
  left_join(site_p, by = "site_id") %>%
  mutate(across(n_animal_obs:n_plant_obs,  ~ replace_na(.x, 0))) %>%
  mutate(obs_total = n_animal_obs + n_plant_obs) %>%
  rename(obs_animal = n_animal_obs,
         obs_plant = n_plant_obs)
#TODO: Decide whether we want filters for number of years (>= 3?), number
# of observations (> 100?), or number of species?

# Table with the number of observations per year?
# Or maybe a figure, with total observations per year, and a few smaller lines 
# for particular sites (or species) of interest?

#- Table: Yearly summaries ----------------------------------------------------#

# No. of sites, No. of observers, no. plant and animal observations by year
yearst <- si %>%
  group_by(yr) %>%
  summarize(sites = length(unique(site_id)),
            observers = length(unique(person_id))) %>% 
  data.frame()
# Find number of observations of animals
year_a <- animal_obs %>%
  group_by(yr) %>%
  summarize(obs_animal = n()) %>%
  data.frame()
# Find number of observations of plants
year_p <- plant_obs %>%
  group_by(yr) %>%
  summarize(obs_plant = n()) %>%
  data.frame()
# Combine and clean up
yearst <- yearst %>%
  left_join(year_a, by = "yr") %>%
  left_join(year_p, by = "yr") %>%
  mutate(across(obs_animal:obs_plant,  ~ replace_na(.x, 0))) %>%
  mutate(obs_total = obs_animal + obs_plant)
#TODO: Decide whether we want filters for number of sites or observations

#- Map(s) ---------------------------------------------------------------------#

# Should be options if there's one site, 1-9 sites, 10 or more sites

# Might also be options we could turn on or off (eg, color gradiation depending
# on the number of years or observations, inset map if a subset are in very
# close proximity)

# For now, see start in site-maps.R

#- Download and process climate data ------------------------------------------#

# Download in download-daymet.R

# Load climate data
daymet_csv <- paste0("climate-data/daymet-", lpp_short,  
                     "-thru", max(all_yrs), ".csv")
daymet_zip <- str_replace(daymet_csv, ".csv", ".zip")
unzip(daymet_zip)
clim_full <- read.csv(daymet_csv)
file.remove(daymet_csv)

# Summarize data for each season (as defined previously by NPN)
  # Spring = Mar-May
  # Summer = Jun-Aug
  # Fall = Sep-Nov
  # Winter = Dec-Feb (assigned to year for Jan-Feb)

# Just summarize for temperature, precipitation
clims <- clim_full %>%
  select(-c(daylength_s, rad_Wm2, vp_Pa, swe_kgm2))

# Assign "seasonyr" (so Dec is associated with Jan, Feb in next calendar year)
clims <- clims %>%
  mutate(date = ymd(date)) %>%
  mutate(seasonyr = if_else(mon == 12, year + 1, year))

# Remove winter seasons where we don't have data from all months and add a 
# season label
clims <- clims %>%
  # Remove data from Jan-Feb in first year with climate data
  filter(!(seasonyr == min(seasonyr) & mon %in% 1:2)) %>%
  # Remove data from Dec in last year with climate data
  filter(seasonyr != max(seasonyr)) %>%
  # Season labels
  mutate(season = case_when(
    mon %in% 3:5 ~ "spring",
    mon %in% 6:8 ~ "summer",
    mon %in% 9:11 ~ "fall",
    .default = "winter"
  )) %>%
  mutate(season = factor(season, 
                         levels = c("winter", "spring", "summer", "fall")))

# Calculate mean daily temperatures
clims <- clims %>%
  mutate(tmean_c = round((tmax_c + tmin_c) / 2, 2))

# Summarize data by season (accumulated precip, mean tmin, tmean, tmax)
clims <- clims %>%
  group_by(site, seasonyr, season) %>%
  summarize(prcp = sum(prcp_mm),
            tmin = mean(tmin_c),
            tmean = mean(tmean_c),
            tmax = mean(tmax_c),
            .groups = "keep") %>%
  data.frame()

#- Explore site locations -----------------------------------------------------#

# Wondering whether any/all sites could be aggregated for summarizing and 
# analyzing the data. For the moment, won't worry about filtering by the number 
# of species or observations.

# Could do this by evaluating the distances between locations
# Could also do this by looking at climate data

# # Start by calculating pairwise distances
# sitesll <- sites %>%
#   rename(lon = longitude, lat = latitude) %>%
#   vect(geom = c("lon", "lat"), crs = "epsg:4326")
# dist.matrix <- geosphere::distm(as.matrix(sites[, c("longitude", "latitude")]))
# 
# # Cluster analysis (using medoids/k-means)?
# # Alternatively, could use cluster or apcluster packages...
# pamk.best <- fpc::pamk(sites[, c("longitude", "latitude")])
# # pamk.best$nc = optimal number of clusters
# # pamk.best$pamobject$clustering = cluster IDs
# Clusters <- factor(pamk.best$pamobject$cluster)
# ggplot() +
#   geom_point(data = sites, aes(x = longitude, y = latitude, 
#                               color = Clusters),
#              alpha = 0.3, size = 2.5)
# 
# # Differences in climate among sites
# clims_c <- clims %>%
#   group_by(site, season) %>%
#   summarize(prcp = mean(prcp),
#             tmin = mean(tmin),
#             tmean = mean(tmean),
#             tmax = mean(tmax),
#             .groups = "keep") %>%
#   pivot_wider(id_cols = site,
#               names_from = season,
#               values_from = c(prcp, tmin, tmean, tmax)) %>%
#   data.frame() %>%
#   mutate(across(prcp_winter:tmax_fall, ~scale(.x)[, 1]))
# pamk.best.c <- fpc::pamk(clims_c[, -1])
# Clusters.c <- factor(pamk.best.c$pamobject$cluster)
# ggplot() +
#   geom_point(data = sites, aes(x = longitude, y = latitude, 
#                                color = Clusters.c),
#              alpha = 0.3, size = 2.5)
# 
# # Pausing this for now.
# rm(clims_c)

#- Aggregate status information within phenophase group -----------------------#

# Aggregate status information within pheno group (doesn't make sense to 
# aggregate intensity values, since they may not always be the same type)

# For now, will use 7 phenophase groups for plants

# Create new data frame without intensity data and other data summaries
plg_status <- plant_obs2 %>% 
  select(-c(contains("intens"), "n_status", "n_yes", "n_observations"))

# Aggregate information across classes within each pheno group
# TODO: make this easier by putting htings into long form first?
for (group in unique(pl_pheno_classes$group_code)) {
  cols <- pl_pheno_classes$stat_cols_c[pl_pheno_classes$group_code == group]
  plg_status[,paste0("sum_", group)] <- apply(as.matrix(plg_status[,cols]), 
                                              1, na_max)
}
# Remove columns associated with classes and rename columns with status data
plg_status <- plg_status %>% select(-contains("status_"))
colnames(plg_status) <- str_replace(colnames(plg_status), "sum", "status")

#- Decide on calendar vs water year -------------------------------------------#

# Identify whether we want to evaluate onset dates within a calendar year
# or within a water year (1 Oct - 30 Sep)
water_year <- FALSE

# if we're using water year, create new yr and doy columns
plg_status <- plg_status %>% rename(doy = day_of_year)
if (water_year) {
  plg_status <- plg_status %>%
    mutate(mon = month(obsdate)) %>%
    mutate(wateryr = if_else(mon %in% 10:12, yr + 1, yr)) %>%
    mutate(wateryrday0 = make_date(year = wateryr - 1, month = 9, day = 30)) %>%
    mutate(doy = as.numeric(obsdate - wateryrday0)) %>%
    select(-c(mon, wateryrday0, yr)) %>%
    rename(yr = wateryr)
}

#- Calculate/extract onset dates ----------------------------------------------#

# For now, just using the first "yes" of the year (or water year).

# Create a unique ID for each plant-year to make various matches easier
plg_status <- plg_status %>%
  mutate(plantyr = paste0(individual_id, "_", yr))

# Create a new dataframe summarizing data available for each plant-year
py <- plg_status %>%
  group_by(plantyr, individual_id, yr, common_name, site_id) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_status_lf = sum(!is.na(status_lf)),
            n_status_lfe = sum(!is.na(status_lfe)),
            n_status_fl = sum(!is.na(status_fl)),
            n_status_flo = sum(!is.na(status_flo)),
            n_status_fle = sum(!is.na(status_fle)),
            n_status_fr = sum(!is.na(status_fr)),
            n_status_frr = sum(!is.na(status_frr)),
            n_lf_yes = sum(!is.na(status_lf) & status_lf == 1),
            n_lfe_yes = sum(!is.na(status_lfe) & status_lfe == 1),
            n_fl_yes = sum(!is.na(status_fl) & status_fl == 1),
            n_flo_yes = sum(!is.na(status_flo) & status_flo == 1),
            n_fle_yes = sum(!is.na(status_fle) & status_fle == 1),
            n_fr_yes = sum(!is.na(status_fr) & status_fr == 1),
            n_frr_yes = sum(!is.na(status_frr) & status_frr == 1),
            .groups = "keep") %>%
  data.frame()

# Remove years when plant was observed only once
py <- py %>%
  filter(n_obs > 1)

# Add columns for onset info (for all pheno groups, even if there was never any
# status information collected?)
# *_firstyes = DOY with first "yes" response
# *-lastno = number of days between prior "no" response and first "yes" response 
# (will be NA if no prior "no")
py$lf_firstyes <- NA
py$lf_lastno <- NA
py$lfe_firstyes <- NA
py$lfe_lastno <- NA
py$fl_firstyes <- NA
py$fl_lastno <- NA
py$flo_firstyes <- NA
py$flo_lastno <- NA
py$fle_firstyes <- NA
py$fle_lastno <- NA
py$fr_firstyes <- NA
py$fr_lastno <- NA
py$frr_firstyes <- NA
py$frr_lastno <- NA

# Loop through each plant-year (note that this can take a few minutes)
grps <- unique(pl_pheno_classes$group_code)
grpsc <- unique(pl_pheno_classes$stat_cols_g)

# Loop through plant-years
for (i in 1:nrow(py)) {
  
  # Extract data for plant-year 
  tmp1 <- plg_status[plg_status$plantyr == py$plantyr[i],]
  
  # Loop through pheno groups
  for (j in 1:length(grps)) {
    if (py[i, paste0("n_", grps[j], "_yes")] == 0) {next}
    
    # Extract only those dates when pheno group status recorded
    tmp2 <- tmp1[!is.na(tmp1[,grpsc[j]]),]
    
    # Identify rows of tmp2 with first "yes" and first "no"
    first_yes_ind <- first(which(tmp2[,grpsc[j]] == 1))
    first_no_ind <- first(which(tmp2[,grpsc[j]] == 0))
    # Extract onset date (first doy with yes)
    py[i, paste0(grps[j], "_firstyes")] <- tmp2$doy[first_yes_ind]
    # Calculate number of days since last no. If there isn't a prior no, then
    # days_lastno = NA
    if (first_yes_ind > first_no_ind & !is.na(first_no_ind)) {
      py[i, paste0(grps[j], "_lastno")] <- 
        py[i, paste0(grps[j], "_firstyes")] - tmp2$doy[first_yes_ind - 1]
    } else {
      py[i, paste0(grps[j], "_lastno")] <- NA
    }
  }
}
# Could look at ranges of first yes doy to see whether we're running up against
# boundaries and should reconsider using calendar year or water year

#- Summarize/visualize onset data ---------------------------------------------#

# Within each year, filter to include observations only if there is a preceding 
# "no" within X days
# Identify what the cutoff (X) will be:
prior_days <- 14

# Easiest to filter the data if we put in long form first
onsets <- py %>%
  select(-c(n_obs, first_obs, last_obs, 
            contains("n_status"), contains("_yes"))) 
onsets <- onsets %>%
  pivot_longer(cols = c(ends_with("firstyes"),ends_with("lastno")),
               names_to = c("phenogroup", ".value"),
               names_sep = "_") %>%
  data.frame()
onsets <- onsets %>%
  filter(!is.na(firstyes)) %>%
  filter(!is.na(lastno) & lastno <= prior_days)
# Check
# sum(onsets$phenogroup == "flo"); sum(py$flo_lastno <= 14, na.rm = TRUE)

# Look at sample sizes
onsets %>% select(common_name, phenogroup) %>% table

# Look at how sample sizes might change if we reduced the window for prior no
onset_ss <- onsets %>%
  group_by(common_name, phenogroup) %>%
  summarize(n_firstyes_14 = n(),
            n_firstyes_7 = sum(lastno <= 7),
            .groups = "keep") %>%
  data.frame() %>%  
  mutate(prop_7 = round(n_firstyes_7 / n_firstyes_14, 2))
onset_ss
# By species
onset_ss %>%
  group_by(common_name) %>%
  summarize(n_obs_14 = sum(n_firstyes_14),
            mean_prop_7 = round(mean(prop_7), 2)) %>%
  data.frame() %>%
  arrange(desc(mean_prop_7))
# By phenophase group
onset_ss %>%
  group_by(phenogroup) %>%
  summarize(n_obs_14 = sum(n_firstyes_14),
            mean_prop_7 = round(mean(prop_7), 2)) %>%
  data.frame() %>%
  arrange(desc(mean_prop_7))

# Create a column in onsets df with nice names of phenophase groups, for plots
onsets <- onsets %>%
  mutate(group_labels = case_when(
    phenogroup == "lf" ~ "Leaves",
    phenogroup == "lfe" ~ "Leaf senescence",
    phenogroup == "fl" ~ "Flowers",
    phenogroup == "flo" ~ "Open flowers",
    phenogroup == "fle" ~ "Flower senescence",
    phenogroup == "fr" ~ "Fruits",
    phenogroup == "frr" ~ "Ripe fruits"
  )) %>%
  mutate(group_labels = factor(group_labels,
                               levels = c("Leaf senescence",
                                          "Ripe fruits",
                                          "Fruits",
                                          "Flower senescence",
                                          "Open flowers",
                                          "Flowers",
                                          "Leaves")))

# Assign colors for each phenophase group
color_vec <- c("#7fbf7b",   # Leaves
               "#e9a3c9",   # Flowers
               "#c51b7d",   # Open flowers
               "#b2beb5",   # Flower senescence
               "#e7d4e8",   # Fruits
               "#af8dc3",   # Ripe fruits
               "#8c510a")   # Leaf senescence
# Name vector so colors are consistent across figures (in case not all levels 
# are present in all figures)
names(color_vec) <- rev(levels(onsets$group_labels))

# Add sample sizes for each species-phenogroup combination to onsets
onsets <- onsets %>%
  group_by(common_name, phenogroup) %>%
  mutate(n_obs_sppgroup = n()) %>%
  ungroup() %>%
  data.frame()

# Set minimum number of observations per species-phenogroup to summarize onset 
# dates
min_obs <- 15

# Plot onset for one phenogroup by species (lumping years, sites together)
ggplot(data = filter(onsets, phenogroup == "flo", n_obs_sppgroup >= min_obs),   
       aes(x = common_name, y = firstyes)) +
  geom_boxplot(varwidth = TRUE) + 
  labs(y = "Onset DOY - Open flowers") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank())
# Note: by using the varwidth argument, width of boxplots are proportional
# to the square root of the number of observations

# Plot species, phenophase group combos with a minimum number of observations
onsets_sub <- onsets %>%
  filter(n_obs_sppgroup >= min_obs)
spp <- sort(unique(onsets_sub$common_name))

# Plot onset dates for each species, 4 at a time
n_plots <- ceiling(length(spp) / 4)
for (i in 1:n_plots) {
  spps <- spp[(i * 4 - 3):(i * 4)]
  assign(paste0("onsets_p", i),
         ggplot(data = filter(onsets_sub, common_name %in% spps),
                aes(x = group_labels, y = firstyes, fill = group_labels)) +
           geom_boxplot(varwidth = TRUE) +
           ylim(min(onsets_sub$firstyes), max(onsets_sub$firstyes)) +
           facet_wrap(~ common_name, ncol = 1) +
           scale_fill_manual(values = color_vec) +
           labs(x = "", y = "Onset day of year") +
           coord_flip() +
           theme(legend.position = "none"))
}
for (i in 1:n_plots) {
  print(get(paste0("onsets_p",i)))
}
# May want to create species boxplots by region (ie, site "cluster") if sites
# are spread out geographically.

#- Trends in onset dates ------------------------------------------------------#

# For now, approaching this in a very simple way: 
  # assuming that the date of the first yes is the onset date
  # using the first yes of the year that was preceded by a no within X days
  # just looking at linear trends in a phenophase onset by spp or functional type

# Set the maximum number of days between first yes and prior no
prior_days_trends <- 14
onsett <- onsets %>% filter(lastno <= prior_days_trends)

# Identify the first year with >= 15 observations (across all species and 
# phenophases, and sites) and remove data before then
yr_obs <- count(onsett, yr) %>% arrange(yr)
min_yr <- yr_obs$yr[which(yr_obs$n >= 15)][1]
onsett <- onsett %>% filter(yr >= min_yr)

# Calculate the number of years we have data for each species-phenophase group 
# combination, and the number of observations we have for each 
# species-phenogroup-yr combination. (Note that if phenology likely differed 
# greatly by sites, then we would also want to evaluate sample sizes for 
# combinations at each site, but not worrying about that for now).
onsett <- onsett %>%
  group_by(common_name, phenogroup) %>%
  mutate(nyrs_sppgroup = n_distinct(yr)) %>%
  ungroup() %>%
  group_by(common_name, phenogroup, yr) %>%
  mutate(n_sppgroupyr = n()) %>%
  ungroup() %>%
  data.frame()

# Set minimum sample sizes (number of years with X observations) for each 
# species-phenophase group.
min_yrs <- 5
min_obs_per_yr <- 3
onsett <- onsett %>%
  group_by(common_name, phenogroup) %>%
  mutate(nyrs_minss = n_distinct(yr[n_sppgroupyr >= min_obs_per_yr])) %>%
  ungroup() %>%
  data.frame()
onsett <- onsett %>%
  filter(nyrs_minss >= min_yrs)

# Add in functional group information
onsett <- onsett %>%
  left_join(select(species, common_name, functional_type), by = "common_name")
onsett %>%
  group_by(functional_type) %>%
  summarize(n = n(),
            n_spp = n_distinct(common_name)) %>%
  data.frame()
# Combining all forbs and grasses so there's more than one species per 
# functional group
onsett <- onsett %>%
  mutate(func_group = case_when(
    functional_type %in% c("Forb", "Graminoid", "Semi-evergreen forb") ~
      "Forb or grass",
    .default = functional_type
  ))

# See how many observations we have for each species, phenogroup
onsett %>% select(common_name, phenogroup) %>% table()

# Evaluating trends in yellow birch
yebi <- filter(onsett, common_name == "yellow birch") %>%
  mutate(individual_id = factor(individual_id),
         site_id = factor(site_id))
count(yebi, yr) # 10 years
count(yebi, site_id) # 12 sites
count(yebi, individual_id) # 29 sites

ggplot(data = yebi, aes(x = yr, y = firstyes, color = site_id)) +
  geom_point() +
  geom_smooth(method = lm,  fullrange = FALSE, aes(group = 1), color = "black") +
  facet_grid(.~phenogroup) + 
  theme(legend.position = "none")

# Note that individual IDs are implicitly nested in sites (each individual_ID
# is only associated with one site), so we don't need to worry about how
# we specify the nestedness in the model formula. However, we may likely 
# run into issues with singular fits...

summary(lmer(firstyes ~ yr + (1|site_id) + (1|individual_id),
             data = filter(yebi, phenogroup == "lf")))
summary(lmer(firstyes ~ yr + (1|individual_id),
             data = filter(yebi, phenogroup == "lf")))
summary(lmer(firstyes ~ yr + (1|site_id),
             data = filter(yebi, phenogroup == "lf")))
# All have issues with singularity

summary(lm(firstyes ~ yr + site_id,
        data = filter(yebi, phenogroup == "lf")))

filter(yebi, phenogroup == "lf" & yr == 2014)
# Weird... 10 sites, 25 trees have onsets in 2014 but all the data look the same
# 15 onsets at doy 104 with last no 6 days prior
# 10 onsets at doy 152 with last no 12 days prior
# Do all volunteers go out on same day?
filter(yebi, phenogroup == "flo" & yr == 2015)
# Even more extreme: 23 open flower onsets, all but one on the same day

# Evaluating trends in bluebead
blue <- filter(onsett, common_name == "bluebead") %>%
  mutate(individual_id = factor(individual_id),
         site_id = factor(site_id))
count(blue, yr) # 10 years
count(blue, site_id) # 34 sites
count(blue, individual_id) # 34 individuals

ggplot(data = blue, aes(x = yr, y = firstyes, color = site_id)) +
  geom_point() +
  geom_smooth(method = lm,  fullrange = FALSE, aes(group = 1), color = "black") +
  facet_grid(.~phenogroup) + 
  theme(legend.position = "none")

summary(lmer(firstyes ~ yr + (1|site_id),
             data = filter(blue, phenogroup == "lf")))
summary(lmer(firstyes ~ yr + (1|site_id),
             data = filter(blue, phenogroup == "fl")))
summary(lmer(firstyes ~ yr + (1|site_id),
             data = filter(blue, phenogroup == "flo")))
summary(lmer(firstyes ~ yr + (1|site_id),
             data = filter(blue, phenogroup == "fr")))
summary(lmer(firstyes ~ yr + (1|site_id),
             data = filter(blue, phenogroup == "frr")))
# None of these are problematic

# Evaluating trends in mountain avens (lowest sample sizes of any species)
moav <- filter(onsett, common_name == "mountain avens") %>%
  mutate(individual_id = factor(individual_id),
         site_id = factor(site_id))
count(moav, yr) # 10 years
count(moav, site_id) # 8 sites
count(moav, individual_id) # 8 individuals

ggplot(data = moav, aes(x = yr, y = firstyes, color = site_id)) +
  geom_point() +
  geom_smooth(method = lm,  fullrange = FALSE, aes(group = 1), color = "black") +
  facet_grid(.~phenogroup) + 
  theme(legend.position = "none")

summary(lmer(firstyes ~ yr + (1|site_id),
             data = filter(moav, phenogroup == "fl"))) # singular fit (n = 28)
summary(lmer(firstyes ~ yr + (1|site_id), 
             data = filter(moav, phenogroup == "flo"))) # n = 36
summary(lmer(firstyes ~ yr + (1|site_id),
             data = filter(moav, phenogroup == "fr"))) # singular fit (n = 34)

# Run analyses for each functional group?
onsett %>% select(func_group, phenogroup) %>% table() 
onsett %>%
  group_by(func_group, common_name, phenogroup) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = phenogroup,
              values_from = n)
# Maybe should combine all trees since there are only 3 evergreen broadleaf spp?

# For tree species, can only do random intercepts (not intercepts and slopes)
tree_lf_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
                     data = filter(onsett, func_group != "Forb or grass",
                                   phenogroup == "lf"))
summary(tree_lf_sppr) # NS
tree_lfe_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
                      data = filter(onsett, func_group != "Forb or grass",
                                   phenogroup == "lfe"))
summary(tree_lfe_sppr) # S (positive)
tree_fl_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
                      data = filter(onsett, func_group != "Forb or grass",
                                    phenogroup == "fl"))
summary(tree_fl_sppr) # NS
tree_flo_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
                      data = filter(onsett, func_group != "Forb or grass",
                                    phenogroup == "flo"))
summary(tree_flo_sppr)
tree_fr_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
                     data = filter(onsett, func_group != "Forb or grass",
                                    phenogroup == "fr"))
summary(tree_fr_sppr) # NS (positive)
tree_frr_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
                      data = filter(onsett, func_group != "Forb or grass",
                                   phenogroup == "frr"))
summary(tree_frr_sppr) # NS (negative)

# Deciduous broadleaf 
# There are 5-6 species for all phenogroups except frr. I'll try random 
# effects, but not sure if that's a great idea.
# db_lf_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
#                    data = filter(onsett, func_group == "Deciduous broadleaf",
#                                  phenogroup == "lf"))
# summary(db_lf_sppr)
# db_lf_sppa <- lmer(firstyes ~ yr + common_name + (1|site_id),
#                    data = filter(onsett, func_group == "Deciduous broadleaf",
#                                  phenogroup == "lf"))
# summary(db_lf_sppa)
# db_lf_sppi <- lmer(firstyes ~ yr * common_name + (1|site_id),
#                    data = filter(onsett, func_group == "Deciduous broadleaf",
#                                  phenogroup == "lf"))
# summary(db_lf_sppi)
# db_lfe_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
#                     data = filter(onsett, func_group == "Deciduous broadleaf",
#                                  phenogroup == "lfe"))
# summary(db_lfe_sppr)
# db_fl_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
#                     data = filter(onsett, func_group == "Deciduous broadleaf",
#                                   phenogroup == "fl"))
# summary(db_fl_sppr)
# db_flo_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
#                     data = filter(onsett, func_group == "Deciduous broadleaf",
#                                  phenogroup == "flo"))
# summary(db_flo_sppr)
# db_fr_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
#                   data = filter(onsett, func_group == "Deciduous broadleaf",
#                                   phenogroup == "fr"))
# summary(db_fr_sppr)
# db_frr_sppa <- lmer(firstyes ~ yr + common_name + (1|site_id),
#                     data = filter(onsett, func_group == "Deciduous broadleaf",
#                                  phenogroup == "frr"))
# summary(db_frr_sppa)
# anova(db_frr_sppa)

# Forbs and grasses
forb_lf_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
                     data = filter(onsett, func_group == "Forb or grass",
                                   phenogroup == "lf"))
summary(forb_lf_sppr) # S (positive)
forb_fl_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
                     data = filter(onsett, func_group == "Forb or grass",
                                   phenogroup == "fl"))
summary(forb_fl_sppr) # NS
forb_flo_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
                     data = filter(onsett, func_group == "Forb or grass",
                                   phenogroup == "flo"))
summary(forb_flo_sppr) # S (negative)
forb_fr_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
                     data = filter(onsett, func_group == "Forb or grass",
                                   phenogroup == "fr"))
summary(forb_fr_sppr) # NS
forb_frr_sppr <- lmer(firstyes ~ yr + (1|common_name) + (1|site_id),
                     data = filter(onsett, func_group == "Forb or grass",
                                   phenogroup == "frr"))
summary(forb_frr_sppr) # NS

# Then run species-level analyses (plot all trends in a functional group
# in same figure in report).

onsett %>%
  group_by(func_group, common_name, phenogroup) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = phenogroup,
              values_from = n)

# Forb/grass species - open flowers
forb_flo_spp <- onsett %>%
  filter(func_group == "Forb or grass" & phenogroup == "flo") %>%
  pull(common_name) %>%
  unique()
for(i in 1:length(forb_flo_spp)) {
  df <- onsett %>%
    filter(common_name == forb_flo_spp[i] & phenogroup == "flo")
  assign(paste0("spp", i),
         lmer(firstyes ~ yr + (1|site_id), data = df))
  if (isSingular(get(paste0("spp", i)))) {
    message("Model for ", forb_flo_spp[i], " is Singular")
  }
}

yrpred <- data.frame(yr = seq(min(onsett$yr), max(onsett$yr), length = 100),
                     site_id = 0)

# PICK UP HERE
# See GLMM FAQ for prediction with CIs, though maybe for plotting species-level
# trends we won't need that (if multiple going on the same figure panel)



#- Yet to work on... ----------------------------------------------------------#

# Would be good to get an understanding of less "seasonal" species/regions for 
# which we're less interested in initial onset, and more interested in when 
# phenophases occur since that could happen over long period or multiple times
# a year. 

#NEXT (and related to point above): Do we want to also create visualizations for 
# weekly proportions of observations in a particular phenoclass/group? 
# (using GAMs or something else?). Yes
