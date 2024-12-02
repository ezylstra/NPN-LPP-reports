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
#TODO: Decide whether we need to keep class_id in there.

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

#TODO: Figure out if/where we want QA/QC (within indiv-date-observer) to happen

# For now, will keep record with more advanced phenophase or higher intensity or
# abundance value. Will do this by sorting observations and keeping only the 
# first
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
# a status = 0 could easily be a detectionissue. If we assume there are no false 
# positives, then we can simply use the maximum status value and maximum 
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
    si_sub_plants %>%
      group_by(common_name, class_id) %>%
      summarize(n_phenophases = length(unique(phenophase_id)),
                .groups = "keep") %>%
      data.frame()
  # For each species, may have more than 1-3 phenophases per pheno group
    si_sub_plants %>%
      group_by(common_name, pheno_group) %>%
      summarize(n_phenophases = length(unique(phenophase_id)),
                .groups = "keep") %>%
      data.frame()
  
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
  
  int_cat_sp <- si_sub_plants %>%
    group_by(common_name, phenophase_id, phenophase_description) %>%
    summarize(int_cats = length(unique(intensity_category_id)),
              .groups = "keep") %>%
    data.frame()
  # yellow birch, 483 = Leaves, 3 int cats
  count(filter(si_sub_plants, common_name == "yellow birch", phenophase_id == "483"),
        intensity_category_id) # 40, 73, or NA
  # No, but I think that's because the counts include intensity categories of NA
  # and because the intensity category may have shifted over time:
  count(filter(si_sub_plants, common_name == "yellow birch", phenophase_id == "483"),
        intensity_category_id, yr) # 40 used in 2014-2015; 73 used in 2016-2023
  test <- si_sub_plants %>%
    group_by(individual_id, phenophase_id, yr) %>%
    summarize(int_cat = length(unique(intensity_category_id[!is.na(intensity_category_id)])),
              .groups = "keep") %>%
    data.frame()
  count(test, int_cat) 
  # In a single year, only one intensity category (if any) for that plant, phenophase

# All this means that we put data in wide form, with 2-3 columns for each 
# plant phenophase: status, intensity (& intensity type?)

# Putting data in wide form:
# Can use phenophase class_id instead of phenophase_id because there's a maximum
# of one phenophase per class for each species.

# Create a table with info/names for each plant phenophase class    
pl_pheno_classes <- data.frame(
  class_id = 1:13,
  class_id2 = str_pad(1:13, width = 2, pad = 0),
  pheno_class = c("initial_leaf", "young_leaf", "leaf", "color_leaf", 
                  "fall_leaf", "flower", "open_flower", "pollen", "end_flower",
                  "fruit", "unripe", "ripe", "drop_fruit"),
  pheno_class_g = c(paste0("leaves", 1:3),
                    paste0("leaf_end", 1:2),
                    paste0("flower", 1),
                    paste0("flower_open", 1:2),
                    paste0("flower_end", 1),
                    paste0("fruit", 1:2),
                    paste0("fruit_ripe", 1:2))
)
pl_pheno_classes <- pl_pheno_classes %>%
  mutate(pheno_group = str_sub(pheno_class_g, end = -2)) 
  # %>%
  # mutate(phenoG = paste0(pheno_group, "G"))

# Put in wide form and label columns appropriately
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
# Add in columns for all pheno classes that aren't already in there
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
first_status_col <- which(grepl("status", colnames(plant_obs)))[1]
plant_obs <- plant_obs %>%
  select(colnames(plant_obs)[1:(first_status_col - 1)], all_of(all_cols))
# Put phenophase class name in column name
for (i in 1:nrow(pl_pheno_classes)) {
  colnames(plant_obs) <- str_replace(colnames(plant_obs), 
                                pl_pheno_classes$class_id2[i], 
                                pl_pheno_classes$pheno_class[i])
}

# Resolve any inconsistencies with status
  # If young leaves = 1, leaves = 1 (ok if young = 1, leaves = NA)
  # If color leaves = 1, leaves = 1
  # If open flower = 1, flower = 1
  # If pollen = 1, flower = 1
  # If unripe fruit = 1, fruit = 1
  # If ripe fruit = 1, fruit = 1
  # For now, leaving instances where more specific category = 1 and general
  # category (leaves, flower, fruit) = NA
  #TODO: probably want to provide this info to others and ask about NAs

# First, see how often these issues come up?
count(plant_obs, status_young_leaf, status_leaf) 
count(plant_obs, status_color_leaf, status_leaf) 
count(plant_obs, status_open_flower, status_flower) 
count(plant_obs, status_pollen, status_flower) 
count(plant_obs, status_unripe, status_fruit) 
count(plant_obs, status_ripe, status_fruit) 

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
    (status_unripe == 1 & status_fruit == 0) ~ 1,
    (status_ripe == 1 & status_fruit == 0) ~ 1, 
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
# observations of that plant on that date:
summary(plant_obs$observer_rank[plant_obs$n_observations > 1])
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

#- Tables: Site summaries ------------------------------------------------------#

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
# for particular sites or species of interest?













#- Aggregate status information within pheno group? ---------------------------#

# Aggregate status information within pheno group (doesn't make sense to 
# aggregate intenstiy values, since they may not always be the same type)


groups <- unique(pl_pheno_classes$pheno_group)
for (i in groups) {
  classes <- pl_pheno_classes$pheno_class[pl_pheno_classes$pheno_group == i]
  sclasses <- paste0("status_", classes)
  if (length(sclasses) == 1) {
    plant_obs2[,paste0("Gstatus_", i)] <- plant_obs2[,sclasses]
  } else {
    plant_obs2[,paste0("Gstatus_", i)] <- apply(plant_obs2[,sclasses], 
                                                1, na_max)
  }
}  





#- OLD STUFF ------------------------------------------------------------------#

#- Quick overview of data -----------------------------------------------------#
# Data dimensions
dim(si) # 528,810
count(si, site_name); sort(count(si, site_name)$n) # 105 sites (5 with < 50 records)
count(si, yr) # 2012-2023, but only 7 and 149 records in 2012, 2013 
count(filter(si, kingdom == "Plantae"), common_name)   # 21 plant spp ('ohi'a lehua with 28, rest with > 1900)
count(filter(si, kingdom == "Animalia"), common_name)  # 7 animal spp (but monarch only 1 record)

# How many phenophases recorded across, within species?
count(si, phenophase_description) # 41 phenophases
count(si, pheno_group)            # 9 phenophase "groups"
spp_ph <- si %>%
  group_by(kingdom, common_name) %>%
  summarize(n_phenophases = length(unique(phenophase_description)),
            n_phenogroups = length(unique(pheno_group)),
            .groups = "keep") %>%
  data.frame()
spp_ph 
# For animals: 7-13 phenophases, 1-3 phenogroups per species (excl monarchs)
# For plants: 6-11 phenophases, 5-6 phenogroups per species

# Looking at examples where we have a ton of different phenophases observed
count(filter(obs, n == n_ph), n, n_ph) 
# Up to 11 phenophases recorded as part of one observation 
head(filter(obs, n == 11, n_ph == 11))
count(filter(obs, n == 11, n_ph == 11), common_name)
# 1 animal (BTBW); 3 plants (maples, birch)

# Looking at birds:
filter(si, 
       common_name == "black-throated blue warbler",
       site_name == "Lonesome Lake Hut (LONE1)",
       obsdate == "2014-06-14")
# Live individuals, Feeding, Fruit/seed consumption, Insect consumption,
# Flower visitation, Calls or song (birds), Singing individuals (birds), 
# Mating (male on top), Nest building (birds), Dead individuals, 
# Individuals at a feeding station

# Note that if they have calls or singing, may still have 0 for live individuals
# For these summaries, can probably create a new "phenophase" that is 
# Activity = 1 if any of the phenophase status (Except for dead) are 1, 0 if not
# Remove records of dead individuals since they don't provide much info

# Looking at trees
head(filter(obs, n == 11, n_ph == 11, kingdom == "Plantae"))
filter(si, 
       common_name == "red maple",
       site_name == "Ethan Pond Trail (EP2)",
       obsdate == "2020-06-16")$phenophase_description
# Breaking leaf buds, Leaves, Increasing leaf size, Colored leaves,
# Falling leaves, Flowers or flower buds, Open flowers, 
# Pollen release (flowers), Fruits, Ripe fruits, Recent fruit or seed drop

