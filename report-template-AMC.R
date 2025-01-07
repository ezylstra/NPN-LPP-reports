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
require(ggforce) # Nice facet options for ggplot
require(terra)
require(lme4)
require(ggeffects) # Predictions & plotting for mixed-effect models
# require(lmerTest) # Don't want to load this (but will call in script)

# require(geosphere) # for optional site clustering
# require(fpc) # for optional site clustering
# require(RColorBrewer)

rm(list = ls())

#- Specify data of interest ---------------------------------------------------#
# Name and nickname of LPP
lpp <- "Appalachian Mountain Club"
lpp_short <- "amc"
# lpp <- "Earthwise Aware"
# lpp_short <- "ewa"

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

# Remove any observations of 'ohi'a lehua if present (would be an error)
si <- si %>%
  filter(common_name != "'ohi'a lehua")

# Append phenophase "groups" (csv created in phenophases-intensities.R)
pheno_list <- read.csv("phenophases.csv")
si <- left_join(si, select(pheno_list, phenophase_id, class_id, pheno_group),
                by = "phenophase_id")
# NOTE: for reports, may refer to "leaf senescence" as "colored leaves"

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
  # For each species, there is only 1 phenophase per pheno class EXCEPT:
  # common lilac: Open flowers (205) and Full flowering (206) both in 
  # phenophase class 7 (Open flowers or pollen cones)
    # phaseperclass <- si_sub_plants %>%
    #   group_by(common_name, class_id) %>%
    #   summarize(n_phenophases = length(unique(phenophase_id)),
    #             .groups = "keep") %>%
    #   data.frame()
    # filter(phaseperclass, n_phenophases > 1)
    # count(filter(si, common_name == "common lilac" & class_id == 7), 
    #       phenophase_id, phenophase_description)
    # filter(pheno_list, class_id == 7)
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
    # count(filter(si_sub_plants, common_name == "American chestnut", phenophase_id == "371"),
    #       intensity_category_id) 
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

# Summary table for phenogroups:
pheno_list <- pheno_list %>%
  mutate(pheno_group = replace(pheno_group, 
                               pheno_group == "Leaf senescence", 
                               "Colored leaves"))
pgs_pl <- pheno_list %>%
  filter(class_id %in% 1:13) %>%
  group_by(pheno_group, class_id, class_name) %>%
  summarize(n_phenophases = n(), .groups = "keep") %>%
  data.frame() %>%
  arrange(class_id)
pgs_an <- pheno_list %>%
  filter(!class_id %in% 1:13) %>%
  group_by(pheno_group) %>%
  summarize(n_classes = n_distinct(class_name),
            n_phenophases = n()) %>%
  data.frame()
# write.table(pgs_pl, "clipboard", sep = "\t", row.names = FALSE)
# write.table(pgs_an, "clipboard", sep = "\t", row.names = FALSE)

# Put in wide form
plant_obs <- si_sub_plants %>%
  filter(phenophase_id != 206) %>%
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
         obs_plant = n_plant_obs) %>%
  arrange(desc(years), desc(obs_total))				 
#TODO: Decide whether we want filters for number of years (>= 3?), number
# of observations (> 100?), or number of species?

sitest_sub <- sitest %>%
  select(-c(site_id, yr_first, yr_last)) %>%
  mutate(site_name = str_trim(site_name, side = "both"))

# write.table(sitest_sub, "clipboard", sep = "\t", row.names = FALSE)

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

# write.table(yearst, "clipboard", sep = "\t", row.names = FALSE)

#- Table: Species summaries ---------------------------------------------------#
# Recreating species table after removing years with few observations

yrs_subset <- yearst$yr[yearst$obs_total > 25]

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

# Use the plant/animal_obs2 dataframes to calculate frequency (number of days
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

# For animals, just using live observations
freq_an <- animal_obs2 %>%
  filter(pheno_group == "Live observation") %>%
  select(site_id, common_name, individual_id, day_of_year, yr) %>%
  distinct() %>%
  filter(yr %in% yrs_subset) %>%
  arrange(common_name, individual_id, yr, day_of_year)
freq_an$interval <- NA
for (i in 2:nrow(freq_an)) {
  if (freq_an$individual_id[i] == freq_an$individual_id[i - 1] &
      freq_an$yr[i] == freq_an$yr[i - 1]) {
    freq_an$interval[i] <- freq_an$day_of_year[i] - freq_an$day_of_year[i - 1]
  } 
}
freqs_an <- freq_an %>%
  group_by(common_name, individual_id, yr) %>%
  summarize(nobs = n(), 
            interval = median(interval, na.rm = TRUE),
            .groups = "keep") %>%
  data.frame() %>%
  group_by(common_name) %>%
  summarize(nobs_mn = round(mean(nobs)),
            interval_mn = round(mean(interval, na.rm = TRUE), 1)) %>%
  data.frame()

freqs <- rbind(freqs_pl, freqs_an) %>%
  rename(nobs_mn_per_yrind = nobs_mn,
         interval_mn_per_yrind = interval_mn)
species2 <- species2 %>%
  left_join(freqs, by = "common_name") %>%
  mutate(functional_type = factor(functional_type, 
                                  levels = c("Deciduous broadleaf",
                                             "Evergreen broadleaf",
                                             "Evergreen conifer",
                                             "Pine",
                                             "Evergreen forb",
                                             "Forb",
                                             "Semi-evergreen forb",
                                             "Graminoid",
                                             "Amphibian",
                                             "Bird",
                                             "Insect", 
                                             "Mammal",
                                             "Reptile"))) %>%
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
# TODO: make this easier by putting things into long form first?
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
    phenogroup == "lfe" ~ "Colored leaves",
    phenogroup == "fl" ~ "Flowers",
    phenogroup == "flo" ~ "Open flowers",
    phenogroup == "fle" ~ "Flower senescence",
    phenogroup == "fr" ~ "Fruits",
    phenogroup == "frr" ~ "Ripe fruits"
  )) %>%
  mutate(group_labels = factor(group_labels,
                               levels = c("Colored leaves",
                                          "Ripe fruits",
                                          "Fruits",
                                          "Flower senescence",
                                          "Open flowers",
                                          "Flowers",
                                          "Leaves")))

# Assign colors for each phenophase group
color_vec <- c("#b2df8a",   # Leaves
               "#fdbf6f",   # Flowers
               "#ff7f00",   # Open flowers
               "#b2beb5",   # Flower senescence
               "#cab2d6",   # Fruits
               "#6a3d9a",   # Ripe fruits
               "#8c510a")   # Colored leaves
# Name vector so colors are consistent across figures (in case not all levels 
# are present in all figures)
names(color_vec) <- rev(levels(onsets$group_labels))

# Add sample sizes for each species-phenogroup combination to onsets
onsets <- onsets %>%
  group_by(common_name, phenogroup) %>%
  mutate(n_obs_sppgroup = n()) %>%
  ungroup() %>%
  data.frame()

# Add in functional group information
onsets <- onsets %>%
  left_join(select(species, common_name, functional_type), by = "common_name")
onsets %>%
  group_by(functional_type) %>%
  summarize(n = n(),
            n_spp = n_distinct(common_name)) %>%
  data.frame()
onsets <- onsets %>%
  mutate(common_name_full = paste0(common_name, " (", 
                                   str_to_lower(functional_type), ")"))
# Combining all forbs and grasses so there's more than one species per 
# functional group
onsets <- onsets %>%
  mutate(func_group = case_when(
    functional_type %in% c("Forb", "Graminoid", 
                           "Semi-evergreen forb", "Evergreen forb") ~
      "Forb or grass",
    .default = functional_type
  ))
onsets <- onsets %>%
  arrange(func_group, functional_type, common_name, .locale = "en")

# Set minimum number of observations per species-phenogroup to summarize onset 
# dates
min_obs <- 15

# Plot onset for one phenogroup by species (lumping years, sites together).
# Functional groups in facets

phenogs <- c("lf", "fl", "flo", "fr", "frr", "lfe")
phenogs2 <- c("leaves", "flowers", "open flowers", "fruits", "ripe fruits",
              "colored leaves")
onsets_plot_width <- 6.5
onsets_plot_heightmax <- 6.5
total_nspp <- n_distinct(onsets$common_name[onsets$n_obs_sppgroup >= min_obs])

for (i in 1:length(phenogs)) {
  
  phenog <- phenogs[i]
  col <- color_vec[str_to_sentence(phenogs2[i])]
  
  # Create nice labels for dates on Y axes in plots
  plotdat <- filter(onsets, phenogroup == phenog, n_obs_sppgroup >= min_obs)
  mindoy <- min(plotdat$firstyes)
  maxdoy <- max(plotdat$firstyes)
  doys <- seq(mindoy, maxdoy)
  plotdates <- as.Date(paste(2023, doys, sep = "-"), "%Y-%j")
  date1ind <- which(day(plotdates) == 1)
  doybreaks <- doys[date1ind]
  datebreaks <- plotdates[date1ind] %>% format("%m/%d")
  datelims <- c(min(doys), max(doys) + 10)
  
  # Function to help add sample sizes to plot (95% of the way across date axis)
  n_fun <- function(x) {
    return(data.frame(y = 0.99 * datelims[2],
                      label = length(x)))
  }
  
  # Calculate height of figure
  plot_nspp <- n_distinct(plotdat$common_name)
  onsets_plot_height <- onsets_plot_heightmax * (plot_nspp / total_nspp) + 0.5
  
  assign(paste0("onsets_plot_", phenog),
   ggplot(data = plotdat, aes(x = common_name, y = firstyes)) +
     geom_boxplot(fill = col) + 
     stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5,
                  family = "sans", size = 3) +
     labs(title = paste0("First observation of ", phenogs2[i])) +
     ggforce::facet_col(~func_group, scales = "free_y", space = "free") +
     scale_y_continuous(limits = datelims, breaks = doybreaks,
                        labels = datebreaks) +
     coord_flip() +
     theme_bw() +
     theme(axis.title = element_blank(),
           axis.text.x = element_text(size = 9), 
           axis.text.y = element_text(size = 9),
           strip.text = element_text(size = 9))
  )

  # ggsave(paste0("output/onsets-plot-", phenog, "-", lpp_short, ".png"),
  #        get(paste0("onsets_plot_", phenog)),
  #        width = onsets_plot_width,
  #        height = onsets_plot_height,
  #        units = "in",
  #        dpi = 600)
  
}
# Saved a ggplot object for each phenogroup: onsets_plot_PG
# and saved each to file

# Plot species, phenophase group combos with a minimum number of observations
# Only include spp that have sufficient sample sizes for 3 or more phenophases
onsets_sub <- onsets %>%
  filter(n_obs_sppgroup >= min_obs) %>%
  group_by(common_name) %>%
  mutate(n_phenogroups = n_distinct(phenogroup)) %>%
  ungroup() %>%
  data.frame() %>%
  filter(n_phenogroups > 2)
  
spp <- unique(onsets_sub$common_name_full)

# Create nice labels for dates on X axes in plots
mindoy <- min(onsets_sub$firstyes)
maxdoy <- max(onsets_sub$firstyes)
doys <- seq(mindoy, maxdoy)
plotdates <- as.Date(paste(2023, doys, sep = "-"), "%Y-%j")
date1ind <- which(day(plotdates) == 1)
doybreaks <- doys[date1ind]
datebreaks <- plotdates[date1ind] %>% format("%m/%d")
datelims <- c(min(doys), max(doys) + 10)

# Plot onset dates for each species, 4 at a time
n_plots <- ceiling(length(spp) / 4)
for (i in 1:n_plots) {
  spps <- spp[(i * 4 - 3):(i * 4)]
  plotdat <- filter(onsets_sub, common_name_full %in% spps)
  assign(paste0("onsets_p", i),
         ggplot(data = plotdat,
                aes(x = group_labels, y = firstyes, fill = group_labels)) +
           geom_boxplot() +
           stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, 
                        family = "sans", size = 3) +
           scale_y_continuous(limits = datelims, breaks = doybreaks, 
                              labels = datebreaks) +
           facet_wrap(~ factor(common_name_full, levels = spps), ncol = 1) +
           scale_fill_manual(values = color_vec) +
           labs(x = "", y = "First day observed") +
           coord_flip() +
           theme_bw() +
           theme(legend.position = "none",
                 axis.title = element_blank(),
                 axis.text.x = element_text(size = 9), 
                 axis.text.y = element_text(size = 9),
                 strip.text = element_text(size = 9))
         )
  # ggsave(paste0("output/onsets-plot-spp-", lpp_short, "-", i, ".png"),
  #        get(paste0("onsets_p", i)),
  #        width = onsets_plot_width,
  #        height = 7.5,
  #        units = "in",
  #        dpi = 600)
}

#- Trends in onset dates ------------------------------------------------------#

# For now, approaching this in a very simple way: 
  # assuming that the date of the first yes in series is the onset date
  # using the first yes of the year that was preceded by a no within X days
  # just looking at linear trends in phenophase onset by spp or functional type

# Set the maximum number of days between first yes and prior no
prior_days_trends <- 14
onsett <- onsets %>% filter(lastno <= prior_days_trends)

# Identify the first year with >= 15 observations (across all species, 
# phenophases, and sites) and remove data before then
yr_obs <- count(onsett, yr) %>% arrange(yr)
min_yr <- yr_obs$yr[which(yr_obs$n >= 15)][1]
onsett <- onsett %>% filter(yr >= min_yr)

# Calculate the number of years we have data for each species-phenophase group 
# combination, and the number of observations we have for each 
# species-phenogroup-yr combination. (Note that if phenology differed greatly
# between sites, then we would also want to evaluate sample sizes for 
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

# Run analyses for each functional group
onsett %>% select(func_group, phenogroup) %>% table() 
onsett %>%
  group_by(func_group, common_name, phenogroup) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = phenogroup,
              values_from = n)

# Running two models for each functional group - phenophase: 
# 1) RS = models that specify random slopes for each species and random 
#    intercepts for each site, species, individual. If these ran without issues
#    we wouldn't run the 2nd model. However, they often have issues with
#    singular fits where the the variance component associated with random 
#    slopes is 0. Still useful to get that overall trend for the functional 
#    group that accounts for repeated measures. Note that we only attempted 
#    these models if there were 4 or more species. If not, we excluded species
#    from the random part of the model.
# 2) FS = models that have yr * species as fixed effects and site, individual 
#    random intercepts.
# Using ggeffects package to make and plot predictions
# https://strengejacke.github.io/ggeffects/articles/introduction_randomeffects.html
# With predict_response(), default is to produce 95% confidence intervals, but
# could get prediction intervals with interval = "prediction". Function
# generates predictions for an average/typical species, site, individual; or, if 
# using the type = "random" argument, get predictions for each level of 
# specified random effect.Can plot results directly, or convert to data.frame 
# and then use ggplot()

# Create nice labels for phenophase groups
pl_pheno_classes <- pl_pheno_classes %>%
  mutate(phenogroup = case_when(
    group_code == "lf" ~ "Leaves",
    group_code == "lfe" ~ "Colored leaves",
    group_code == "fl" ~ "Flowers",
    group_code == "flo" ~ "Open flowers",
    group_code == "fle" ~ "Flower senescence",
    group_code == "fr" ~ "Fruit",
    group_code == "frr" ~ "Ripe fruit"
  ))
pl_pheno_match <- pl_pheno_classes %>%
  select(group_code, phenogroup) %>%
  rename(phenogroup_f = phenogroup) %>%
  rename(phenogroup = group_code) %>%
  mutate(phenogroup_f = factor(phenogroup_f,
                               levels = c("Leaves",
                                          "Flowers",
                                          "Open flowers",
                                          "Flower senescence",
                                          "Fruit",
                                          "Ripe fruit",
                                          "Colored leaves"))) %>%
  distinct()
onsett <- onsett %>%
  left_join(pl_pheno_match, by = "phenogroup")

func_groups_df <- data.frame(
  func_group = unique(onsett$func_group),
  func_group2 = str_replace_all(unique(onsett$func_group), " ", "_"))

# Create function to identify models that failed to converge
is_warning_generated <- function(m) {
  df <- summary(m)
  !is.null(df$optinfo$conv$lme4$messages) && 
    sum(grepl('failed to converge', df$optinfo$conv$lme4$messages)) > 0
}

for (i in 1:nrow(func_groups_df)) {
  fgroup <- func_groups_df$func_group[i]
  phenogs <- sort(unique(onsett$phenogroup[onsett$func_group == fgroup]))

  # Run mixed models
  for (phenog in phenogs) {
    dat <- filter(onsett, func_group == fgroup, phenogroup == phenog)
    if (n_distinct(dat$common_name) >= 4) {
      suppressWarnings(
        RSm <- lmer(firstyes ~ yr + (0 + yr|common_name) + (1|common_name) + 
                    (1|site_id) + (1|individual_id), data = dat)
      )
      if (is_warning_generated(RSm)) {
        RSm <- lmer(firstyes ~ yr + (0 + yr|common_name) + (1|common_name) + 
                      (1|site_id) + (1|individual_id), data = dat,
                    control = lmerControl(optimizer = "Nelder_Mead"))
        if (is_warning_generated(RSm)) {
          print(paste0("Model failed to converge for RSm: ", fgroup, ", ", phenog))
        }
      }
    } else {
      suppressWarnings(
        RSm <- lmer(firstyes ~ yr + (1|site_id) + (1|individual_id), data = dat)
      )
      if (is_warning_generated(FSm)) {
        RSm <- lmer(firstyes ~ yr + (1|site_id) + (1|individual_id), data = dat,
                    control = lmerControl(optimizer = "Nelder_Mead"))
        if (is_warning_generated(RSm)) {
          print(paste0("Model failed to converge for RSm: ", fgroup, ", ", phenog))
        }
      }
    }
    suppressWarnings(
      FSm <- lmer(firstyes ~ yr * common_name + (1|site_id) + (1|individual_id),
                  data = dat)
    )
    if (is_warning_generated(FSm)) {
      FSm <- lmer(firstyes ~ yr * common_name + (1|site_id) + (1|individual_id),
                  data = dat, control = lmerControl(optimizer = "Nelder_Mead"))
    }
    RSp <- data.frame(predict_response(RSm, terms = "yr")) %>%
      mutate(group = "mean", 
             func_group = fgroup,
             phenogroup = phenog) %>%
      rename(yr = x, spp = group) %>%
      mutate(modeled_spp = if_else(n_distinct(dat$common_name) >= 4, 
                                   "random", "none"))
    FSp <- data.frame(predict_response(FSm, terms = c("yr", "common_name"))) %>%
      mutate(func_group = fgroup,
             phenogroup = phenog) %>%
      rename(yr = x, spp = group) %>%
      mutate(modeled_spp = "fixed")
    preds_new <- rbind(RSp, FSp)
    if (phenog == phenogs[1]) {
      preds <- preds_new
    } else {
      preds <- rbind(preds, preds_new)
    }
    # Remove years with no observations
    sppdat <- dat %>%
      group_by(common_name) %>%
      summarize(minyr = min(yr),
                maxyr = max(yr)) %>%
      data.frame()
    for (spp in sppdat$common_name) {
      preds$predicted[preds$spp == spp & preds$yr < sppdat$minyr[sppdat$common_name == spp]] <- NA 
      preds$conf.low[preds$spp == spp & preds$yr < sppdat$minyr[sppdat$common_name == spp]] <- NA 
      preds$conf.high[preds$spp == spp & preds$yr < sppdat$minyr[sppdat$common_name == spp]] <- NA
      preds$std.error[preds$spp == spp & preds$yr < sppdat$minyr[sppdat$common_name == spp]] <- NA 
      preds$predicted[preds$spp == spp & preds$yr > sppdat$maxyr[sppdat$common_name == spp]] <- NA 
      preds$conf.low[preds$spp == spp & preds$yr > sppdat$maxyr[sppdat$common_name == spp]] <- NA 
      preds$conf.high[preds$spp == spp & preds$yr > sppdat$maxyr[sppdat$common_name == spp]] <- NA 
      preds$std.error[preds$spp == spp & preds$yr > sppdat$maxyr[sppdat$common_name == spp]] <- NA 
    }
  
    # Extract estimate of trends for functional groups
    suppressWarnings(RSmt <- lmerTest::as_lmerModLmerTest(RSm))
    trend <- summary(RSmt)$coefficients["yr", ]
    re_names <- names(summary(RSm)$varcor)
    # est, SE, df, t, P  
    
    trends_new <- data.frame(
      func_group = fgroup,
      phenogroup = phenog,
      nobs = nrow(dat),
      nspp = n_distinct(dat$common_name),
      minyr = min(dat$yr),
      maxyr = max(dat$yr),
      random_slopes = ifelse("common_name.1" %in% re_names, 
                             "Yes", "No"),
      random_ints = str_c(str_subset(re_names, ".1", negate = TRUE), 
                          collapse = ", "),
      beta = round(trend["Estimate"], 5),
      se = round(trend["Std. Error"], 5),
      df = round(trend["df"]),
      t = round(trend["t value"], 3),
      P = round(trend["Pr(>|t|)"], 4),
      row.names = NULL
    )
    if (phenog == phenogs[1]) {
      trends <- trends_new
    } else {
      trends <- rbind(trends, trends_new)
    }
  }
  # Will likely get some warnings about singular fits (and maybe other stuff,
  # but hopefully none of it is prohibitive)

  preds <- preds %>%
    left_join(pl_pheno_match, by = "phenogroup")
  
  # Create nice labels for dates on Y axes in plots
  plotdat <- filter(onsett, func_group == fgroup)
  mindoy <- min(plotdat$firstyes)
  maxdoy <- max(plotdat$firstyes)
  doys <- seq(mindoy, maxdoy)
  plotdates <- as.Date(paste(2023, doys, sep = "-"), "%Y-%j")
  date1ind <- which(day(plotdates) == 1)
  doybreaks <- doys[date1ind]
  datebreaks <- str_glue("{month(plotdates[date1ind])}/{day(plotdates[date1ind])}") %>%
    as.character()

  # Create text for each panel (beta estimates and P-values)
  plotdat_summary <- plotdat %>%
    group_by(phenogroup) %>%
    summarize(mindoy = min(firstyes),
              maxdoy = max(firstyes)) %>%
    data.frame() %>%
    mutate(rangedoy = maxdoy - mindoy,
           label1y = maxdoy,
           label2y = label1y - (0.065 * rangedoy)) %>%
    select(-c(mindoy, maxdoy, rangedoy))
  ann_text <- data.frame(yr = max(preds$yr),
                         phenogroup = trends$phenogroup) %>%
    left_join(distinct(select(preds, phenogroup, phenogroup_f)), 
              by = "phenogroup") %>%
    left_join(plotdat_summary, by = "phenogroup")
  ann_text <- ann_text %>%
    left_join(select(trends, phenogroup, beta, P), by = "phenogroup") %>%
    mutate(beta2 = formatC(beta, digits = 2, format = "f"),
           P2 = formatC(P, digits = 2, format = "f"),
           Pstr = if_else(P2 == "0.00", " < 0.01", paste0(" = ", P2))) %>%
    mutate(label1 = paste0("slope = ", beta2, " days/yr"),
           label2 = paste0("P", Pstr))
  
  text_size <- 8
  ci_alpha <- 0.17
  
  # Save ggplot object
  assign(paste0("trendsplot_", func_groups_df$func_group2[i]),
         ggplot(data = filter(preds, spp != "mean"), aes(x = yr, y = predicted)) +
           geom_line(aes(color = spp)) +
           geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = spp), 
                       alpha = ci_alpha) +
           geom_point(data = plotdat,
                      aes(x = yr, y = firstyes, color = common_name),
                      position = position_dodge(width = 0.3),
                      size = 0.6, alpha = 0.4) +
           geom_line(data = filter(preds, spp == "mean"),
                     aes(x = yr, y = predicted), linewidth = 1) +
           facet_wrap(.~phenogroup_f, scales = "free_y") +
           geom_text(data = ann_text,
                     aes(y = label1y, x = yr, label = label1, hjust = 1, vjust = 1,
                         family = "sans"), size = text_size/.pt) +
           geom_text(data = ann_text,
                     aes(y = label2y, x = yr, label = label2, hjust = 1,
                         vjust = 1, family = "sans"), size = text_size/.pt) +
           labs(y = "First observation date", fill = "", color = "") +
           scale_y_continuous(breaks = doybreaks, labels = datebreaks) +
           theme_bw() +
           theme(axis.title.x = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 legend.position = "bottom",
                 axis.text.x = element_text(size = text_size), 
                 axis.text.y = element_text(size = text_size),
                 strip.text = element_text(size = text_size),
                 axis.title.y = element_text(size = text_size + 1),
                 legend.title = element_text(size = text_size),
                 legend.text = element_text(size = text_size))
  )
  # ggsave(paste0("output/trends-plot-", func_groups_df$func_group2[i],
  #               "-", lpp_short, ".png"),
  #        get(paste0("trendsplot_", func_groups_df$func_group2[i])),
  #        width = 6.5,
  #        height = 6.5,
  #        units = "in",
  #        dpi = 600)
  
  assign(paste0("preds_", func_groups_df$func_group2[i]), preds)
  assign(paste0("trends_", func_groups_df$func_group2[i]), trends)
}

# For each functional group, have:
# trendsplot_FUNC = ggplot object with a panel for each phenophase
# preds_FUNC = table with trend predictions (and 95% CIs)
# trends_FUNC = table with trend estimates

fg_trends <- get(paste0("trends_", func_groups_df$func_group2[1]))
for (i in 2:nrow(func_groups_df)) {
  fg_trends <- rbind(fg_trends, 
                     get(paste0("trends_", func_groups_df$func_group2[i])))
}
fg_trends <- fg_trends %>%
  mutate(random_ints = str_remove_all(random_ints, "_id")) %>%
  mutate(random_ints = str_replace_all(random_ints, "common_name", "species"))
# write.table(fg_trends, "clipboard", sep = "\t", row.names = FALSE)

#- Characterizing phenophase periods ------------------------------------------#
# We're not always just interested in onsets....

df <- plg_status

# Add week to the data so we can calculate weekly proportions
# We'll also create wk_doy columns to assign each week with a day of the year:
  # wk_doy1 = start of each week (eg, date for week 1 would be Jan 1)
  # wk_doy4 = middle of each week (eg, date for week 1 would be Jan 4)
df <- df %>%
  # remove phenophase that is always NA
  select(-status_fle) %>%
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

# Summarize data available across all plants, years, phenophases
propdat_wk <- df %>%
  group_by(common_name, phenogroup, wk) %>%
  summarize(n_obs = n(), 
            n_indiv = n_distinct(individual_id),
            n_yrs = n_distinct(yr),
            n_sites = n_distinct(site_id),
            .groups = "keep") %>%
  data.frame()
propdat_wk_spp <- propdat_wk %>%
  group_by(common_name, phenogroup) %>%
  summarize(mn_obs_per_wk = round(mean(n_obs), 1),
            max_n_yrs = max(n_yrs),
            max_n_sites = max(n_sites),
            .groups = "keep") %>%
  data.frame()

# Exclude any species, phenophase combinations where the mean number of
# observations per week is < 10
spp_ph_keep <- propdat_wk_spp %>%
  filter(mn_obs_per_wk >= 10) %>%
  mutate(spp_ph = paste0(common_name, "_", phenogroup))
df <- df %>%
  mutate(spp_ph = paste0(common_name, "_", phenogroup)) %>%
  filter(spp_ph %in% spp_ph_keep$spp_ph)

# Calculate weekly proportions of yeses for each species, phenophase
wkprops <- df %>%
  group_by(common_name, phenogroup, wk, wk_date1, wk_doy1, 
           wk_date4, wk_doy4) %>%
  summarize(n_obs = n(),
            n_yes = sum(status),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(prop_yes = n_yes / n_obs) %>%
  left_join(distinct(select(onsett, phenogroup, phenogroup_f)), 
            by = "phenogroup")

# Remove any weekly proportions that are based on < 5 observations (over 
# all years, individuals)
wkprops <- wkprops %>%
  filter(n_obs >= 5)

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

  # Subset if LPP only monitored over part of the year
  lpp_month1 <- month(min(wkprops$wk_date4))
  lpp_month2 <- month(max(wkprops$wk_date4))
  x_lab <- x_lab[month(x_lab) %in% lpp_month1:lpp_month2]
  tickind <- which(month(x_tick) %in% lpp_month1:lpp_month2)
  if (last(tickind) != 13) {
    tickind <- c(tickind, 1 + last(tickind))
  }
  x_tick <- x_tick[tickind]

  x_lab_doy <- yday(x_lab)
  x_tick_doy <- yday(x_tick)
  n_x_tick <- length(x_tick)
  month_labels <- month.abb[month(x_lab)]
  
  color_vec2 <- color_vec 
  names(color_vec2)[5:6] <- c("Fruit", "Ripe fruit")
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
  darker_gr <- rgb2col(round(col2rgb(color_vec2["Leaves"]) * 0.8))
  color_vec2["Leaves"] <- darker_gr
  
wkprops <- wkprops %>% 
  left_join(select(species, common_name, functional_type), 
            by = "common_name") %>%
  mutate(func_group = case_when(
    functional_type %in% c("Forb", "Graminoid", 
                           "Semi-evergreen forb", "Evergreen forb") ~
      "Forb or grass",
    .default = functional_type
  )) %>%
  mutate(common_name_full = paste0(common_name, " (", 
                                   str_to_lower(functional_type), ")")) %>%
  arrange(func_group, functional_type, common_name)

spp_list <- unique(wkprops$common_name_full)

# 4 species in each multi-panel plot
n_plots <- ceiling(length(spp_list) / 4)

pt_alpha <- 0.3

for (i in 1:n_plots) {
  spps <- spp_list[(i * 4 - 3):(i * 4)]
  props4 <- filter(wkprops, common_name_full %in% spps)
  
  # Get smoothed spline for each species, phenophase
  for (k in 1:4) {
    props <- filter(props4, common_name_full == spps[k])
    pg_list <- unique(props$phenogroup)
    
    for (j in pg_list) {
      
      props_pg <- props %>% filter(phenogroup == j)
      pg_preds <- data.frame(spline(x = props_pg$wk_doy4, 
                                    y = props_pg$prop_yes,
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
               aes(x = wk_doy4, y = prop_yes,
                   size = n_obs, color = phenogroup_f), alpha = pt_alpha) +
    scale_size_continuous(range = c(0.5, 4)) +
    geom_line(data = pg_preds4, 
              aes(x = x, y = y01, color = phenogroup_f)) +
    scale_color_manual(values = color_vec2) +
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
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
  )
  # ggsave(paste0("output/phenophase-periods-spp-", lpp_short, "-", i, ".png"),
  #        get(paste0("phenophase_periods_p", i)),
  #        width = 6.5,
  #        height = 7.5,
  #        units = "in",
  #        dpi = 600)
}

  