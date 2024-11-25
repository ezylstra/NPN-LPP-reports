# Example LPP report
# Erin Zylstra
# erinz@usanpn.org
# 2024-11-21

require(rnpn)
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

# Append phenophase "groups"
pheno_list <- read.csv("phenophases.csv")
si <- left_join(si, select(pheno_list, phenophase_id, pheno_group),
                by = "phenophase_id")

# Append intensity midpoints (or single values)
intensity_list <- read.csv("intensities.csv")
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
# Any phenophase_status == -1 when intensity/abundance value present?
count(si, phenophase_status, is.na(intensity_value), is.na(abundance_value))
# Apparently, this does happen.
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

# Check that phenophase status is 1 if there is an intenstity/abundance value
count(si, phenophase_status, is.na(intensity), is.na(abundance_value))

#- Dealing with multiple observations by the same person ----------------------#
# Same plant or animal (same site), same phenophase, same date - same observer

#TODO: Figure out if/where we want QA/QC (within indiv-date-observer) to happen

# For now, will keep record with more advanced phenophase or higher intensity or
# abundance value. Will do this by sorting observations and keeping only the 
# first
idop <- si %>%
  group_by(common_name, individual_id, obsdate, person_id, phenophase_id) %>%
  summarize(n_obs = n(),
            .groups = "keep") %>%
  data.frame()
idop$obsnum <- 1:nrow(idop)

si <- si %>%
  arrange(person_id, individual_id, obsdate, phenophase_id, 
          desc(phenophase_status), desc(intensity), desc(abundance_value)) %>%
  left_join(select(idop, -c(n_obs, common_name)), 
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

#- Dealing with multiple observations by different people ---------------------#
# Same plant or animal, same phenophase, same date - different observers

# For plants, multiple observations of the same individual and phenophase on the 
# same date can tell us about inter-observer variation (assuming phenophase 
# status and intensity doesn't change).
# For animals, multiple observations of the same species at a site on the same
# date provides a confounded measure of inter-observer variation and variation 
# in abundance or activity of the species.

# To summarize sampling intensity (when and where observations were made), 
# probably best to keep all the observations in the dataset. 

# To analyze phenophase onsets or animal activity/occurrence, probably best
# to subset observations so there is a maximum of one record for each 
# phenophase associated with a plant/animal at a given site each day.
ido <- si %>%
  group_by(common_name, individual_id, obsdate, phenophase_id) %>%
  summarize(n_obs = n(),
            .groups = "keep") %>%
  data.frame()
ido$obsnum <- 1:nrow(ido)

# For an animal species, multiple observations doesn't seem too problematic 
# because a status = 0 could easily be a detection issue. If we assume there are
# no false positives, then we can simply use the maximum status value and 
# maximum intensity or abundance values.

# For an individual plant, we have some choices about how to characterize
# the phenological state and intensity (if recorded) when multiple observations
# were made: 
# 1) retain one observation of the plant per day.
# 2) retain all data, use status = 1 observations, and for intensity, average 
#    over values

si_sub_animals <- si %>%
  filter(kingdom == "Animalia") %>%
  arrange(common_name, individual_id, obsdate, phenophase_id, 
          desc(phenophase_status), desc(abundance_value), desc(intensity)) %>%
  left_join(select(ido, -c(n_obs, common_name)),
            by = c("individual_id", "obsdate", "phenophase_id")) %>%
  mutate(dups = sequence(rle(as.character(obsnum))$lengths)) %>%
  filter(dups == 1) %>%
  select(-c(obsnum, dups)) %>%
  arrange(common_name, obsdate, person_id, phenophase_id)
  

#- PICK UP HERE ---------------------------------------------------------------#

#TODO: After removing multiple observations of the same phenophase, we'll 
# aggregate information within pheno group 



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

