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
si <- si %>% distinct(across(-observation_id), .keep_all = TRUE)
#TODO: Probably don't need to keep observation_id with .keep_all arg

# Was going to remove records with unknown phenophase_status (-1) here, but
# I'll wait since it's possible some values could be inferred from reports of 
# other phenophases.

# Append phenophase "groups"
pheno_list <- read.csv("phenophases.csv")
si <- left_join(si, select(pheno_list, phenophase_id, pheno_group),
                by = "phenophase_id")

# Append intensity midpoints (or single values)


# Append abundance midpoints (or single values)



#- Quick overviews of data ----------------------------------------------------#
# Data dimensions
dim(si) # 528,810
count(si, site_name); sort(count(si, site_name)$n) # 105 sites (5 with < 50 records)
count(si, yr) # 2012-2023, but only 7 and 149 records in 2012, 2013 
count(filter(si, kingdom == "Plantae"), common_name)   # 21 plant spp ('ohi'a lehua with 28, rest with > 1900)
count(filter(si, kingdom == "Animalia"), common_name)  # 7 animal spp (but monarch only 1 record)
count(si, phenophase_description) # 41

# How many phenophases were recorded for each species?
spp_ph <- si %>%
  group_by(kingdom, common_name) %>%
  summarize(n_phenophases = length(unique(phenophase_description)),
            .groups = "keep") %>%
  data.frame()
spp_ph # Except for monarch, 6-13 phenophases per spp
# How many species for each phenophase?
ph_spp <- si %>%
  group_by(kingdom, phenophase_description) %>%
  summarize(n_spp = length(unique(common_name)),
            .groups = "keep") %>%
  data.frame()
ph_spp # 1-21 species per phenophase

#- Where possible, impute phenophase status -----------------------------------#
# Any phenophase_status == -1 when intensity/abundance value present?
count(si, phenophase_status, is.na(intensity_value), is.na(abundance_value))
# Apparently, this does happen.
#TODO: QA/QC fix. Can this be prevented in the app?

# Change unknown phenophase_status to 1 if an intensity/abundance value reported
si <- si %>%
  mutate(phenophase_status = case_when(
    !is.na(intensity_value) ~ 1,
    !is.na(abundance_value) ~ 1,
    .default = phenophase_status
  ))



#- PICK UP HERE ---------------------------------------------------------------#
# Want to summarize the number of observations. Here, defining an observation as 
# all the data collected on a plant or an animal species by an observer on a 
# particular date. 

# Summarize information associated with each observation
obs <- si %>%
  group_by(kingdom, common_name, individual_id, obsdate, yr, person_id) %>%
  summarize(n = n(),
            n_ph = length(unique(phenophase_description)),
            n_ph_unk = sum(phenophase_status == -1),
            .groups = "keep") %>%
  data.frame()

count(obs, n, n_ph, n_ph_unk, name = "count")
count(obs, kingdom, n == n_ph, n_ph_unk > 0, name = "count") 
# 50 times an observer monitored a plant multiple times in one day
# Check an example where this occurred when phenophase status was always known
head(filter(obs, n > n_ph, n_ph_unk == 0))
si %>%
  filter(common_name == "Bigelow's sedge",
         individual_id == 59773,
         obsdate == "2018-06-16", 
         person_id == 42888) %>%
  select(observation_id, common_name, phenophase_description, 
         phenophase_status, intensity_value) %>%
  arrange(phenophase_description)
# In this case, reports of flower heads with different status & reports of 
# open flowers & fruit with different intensity values
#TODO: Figure out if/where we want QA/QC (within observer-individual-date) to happen


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

