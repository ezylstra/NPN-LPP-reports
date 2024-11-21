# Example LPP report
# Erin Zylstra
# ezylstra@usanpn.org
# 2024-11-21

require(rnpn)
require(dplyr)
require(lubridate)
require(stringr)
require(ggplot2)

rm(list = ls())

# Name and nickname of LPP
lpp <- "Appalachian Mountain Club"
lpp_short <- "amc"
# lpp <- "Earth"
# lpp_short <- 

# Specify years of interest
all_yrs <- 2007:2023

# Check if the status/intensity data have already been downloaded and saved. 
# If so, skip this step and load csv
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
    climate_data = FALSE,  # TODO: see whether this makes sense or whether I should be downloaded data separately
    additional_fields = "site_name"
  )
  
  # Remove some fields that are irrelevant (to make file smaller)
  si <- si %>%
    select(-c(update_datetime, species_id, species, phenophase_id,
              intensity_category_id))
  
  write.csv(si, si_csv_name, row.names = FALSE)
} 

# Load status/intensity data
si <- read.csv(si_csv_name)

# Format date
si <- si %>%
  mutate(obsdate = ymd(observation_date)) %>%
  mutate(yr = year(observation_date)) %>%
  select(-observation_date)

# Quick glance at data dimensions for AMC
dim(si) #529,642
count(si, site_name); sort(count(si, site_name)$n) # 105 sites (5 with < 50 records)
count(si, yr) # 2012-2023, but only 7 and 149 records in 2012, 2013 
count(filter(si, kingdom == "Plantae"), common_name)   # 21 plant spp
count(filter(si, kingdom == "Animalia"), common_name)  # 7 animal spp (but monarch only 1 record)



