require(daymetr)

#- Create filenames -----------------------------------------------------------#

daymet_csv <- paste0("climate-data/daymet-", lpp_short,  
                     "-thru", max(all_yrs), ".csv")
daymet_zip <- str_replace(daymet_csv, ".csv", ".zip")

#- Explore whether it's worth filtering sites (max of one per daymet cell) ----#
 
# # Find daymet tile(s) that contain LPP points
# tiles <- vect(tile_outlines)
# sitesll <- sites %>%
#   rename(lon = longitude, lat = latitude) %>%
#   vect(geom = c("lon", "lat"), crs = "epsg:4326")
# tilesllp <- extract(tiles, sitesll)
# sites$daymet_tile <- tilesllp$TileID
# tilesllp <- unique(sites$daymet_tile)
# 
# # Download tiled data
# download_daymet_tiles(tiles = tilesllp,
#                       start = 2020, end = 2020, param = "tmin")
# # Convert nc files to geotiff
# nc2tif(tempdir())
# newpath <- tempdir() %>%
#   str_replace_all("\\\\", "/")
# tilefiles <- paste0(newpath, "/tmin_2020_", tilesllp, ".tif")
# 
# r_list <- list()
# for (i in 1:length(tilesllp)) {
#   r_list[[i]] <- terra::rast(tilefiles[i])[[1]]
# }
# r_merge <- sprc(r_list)
# rm(r_list)
# r_merge <- terra::merge(r_merge)
# r_merge <- terra::project(r_merge, "epsg:4326")
# sites$daymet_cell <- extract(r_merge, sitesll, cells = TRUE)[, "cell"]

# For AMC, this cuts the amount of data we need to extract in half... (50 v 105)
# Is this worth the effort?

#- Download data --------------------------------------------------------------#

# Write location data to temporary csv file for batch download
locs <- sites %>%
  select(site_id, latitude, longitude) %>%
  rename(lat = latitude,
         lon = longitude)
write.csv(locs, paste0(tempdir(), "/locs.csv"), row.names = FALSE)
#TODO: See whether there are duplicate geographic locations? If so, could remove
# prior to download and attach to both sites later.

# Download data
# (One year prior to first year of observations through last year)
start <- min(yearst$yr) - 1
end <- max(yearst$yr)
clim_data <- download_daymet_batch(
  file_location = paste0(tempdir(), "/locs.csv"),
  start = start,
  end = end,
  internal = TRUE,
  simplify = TRUE
)

# Clean up data and put in wide format (one row per site-day)
clim_data <- clim_data %>%
  select(-c(tile, latitude, longitude, altitude)) %>%
  mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j")) %>%
  mutate(date = ymd(date)) %>%
  mutate(mon = month(date),
         day = day(date)) %>%
  select(site, date, year, mon, day, yday, measurement, value) %>%
  pivot_wider(names_from = measurement,
              values_from = value) %>%
  rename(daylength_s = dayl..s.,
         prcp_mm = prcp..mm.day.,
         rad_Wm2 = srad..W.m.2.,
         swe_kgm2 = swe..kg.m.2.,
         tmax_c = tmax..deg.c.,
         tmin_c = tmin..deg.c.,
         vp_Pa = vp..Pa.) %>%
  data.frame()

# Save daily data to zip file (need to save, then remove csv file)
write.csv(clim_data, file = daymet_csv, row.names = FALSE)
zip(daymet_zip, files = daymet_csv)
file.remove(daymet_csv)

# Remove dataframe
rm(clim_data)
 