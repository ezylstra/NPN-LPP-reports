require(daymetr)

# Create filenames ------------------------------------------------------------#

daymet_csv <- paste0("climate-data/daymet-", lpp_short,  
                     "-thru", max(all_yrs), ".csv")
daymet_zip <- str_replace(daymet_csv, ".csv", ".zip")

# Explore whether it's worth filtering sites (max of one per daymet cell) -----#

# Find daymet tile(s) that contain LPP points
tiles <- vect(tile_outlines)
sitesll <- sites %>%
  rename(lon = longitude, lat = latitude) %>%
  vect(geom = c("lon", "lat"), crs = "epsg:4326")
tilesllp <- extract(tiles, sitesll)
sites$daymet_tile <- tilesllp$TileID
tilesllp <- unique(sites$daymet_tile)

# Download tiled data
download_daymet_tiles(tiles = tilesllp,
                      start = 2020, end = 2020, param = "tmin")
# Convert nc files to geotiff
nc2tif(tempdir())
newpath <- tempdir() %>%
  str_replace_all("\\\\", "/")
tilefiles <- paste0(newpath, "/tmin_2020_", tilesllp, ".tif")

r_list <- list()
for (i in 1:length(tilesllp)) {
  r_list[[i]] <- terra::rast(tilefiles[i])[[1]]
}
r_merge <- sprc(r_list)
rm(r_list)
r_merge <- terra::merge(r_merge)
r_merge <- terra::project(r_merge, "epsg:4326")
sites$daymet_cell <- extract(r_merge, sitesll, cells = TRUE)[, "cell"]
# For AMC, this cuts the amount of data we need to extract in half... (50 v 105)

# Is this worth the effort?

# PICK UP HERE #############

# Download data ---------------------------------------------------------------#

get_daymet <- function(i) {
  
  temp_lat <- sites[i, ] %>% pull(lat)
  temp_lon <- sites[i, ] %>% pull(lon)
  temp_site <- sites[i, ] %>% pull(site_id)
  
  temp_daymet <- download_daymet(
    lat = temp_lat,
    lon = temp_lon, 
    start = first(yrs),
    end = last(yrs)
  ) %>%
    #--- just get the data part ---#
    .$data %>% 
    #--- assign site_id so we can match with onset data ---#
    mutate(site_id = temp_site) %>%
    #--- get date from day of the year ---#
    mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"))
  
  return(temp_daymet)
}

daymet_df <- lapply(1:nsites, get_daymet) %>%
  bind_rows()
