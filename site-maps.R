library(ggmap)
library(tmaptools)

xextend <- 1.0 # 0.3 for AMC, 1 for EWA
yextend <- 0.3
bb <- ggmap::make_bbox(lon = longitude, lat = latitude, data = sites, 
                       f = c(xextend, yextend))

# Map type, best options I identified below:
  # alidade_smooth (grayscale)
  # outdoors (nice with color = "bw")
  # stamen_terrain
maptype <- "outdoors"
# Map color (bw or color)
mapcol <- "bw"

# Note: Need to have an API key for Stadia maps first
# Convenient to save this to the .Renviron file
map <- ggmap::get_stadiamap(bbox = bb, 
                            maptype = maptype, 
                            color = mapcol,
                            force = TRUE, 
                            zoom = 13)
# Note: can save a ggmap object to a RData file and then read it back with load
# save(map, file = "xxx.RData")
# load(file = "xxx.RData")
#TODO: Figure out how to automate the zoom argument based on the scale of 
# observations. zoom = 7 worked for AMC; zoom 13 worked for EWA.

sites_map <- ggmap(map) +
  geom_point(data = sites, aes(x = longitude, y = latitude),
             color = "blue", alpha = 0.6, size = 1.5) +
  theme_void() # gets rid of axes completely
sites_map
# If we want to save map object without any surrounding white space by 
# specifying either the height or width of an object, we need to know the ratio
# (width to height). Using tmaptools to do that
asp <- tmaptools::get_asp_ratio(bb)

width <- 5
ggsave(paste0("output/map-", lpp_short, ".png"),
       sites_map,
       width = width,
       height = width / asp,
       units = "in", 
       dpi = 600)


#- Playing around with tmap package below -------------------------------------#

# extent <- 0.10 # x * 100 % of range
# lat_range <- max(sites$latitude) - min(sites$latitude)
# lat_min <- min(sites$latitude) - extent * lat_range
# lat_max <- max(sites$latitude) + extent * lat_range
# lng_range <- max(sites$longitude) - min(sites$longitude)
# lng_min <- min(sites$longitude) - extent * lng_range
# lng_max <- max(sites$longitude) + extent * lng_range
# center_lat <- mean(range(sites$latitude))
# center_lng <- mean(range(sites$longitude))
# bbox <- c(lng_min, lat_min, lng_max, lat_max)

library(tmap)
library(tmaptools)

coords <- st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326)
osm <- read_osm(coords, type = "osm", ext = 1.3)
background <- tm_shape(osm) + 
  tm_rgb()
amc_map <- background +
  tm_shape(coords) + 
  tm_dots(size = 0.4, 
          fill = "blue",
          fill_alpha = 0.5)
amc_map
tmap_save(amc_map, 
          filename = "C:/Users/erin/Desktop/test-map-esri.png",
          width = 5,
          units = "in",
          dpi = 600)
# This looks okay, but the map text is really small. Wondering whether this is
# something that can be adjusted, or whether it's best to use an interactive 
# map (via tmap or leaflet) and then print or save an image when zoomed 
# appropriately. 

tm <- qtm(coords)
map <- tmap_leaflet(tm)
map
# Interactive

tmap_mode("view")
tm_shape(coords) + tm_dots(size = 0.4, 
                           fill = "blue",
                           fill_alpha = 0.5)
tmap_mode("plot")

