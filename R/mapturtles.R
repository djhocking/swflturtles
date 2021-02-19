#############################
# Analysis of Southwestern Florida Box Turtle Radiotelemetry Data
# Daniel J. Hocking
# Collaborators: Jordan Donini
#############################

#----- Set checkpoint -----

library(checkpoint)
checkpoint("2021-02-15")

#----- Load packages -----
library(tidyverse)
library(adehabitatHR) # for mcp function
library(scales)
library(ggplot2)
suppressPackageStartupMessages(library(ggmap))
library(broom)
library(ggsn)
library(raster)
library(lubridate)
# library(lintr) # code linting
# library(raster) # raster handling (needed for relief)
# library(viridis) # viridis color scale
# library(cowplot) # stack ggplots

if(!dir.exists("analysis/figures")) dir.create("analysis/figures")

#----- Custom functions (move and then source) ------

# consider removing lat and lon axis labels to protect the turtles and because unnecessary
map_kernel <- function(basemap, raster_df, xmin = -81.79898, xmax = -81.79702, ymin = 26.16420, ymax = 26.16769) {
  map_k90 <- ggmap(basemap, extent = "normal") +
    geom_raster(data = raster_df, aes(x = lon, y = lat, fill = Density), alpha = 0.5) +
    coord_cartesian(xlim = c(xmin, xmax),
                    ylim = c(ymin, ymax)) +
    # coord_fixed() +
    scale_fill_viridis(option = "magma") +
    facet_wrap(vars(season)) + # make wrap variable an option
    # theme(legend.position = "bottom") + # c(0.90, 0.72)) +
    labs(x = "Longitude", y = "Latitude") +
    ggsn::scalebar(raster_df, dist = 25, st.size = 4, height=0.01, dist_unit = "m", transform = TRUE, model = "WGS84", st.color = "white", location = "bottomright") +
    theme_bw() +
    theme(legend.justification = c(1, 1),
          legend.position = c(1, 1),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("Individual: ", unique(raster_df$id), " (", unique(raster_df$sex), ")"))

  return(map_k90)
}

#----- Load data -----
# turtles_df <- readRDS("analysis/data/derived_data/turtle_telemetry.rds") # data recorded on GPS in WGS84
turtles_df <- read_csv("analysis/data/raw_data/bauri_telemetry_data_Naples_Preserve_2019-2020.csv")
individuals_df <- read_csv("analysis/data/raw_data/bauri_individual_info.csv") %>% # individual covariates including sex
  mutate(id = as.character(id))

turtles_df %>%
  dplyr::filter(is.na(id))

#----- data summary -----
str(turtles_df) # structure of the data
summary(turtles_df) # general summary
(n_turtles <- length(unique(turtles_df$id))) # number of turtles tracked
(n_turtles <- length(unique(turtles_df$pit)))

# Remove missing points
turtles_df <- turtles_df %>%
  dplyr::filter(!is.na(lat),
                !is.na(lon))

# number of observations (non-unique locations) per turtle
n_obs <- turtles_df %>%
  dplyr::mutate(id = paste0(id, "_", pit)) %>%
  group_by(id) %>%
  summarise(obs = n())
n_obs

n_obs <- turtles_df %>%
  group_by(id) %>%
  summarise(obs = n())
n_obs

# proportion of locations in each habitat

# proportion of each behavior

#------ prep data spatial ------
# currently using the sp package rather than sf to work with adehabitatHR
points_sp <- turtles_df %>%
  dplyr::select(id, x = lon, y = lat, season) %>%
  mutate(x = -1 * x) %>%
  as.data.frame()

coordinates(points_sp) <- c("x", "y")

# set coordinate system and projection
proj4string(points_sp) <- CRS( "+proj=longlat +datum=WGS84") # +units=m +no_defs" )

points_utm <- spTransform(points_sp, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m"))

points_utm1 <- points_utm
points_utm1 <- points_utm[ , "id"]

#----- Calculate MCPs -----
turtles_mcp <- mcp(points_utm1, percent = 100, unin = "m", unout = "ha")
turtles_mcp_100 <- turtles_mcp

# View results
turtles_mcp
summary(turtles_mcp)

# Plot the home range size by percent MCP
mcp_ranges <- mcp.area(points_utm1, percent = seq(50, 100, by = 5),
                       unin = c("m"),
                       unout = c("ha"), plotit = TRUE)

# Home ranges can be tricky, especially with turtles and with few relocations or movements. Consider if any migrations (e.g. female breeding, season switching for temperature, moisture, reproductive opportunity, or food resources) or dispersal (e.g. permanent movement to new home range area, often juveniles but not always) occur.

# basic plot to make sure things look ok
plot(points_utm, col = as.factor(points_utm@data$id), pch = 20)
plot(turtles_mcp, col = alpha(1:n_turtles, 0.5), add = TRUE)

#---- plot home ranges on base map -----
# convert to lat-lon for human readable scale for map
points_latlon <- spTransform(points_sp, CRS("+proj=longlat"))
mcp_latlon <- spTransform(turtles_mcp, CRS("+proj=longlat"))

# ggmap now requires registration with google - consider other options
basemap_turtles <- get_map(location = c(lon = mean(points_latlon@coords[ , 1]),
                                        lat = mean(points_latlon@coords[ , 2])),
                           maptype = "satellite",
                           zoom = 17,
                           source = "google")

# make dataframe for use in ggplot/ggmap
turtles_sdf <- data.frame(id = as.character(points_latlon@data$id),
                          points_latlon@coords) %>%
  left_join(individuals_df)

polys <- as.data.frame(broom::tidy(mcp_latlon)) %>%
  mutate(id = as.character(id)) %>%
  left_join(individuals_df)

map_turtles1 <- ggmap(basemap_turtles, extent = "panel") +
  geom_polygon(data = polys,
               aes(x = long, y = lat, fill = id, colour = id),
               alpha = 0.3) +
  # geom_point(data = turtles_sdf,
  #            aes(x = x, y = y, colour = id))  +
  coord_cartesian(xlim = c(min(points_latlon@coords[ , 1])-0.0002, max(points_latlon@coords[ , 1])+0.0002),
                  ylim = c(min(points_latlon@coords[ , 2])-0.0002, max(points_latlon@coords[ , 2])+0.0002)) +
  theme(legend.position = c(0.94, 0.72)) +
  labs(x = "Longitude", y = "Latitude") +
  ggsn::scalebar(polys, dist = 25, st.size=3, height=0.01, dist_unit = "m", transform = TRUE, model = "WGS84")

# use consistent colors by turtle so can compare across maps more easily
# scale_fill_manual(name = "Turtle ID",
#                   values = num,
#                   breaks = id) +
#   scale_colour_manual(name = "Turtle ID",
#                       values = num,
#                       breaks = id) +

# label the mis-marked point
# ggplot(turtles_df, aes(x = -1*lon, y = lat)) + geom_point() + geom_text(data = filter(turtles_df, id == 3), aes(label = date), hjust = 0, vjust = 0)

ggsave(map_turtles1, file = "analysis/figures/tbauri_mcp_100.pdf", width = 4, units = "in")
ggsave(map_turtles1, file = "analysis/figures/tbauri_mcp_100.png")
ggsave(map_turtles1, file = "analysis/figures/tbauri_mcp_100.tiff", dpi = 1000, width = 4, units = "in")

# Maybe get a raster manually and use sf and ggplot2 rather than ggmap in the future given new api restrictions by google. use library sf in place of sp now - translate above to sf if possible later

# separate for males and females

map_turtles2 <- map_turtles1 +
  geom_point(data = turtles_sdf,
             aes(x = x, y = y, colour = id))  +
  facet_wrap(~sex)


ggsave(map_turtles2, file = "analysis/figures/tbauri_mcp_100_sex.pdf", width = 8, units = "in")
ggsave(map_turtles2, file = "analysis/figures/tbauri_mcp_100_sex.png")
ggsave(map_turtles2, file = "analysis/figures/tbauri_mcp_100_sex.tiff", dpi = 1000, width = 8, units = "in")

#----------- MCP by season and sex ------------

# trick the stupid mcp function by making id = id_season, then break back apart for ggplot
points_sp <- turtles_df %>%
  dplyr::mutate(id = paste0(id, "_", season)) %>%
  dplyr::select(id, x = lon, y = lat, season) %>%
  mutate(x = -1 * x) %>%
  as.data.frame()

coordinates(points_sp) <- c("x", "y")
proj4string(points_sp) <- CRS( "+proj=longlat +datum=WGS84") # +units=m +no_defs" )
points_utm <- spTransform(points_sp, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m"))

points_utm1 <- points_utm
points_utm1 <- points_utm[ , "id"]

#----- Calculate MCPs -----
turtles_mcp <- mcp(points_utm1, percent = 100, unin = "m", unout = "ha")
turtles_mcp_sex_season <- turtles_mcp

# View results
turtles_mcp
summary(turtles_mcp)

# Plot the home range size by percent MCP
mcp_ranges <- mcp.area(points_utm1, percent = seq(50, 100, by = 5),
                       unin = c("m"),
                       unout = c("ha"), plotit = FALSE)

points_latlon <- spTransform(points_sp, CRS("+proj=longlat"))
mcp_latlon <- spTransform(turtles_mcp, CRS("+proj=longlat"))

# ggmap now requires registration with google - consider other options
# basemap_turtles <- get_map(location = c(lon = mean(points_latlon@coords[ , 1]),
#                                         lat = mean(points_latlon@coords[ , 2])),
#                            maptype = "satellite",
#                            zoom = 17,
#                            source = "google")

# make dataframe for use in ggplot/ggmap
turtles_sdf <- data.frame(id = as.character(points_latlon@data$id),
                          points_latlon@coords) %>%
  dplyr::mutate(id2 = id) %>%
  tidyr::separate(id, into = c("id", "season"), sep = "_") %>%
  left_join(individuals_df)

polys <- as.data.frame(broom::tidy(mcp_latlon)) %>%
  mutate(id2 = id) %>%
  tidyr::separate(id, into = c("id", "season"), sep = "_") %>%
  mutate(id = as.character(id)) %>%
  left_join(individuals_df)

map_turtles1 <- ggmap(basemap_turtles, extent = "panel") +
  geom_polygon(data = polys,
               aes(x = long, y = lat, fill = id, colour = id),
               alpha = 0.2) +
  # geom_point(data = turtles_sdf,
  #            aes(x = x, y = y, colour = id))  +
  coord_cartesian(xlim = c(min(points_latlon@coords[ , 1])-0.0002, max(points_latlon@coords[ , 1])+0.0002),
                  ylim = c(min(points_latlon@coords[ , 2])-0.0002, max(points_latlon@coords[ , 2])+0.0002)) +
  # facet_grid(rows = season, cols = sex) +
  theme(legend.position = "bottom") + # c(0.90, 0.72)) +
  labs(x = "Longitude", y = "Latitude") +
  ggsn::scalebar(polys, dist = 25, st.size=3, height=0.01, dist_unit = "m", transform = TRUE, model = "WGS84")

# use consistent colors by turtle so can compare across maps more easily
# scale_fill_manual(name = "Turtle ID",
#                   values = num,
#                   breaks = id) +
#   scale_colour_manual(name = "Turtle ID",
#                       values = num,
#                       breaks = id) +

# separate for males and females

map_turtles2 <- map_turtles1 +
  geom_point(data = turtles_sdf,
             aes(x = x, y = y, colour = id), alpha = 0.8, size = 0.5)  +
  facet_grid(rows = vars(season), cols = vars(sex))

saveRDS(map_turtles2, file = "analysis/figures/map_turtles_sex_season.rds")

ggsave(map_turtles2, file = "analysis/figures/tbauri_mcp_100_sex_season.pdf", width = 8, units = "in")
ggsave(map_turtles2, file = "analysis/figures/tbauri_mcp_100_sex_season.png")
ggsave(map_turtles2, file = "analysis/figures/tbauri_mcp_100_sex_season.tiff", dpi = 1000, width = 8, units = "in")


#----------- MCP by season and sex -----------


#----------- MCP by year and season ----------


# color density for overlap density to see hot spot areas?

#------------- split by year and sex--------------

# need to do mcps by those groups then combine into one

#----- Kernal density estimates and maps -----

############ problem with kernals is doesn't limit by fence
library(sf)
library(rgeos)
library(leaflet)

# bivariate normal utilization distribtution for prob density of finding animal at a point

# base kernel with default ad hoc smoothing factor
turtle_kernel <- kernelUD(points_utm1, h = "href")  # href = the reference bandwidth

# kernel using Least Square Cross Validation for smoother
turtle_kernel_lscv <- kernelUD(points_utm1, h = "LSCV")

# plotLSCV(turtle_kernel_lscv)

class(turtle_kernel)

# convert kernel to spatial polygon
kernel_sp <- getverticeshr(turtle_kernel_lscv, percent = 90, unin = "m", unout = "ha")
data.frame(kernel_sp)

home_ranges <- distinct(points_utm@data) %>%
  left_join(data.frame(turtles_mcp)) %>%
  rename(mcp_100 = area) %>%
  left_join(data.frame(kernel_sp)) %>%
    rename(kernel_90 = area) %>%
  dplyr::mutate(id = stringr::str_extract(id, "[^_]+")) %>%
  left_join(individuals_df) %>%
  dplyr::arrange(season, sex, num, id)

write_csv(home_ranges, "analysis/data/derived_data/home_ranges.csv")


#---------- Plot Kernel Densities -------------

kernel_vol <- getvolumeUD(turtle_kernel_lscv)

names(kernel_vol)
raster_dry_1 <- raster(as(kernel_vol$`1_Dry`,"SpatialPixelsDataFrame")) # convert to raster

raster_dry_1_latlon <- raster::projectRaster(raster_dry_1, crs = "+proj=longlat +datum=WGS84") # reproject raster to latlon
raster_dry_1_pts_latlon <- raster::rasterToPoints(raster_dry_1_latlon, spatial = TRUE) # convert to points

# convert to dataframe and filter to 90% kernel density max
raster_dry_1_df <- as.data.frame(raster_dry_1_pts_latlon, xy = TRUE) %>%
  dplyr::filter(n <= 90) %>%
  dplyr::select(lon = x, lat = y, Density = n)

for(i in 1:length(names(kernel_vol))) {

  raster_dry_1 <- raster(as(kernel_vol[[names(kernel_vol)[i]]],"SpatialPixelsDataFrame")) # convert to raster

  raster_dry_1_latlon <- raster::projectRaster(raster_dry_1, crs = "+proj=longlat +datum=WGS84") # reproject raster to latlon
  raster_dry_1_pts_latlon <- raster::rasterToPoints(raster_dry_1_latlon, spatial = TRUE) # convert to points

  # convert to dataframe and filter to 90% kernel density max
  tmp <- as.data.frame(raster_dry_1_pts_latlon, xy = TRUE) %>%
    dplyr::filter(n <= 90) %>%
    dplyr::select(lon = x, lat = y, Density = n) %>%
    dplyr::mutate(grp = names(kernel_vol)[i]) %>%
    tidyr::separate(grp, into = c("id", "season"), sep = "_") %>%
    mutate(id = as.character(id)) %>%
    left_join(individuals_df)

  if(i == 1) {
    raster_dry_1_df <- tmp
  } else {
    raster_dry_1_df <- dplyr::bind_rows(raster_dry_1_df, tmp)
  }

}


for(i in 1:length(unique(raster_dry_1_df$id))) {
  df1 <- raster_dry_1_df %>%
    dplyr::filter(id == unique(raster_dry_1_df$id)[i])
  m2 <- map_kernel(basemap = basemap_turtles, raster_df = df1)

  ggsave(filename = paste0("analysis/figures/kernel_map_season_", unique(df1$id), ".tiff"), plot = m2, height = 10, width = 14, units = "in", dpi = 100)
}




df1 <- raster_dry_1_df %>%
  dplyr::filter(id == 1)
m1 <- map_kernel(basemap = basemap_turtles, raster_df = df1)



#----- annimate movements -----
# visualize straight-line movements between relocation points




#----- clean up -----

rm(list = ls())
gc()
