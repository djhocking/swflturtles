#############################
# Analysis of Southwestern Florida Box Turtle Radiotelemetry Data
# Daniel J. Hocking
# Collaborators: Jordan Donini
#############################

#----- Load packages -----
library(tidyverse)
library(adehabitatHR) # for mcp function
library(scales)
suppressPackageStartupMessages(library(ggmap))
library(broom)
library(ggsn)
library(raster)
# library(lintr) # code linting
# library(raster) # raster handling (needed for relief)
# library(viridis) # viridis color scale
# library(cowplot) # stack ggplots


#----- Load data -----
turtles_df <- read_csv("analysis/data/raw_data/bauri_telemetry_data_Naples_Preserve_2019-2020.csv") # data recorded on GPS in WGS84
individuals_df <- read_csv("analysis/data/raw_data/bauri_individual_info.csv") %>% # individual covariates including sex
  mutate(id = as.character(id))

#----- data summary -----
str(turtles_df) # structure of the data
summary(turtles_df) # general summary
(n_turtles <- length(unique(turtles_df$id))) # number of turtles tracked
(n_turtles <- length(unique(turtles_df$pit)))

# number of observations (non-unique locations) per turtle
n_obs <- turtles_df %>%
  group_by(id) %>%
  summarise(obs = n())
n_obs

# proportion of locations in each habitat

# proportion of each behavior

#------ prep data spatial ------
# currently using the sp package rather than sf to work with adehabitatHR
points_sp <- turtles_df %>%
  dplyr::select(id, x = lon, y = lat) %>%
  mutate(x = -1 * x) %>%
  as.data.frame()

coordinates(points_sp) <- c("x", "y")

# set coordinate system and projection
proj4string(points_sp) <- CRS( "+proj=longlat +datum=WGS84") # +units=m +no_defs" )

points_utm <- spTransform(points_sp, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m"))

#----- Calculate MCPs -----
turtles_mcp <- mcp(points_utm, percent = 100, unin = "m", unout = "ha")

# View results
turtles_mcp
summary(turtles_mcp)

# Plot the home range size by percent MCP
mcp_ranges <- mcp.area(points_utm, percent = seq(50, 100, by = 5),
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


# color density for overlap density to see hot spot areas?


#----- Kernal density estimates and maps -----
library(sf)
library(rgeos)
library(leaflet)

# bivariate normal utilization distribtution for prob density of finding animal at a point

# base kernel with default ad hoc smoothing factor
turtle_kernel <- kernelUD(points_utm, h = "href", )  # href = the reference bandwidth

# kernel using Least Square Cross Validation for smoother
turtle_kernel_lscv <- kernelUD(points_utm, h = "LSCV")

# plotLSCV(turtle_kernel_lscv)

class(turtle_kernel)

# convert kernel to spatial polygon
kernel_sp <- getverticeshr(turtle_kernel_lscv, percent = 90, unin = "m", unout = "ha")
data.frame(kernel_sp)

home_ranges <- turtles_mcp %>%
  data.frame() %>%
  rename(mcp_100 = area) %>%
  left_join(data.frame(kernel_sp)) %>%
    rename(kernel_90 = area)
write_csv(home_ranges, "analysis/data/derived_data/home_ranges.csv")

kernel_latlon <- spTransform(kernel_sp, CRS("+proj=longlat"))

# make dataframe for use in ggplot/ggmap
polys_kernel <- as.data.frame(broom::tidy(kernel_latlon)) %>%
  mutate(id = as.character(id)) %>%
  left_join(individuals_df)

map_turtles_kernels <- ggmap(basemap_turtles, extent = "panel", maprange = FALSE) +
  geom_polygon(data = polys_kernel,
               aes(x = long, y = lat, fill = id, colour = id),
               alpha = 0.3)  # +
  # coord_cartesian(xlim = c(min(points_latlon@coords[ , 1])-0.0002, max(points_latlon@coords[ , 1])+0.0002),
  #                 ylim = c(min(points_latlon@coords[ , 2])-0.0002, max(points_latlon@coords[ , 2])+0.0002)) +
  # theme(legend.position = c(0.94, 0.72)) +
  # labs(x = "Longitude", y = "Latitude") +
  # ggsn::scalebar(polys, dist = 25, st.size=3, height=0.01, dist_unit = "m", transform = TRUE, model = "WGS84")
map_turtles_kernels

ggsave(map_turtles2, file = "analysis/figures/tbauri_kernel_90.pdf", width = 8, units = "in")


kernel_vol <- getvolumeUD(turtle_kernel_lscv)

image(kernel_vol[[1]])
# plot(kernel_vol[[2]], add = T)
xyzv <- as.image.SpatialGridDataFrame(kernel_vol[[1]])
contour(xyzv, add=TRUE)

image(kernel_vol[[2]])
contour(as.image.SpatialGridDataFrame(kernel_vol[[2]]), add=TRUE)

turtle_kernel_lscv2 <- kernelUD(points_utm, h = "LSCV", same4all = TRUE)
# kernel_vol2 <- getvolumeUD(turtle_kernel_lscv, same4all = TRUE)
foo <- estUDm2spixdf(turtle_kernel_lscv2)
plot(foo) # potential to separate them or add them all on the basemap?
# image(foo)

crs(kernel_vol[[1]]) <- "+proj=utm +zone=17 +datum=WGS84 +units=m"
kernel_vol[[1]] <- projectRaster(raster(kernel_vol[[1]]), crs = "+proj=longlat")
kernal_rasters <- list()
# kernal_rasters[[1]] <- as.data.frame(broom::tidy(kernel_vol[[1]]))


kernal_rasters[[1]] <- as.data.frame.estUD(kernel_vol[[1]]) %>%
  rename(lon = x, lat = y)
# projectRaster(kernal_rasters[[1]], crs = "+proj=longlat")
# crs(kernal_rasters[[1]]) <- "+proj=utm +zone=17 +datum=WGS84 +units=m"

rtp <- rasterToPolygons(kernel_vol[[1]])

bm <- ggmap(basemap_turtles, extent = "panel")
bm +
  geom_polygon(data = rtp,
               aes(x = long, y = lat,
                   fill = rep(rtp$n, each = 5)),
               size = 0,
               alpha = 0.5)  +
  scale_fill_gradientn("RasterValues", colors = topo.colors(255))

ggmap(basemap_turtles, extent = "panel") +
  inset_raster(raster(kernel_vol[[1]]))

ggmap(basemap_turtles, extent = "panel") +
  inset_raster(raster(xyzv))



ggmap(basemap_turtles, extent = "panel") +
  geom_raster(kernal_rasters[[1]], aes(lon, lat))


ggmap(basemap_turtles, extent = "panel") +
  # ggplot(data = xyzv, aes(x, y)) +
  geom_tile(xyzv, aes(x, y, fill = z))


ggplot(kernel_vol[[1]], aes(lon, lat, fill = n))


kernel_latlon %>%
  st_as_sf() %>%
  leaflet() %>%
  addTiles() %>%
  addPolygons(fillColor = id, group = id)

kernel_latlon %>%
  st_as_sf() %>%
  leaflet() %>%
  addTiles() %>%
  # addRasterImage(raster(kernal_rasters[[1]]))
  addRasterImage(raster(kernel_vol[[1]]), opacity = 0.5)
# addRasterImage(xyzv)

str(raster(kernel_vol[[1]]))
summary(raster(kernel_vol[[1]])@data@values)
hist(raster(kernel_vol[[1]])@data@values)


foo <- raster(kernel_vol[[1]]) # filter to values <90?

kernel_latlon %>%
  st_as_sf() %>%
  leaflet() %>%
  addTiles() %>%
  # addRasterImage(raster(kernal_rasters[[1]]))
  addRasterImage(raster(kernel_vol[[1]]), opacity = 0.5)



#----- annimate movements -----
# visualize straight-line movements between relocation points




#----- clean up -----

rm(list = ls())
