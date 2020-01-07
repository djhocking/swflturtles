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
library(lintr) # code linting
library(raster) # raster handling (needed for relief)
library(viridis) # viridis color scale
library(cowplot) # stack ggplots


#----- Load data -----
turtles_df <- read_csv("analysis/data/raw_data/bauri_telemetry_data_Naples_Preserve_2019-2020.csv") # data recorded on GPS in WGS84

#----- data summary -----
str(turtles_df) # structure of the data
summary(turtles_df) # general summary
(n_turtles <- length(unique(turtles_df$id))) # number of turtles tracked

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
# convert to lat-lon for human readable scale
points_latlon <- spTransform(points_sp, CRS("+proj=longlat"))
mcp_latlon <- spTransform(turtles_mcp, CRS("+proj=longlat"))

# ggmap now requires registration with google - consider other options
basemap_turtles <- get_map(location = c(lon = mean(points_latlon@coords[ , 1]),
                                        lat = mean(points_latlon@coords[ , 2])),
                           maptype = "hybrid",
                           zoom = 18,
                           source = "google")

# make dataframe for use in ggplot/ggmap
turtles_sdf <- data.frame(id = as.character(points_latlon@data$id),
                          points_latlon@coords)

map_turtles1 <- ggmap(basemap_turtles) +
  geom_polygon(data = tidy(mcp_latlon),
               aes(long, lat, colour = id, fill = id),
               alpha = 0.3) + # alpha sets the transparency
  geom_point(data = turtles_sdf,
             aes(x = x, y = y, colour = id))  +
  theme(legend.position = c(0.9, 0.70)) +
  labs(x = "Longitude", y = "Latitude")

map_turtles1

ggsave(map_turtles1, file = "analysis/figures/tbauri_mcp_100.pdf")



# Maybe get a rasta manually and use sf and ggplot2 rather than ggmap in the future given new api restrictions by google. use library sf in place of sp now - translate above to sf if possible later


#----- annimate movements -----
# visualize straight-line movements between relocation points

