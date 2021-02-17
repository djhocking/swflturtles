library(readxl)
library(dplyr)

turtles <- read_xlsx(here::here("analysis", "data", "raw_data", "Napes Preserve Box Turtle Data for GIS.xlsx"), sheet = "Naples Preserve Telemetry ",
                     col_types = c("numeric", "text", "date", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text"),
                     trim_ws = TRUE)

turtles <- turtles %>%
  dplyr::rename(id = `SCUTE ID`,
                pit = `PIT (Last 4)`,
                date = Date,
                season = Season,
                lat = Lat,
                lon = Long,
                sex = SEX,
                habitat = `Habitat: Oak Rosemary Scrub (ORS), Pine Flatwoods (PFW), Prairie Meadow (PM), Former Wetland (FWL), Other (O)`,
                substrate = `Substrate: Leaf Litter (LL), Pine Litter (PL), Sand (S), Other (O)`,
                buried = `Buried? Y/N`,
                burrow = `In Burrow? Y/N`,
                depth = `Buried Depth (CM)`,
                veg_type = `Vegetation Species`) %>%
  dplyr::mutate(lat = stringr::str_trim(lat),
                lon = stringr::str_trim(lon)) %>%
  dplyr::mutate(lat = as.numeric(lat),
                lon = as.numeric(lon),
                sex = toupper(sex),
                buried = toupper(sex),
                burrow = toupper(burrow))

summary(turtles)

if(!dir.exists("analysis/data/derived_data")) dir.create("analysis/data/derived_data")
saveRDS(turtles, here::here("analysis", "data", "derived_data", "turtle_telemetry.rds"))
