library(terra)
library(dplyr)
library(rnaturalearth)
library(stringr)


# set.seed(0)

df <- read.csv(
  "data/vgp_iucn_loc_metadata.csv",
  stringsAsFactors = FALSE
)

# species_list <- df %>%
#   filter(lineage == "Invertebrates") %>%
#   distinct(iucn_name) %>%
#   pull(iucn_name)
# 
# writeLines(species_list, "/global/scratch/users/minoli/maps/data/invertebrates_list.txt")
# 
# shape_dir <- "/global/scratch/users/minoli/maps/data/shapes"
# gpkg_out  <- file.path(shape_dir, "fishes_combined.gpkg")
# 
# shps <- list.files(
#   path       = shape_dir,
#   #pattern    = "REPTILES_PART[12]\\.shp$",
#   #pattern    = "MAMMALS_PART[12]\\.shp$",
#   #pattern    = "AMPHIBIANS_PART[12]\\.shp$",
#   pattern    = "^(MARINEFISH_PART[0-9]+|CROAKERS_DRUMS|EELS|GROUPERS|HAGFISH|SALMONIDS|SEABREAMS_SNAPPERS_GRUNTS|SHARKS_RAYS_CHIMAERAS|STURGEONS_PADDLEFISHES|SYNGNATHIFORM_FISHES|TUNAS_BILLFISHES_SWORDFISH|WRASSES_PARROTFISHES)\\.shp$",
#   full.names = TRUE
# )
# 
# v <- terra::vect(shps)    
# terra::writeVector(v, gpkg_out, filetype="GPKG", overwrite=TRUE)
# 
# 
# goode_crs <- "+proj=igh +datum=WGS84 +no_defs"
# r_ll <- rast(xmin = -180, xmax = 180, ymin = -90, ymax = 90,
#              crs = "EPSG:4326", resolution = 1)
# values(r_ll) <- 0
# tmpl_goode <- project(
#   r_ll, goode_crs, method = "near",
#   res = 10000,
#   filename  = "/global/scratch/users/minoli/maps/data/goode_template.tif"
# )

countries <- ne_countries(scale = "medium", returnclass = "sf")$name_long
vague_terms <- c("Sea", "Ocean", "Bay", "Gulf", "Channel", "Coast", "Territory", "Province", "State")

df <- df %>%
  mutate(
    vague_location_flag = case_when(
      # coordinates exist (not vague)
      !is.na(latlon) ~ FALSE,
      
      # no coords + geo_location is kinda vague
      (is.na(latlon)) & (geo_location %in% countries |
                           str_detect(geo_location, str_c(vague_terms, collapse = "|"))) ~ TRUE,
      
      # geolocation is specific names but no coords
      TRUE ~ FALSE
    ),
    missing_sampling_location_flag = ifelse(
      is.na(latlon) & (is.na(geo_location) | geo_location == ""),
      TRUE, FALSE
    )
  )

write.csv(df, "data/vgp_iucn_loc_metadata.csv")

