library(dplyr)
library(readr)
library(fs)
library(terra)
library(sf)
library(stringr)

df <- read.csv("data/vgp_iucn_loc_metadata.csv")
summary_files <- dir_ls("results/species", recurse = TRUE,
                        glob = "*/species_summary.csv")

all_summaries <- summary_files |>
  lapply(\(p) read_csv(p, show_col_types = FALSE)) |>   
  bind_rows() |>
  distinct()

# few_summaries <- all_summaries %>%  filter(species %in% c("Rhinolophus luctus", 
#                                                           "Acridotheres tristis", 
#                                                           "Pseudophryne corroboree", 
#                                                           "Poromitra crassiceps",
#                                                           "Alca torda"))

.outside_check_impl <- function(lon, lat, rpath, buffer_km = 50) {
  r <- terra::rast(rpath)
  
  mask <- r > 0
  mask <- terra::ifel(is.na(mask), 0, mask)
  
  poly <- try(terra::as.polygons(mask, dissolve = TRUE), silent = TRUE)
  if (inherits(poly, "try-error")) return(list(flag = TRUE,  reason = "mask_to_polygon_failed"))
  poly <- poly[poly[[1]] == 1, , drop = FALSE]
  if (nrow(poly) == 0)          return(list(flag = TRUE,  reason = "empty_mask_polygon"))
  
  ea_crs   <- "EPSG:6933"
  poly_ea  <- try(terra::project(poly, ea_crs), silent = TRUE)
  if (inherits(poly_ea, "try-error")) return(list(flag = TRUE,  reason = "poly_project_failed"))
  
  poly_buf <- try(terra::buffer(poly_ea, width = buffer_km * 1000), silent = TRUE)
  if (inherits(poly_buf, "try-error")) return(list(flag = TRUE,  reason = "buffer_failed"))
  
  pt_wgs <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326)
  pt_ea  <- try(sf::st_transform(pt_wgs, ea_crs), silent = TRUE)
  if (inherits(pt_ea, "try-error")) return(list(flag = TRUE,  reason = "point_transform_failed"))
  
  inside <- as.logical(sf::st_intersects(pt_ea, sf::st_as_sf(poly_buf), sparse = FALSE))
  
  if (inside) list(flag = FALSE, reason = "inside_after_buffer")
  else        list(flag = TRUE,  reason = "outside_after_buffer")
}

check_outside_flag <- function(lon, lat, rpath, buffer_km = 50) {
  if (!is.finite(lon) || !is.finite(lat))  return(NA)                # no coords
  if (is.na(rpath) || !file.exists(rpath)) return(NA)                # no range
  .outside_check_impl(lon, lat, rpath, buffer_km)$flag
}

check_outside_reason <- function(lon, lat, rpath, buffer_km = 50) {
  if (!is.finite(lon) || !is.finite(lat))  return("no_coords")
  if (is.na(rpath) || !file.exists(rpath)) return("no_range")
  .outside_check_impl(lon, lat, rpath, buffer_km)$reason
}

# check_outside_flag   <- function(lon, lat, rpath, buffer_km = 50) .outside_check_impl(lon, lat, rpath, buffer_km)$flag
# check_outside_reason <- function(lon, lat, rpath, buffer_km = 50) .outside_check_impl(lon, lat, rpath, buffer_km)$reason

all_summaries <- all_summaries %>%
  mutate(
    sampling_lon = suppressWarnings(as.numeric(sampling_lon)),
    sampling_lat = suppressWarnings(as.numeric(sampling_lat))
  ) %>%
  rowwise() %>%
  mutate(
    outside_range_flag   = check_outside_flag(sampling_lon, sampling_lat, raster_path_ea),
    outside_range_reason = check_outside_reason(sampling_lon, sampling_lat, raster_path_ea)
  ) %>%
  ungroup()

# few_summaries <- few_summaries %>%
#   rowwise() %>%
#   mutate(
#     outside_range_flag = check_outside_flag(sampling_lon, sampling_lat, raster_path_ea),
#     outside_reason     = check_outside_reason(sampling_lon, sampling_lat, raster_path_ea)
#   ) %>%
#   ungroup()

zoo_terms <- c("Zoo", "Aquarium", "Museum", "Collection")

domesticated_species <- c(
  "Bos taurus",            # cattle
  "Gallus gallus",         # chicken
  "Canis lupus",           # dog
  "Felis catus",           # cat
  "Equus caballus",        # horse
  "Sus scrofa",            # pig
  "Capra hircus",          # goat
  "Ovis aries"             # sheep
)

# all_species <- (df %>% distinct(iucn_name, .keep_all = TRUE)) %>%
#   left_join(
#     (all_summaries %>% distinct(species, .keep_all = TRUE)) %>%
#       dplyr::select(species, raster_path_ea, raster_path_goode,
#                     centroid_lon, centroid_lat, 
#                     outside_range_flag),
#     by = c("iucn_name" = "species")
#   ) %>%
#   mutate(
#     missing_range_flag = is.na(raster_path_ea),
#     zoo_flag = !is.na(geo_location) &
#       str_detect(geo_location, str_c(zoo_terms, collapse = "|")),
#     domesticated_flag = iucn_name %in% domesticated_species
#   )

all_species <- (df%>% distinct(scientific_name, .keep_all = TRUE)) %>%
  left_join(
    (all_summaries %>% distinct(species, .keep_all = TRUE)) %>%
      dplyr::select(species, raster_path_ea, raster_path_goode,
                    centroid_lon, centroid_lat, 
                    outside_range_flag),
    by = c("iucn_name" = "species")
  ) %>%
  mutate(
    missing_range_flag = is.na(raster_path_ea),
    zoo_flag = !is.na(geo_location) &
      str_detect(geo_location, str_c(zoo_terms, collapse = "|")),
    domesticated_flag = iucn_name %in% domesticated_species
  ) %>%
  dplyr::select(-c(X, X.1))


write_csv(all_species, "results/summary/all_species_summary.csv")
