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

few_summaries <- all_summaries %>%  filter(species %in% c("Rhinolophus luctus", 
                                                          "Acridotheres tristis", 
                                                          "Pseudophryne corroboree", 
                                                          "Poromitra crassiceps",
                                                          "Alca torda"))

# check_outside_flag <- function(lon, lat, rpath, buffer_km = 50, debug = FALSE) {
#   flag <- NA
#   
#   if (!is.na(lon) && !is.na(lat) && !is.na(rpath) && file.exists(rpath)) {
#     try({
#       r <- terra::rast(rpath)
#       r_mask <- r > 0
#       
#       res_km  <- mean(res(r)) / 1000
#       n_cells <- max(1, ceiling(buffer_km / res_km))
#       n_cells <- min(n_cells, floor(min(nrow(r), ncol(r)) / 2))  
#       
#       w <- matrix(1, nrow = 2 * n_cells + 1, ncol = 2 * n_cells + 1)
#       r_dilated <- terra::focal(r_mask, w = w, fun = max)
#       
#       pt <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326)
#       pt_r <- sf::st_transform(pt, crs(r))
#       coords <- sf::st_coordinates(pt_r)
#       
#       ext <- ext(r_dilated)
#       inside <- coords[1] >= ext[1] && coords[1] <= ext[2] &&
#         coords[2] >= ext[3] && coords[2] <= ext[4]
#       
#       if (!inside) {
#         flag <- TRUE  
#       } else {
#         val_dil <- terra::extract(r_dilated, coords)[,1]
#         if (is.na(val_dil)) {
#           flag <- TRUE
#         } else {
#           flag <- val_dil == 0
#         }
#       }
#       
#       if (debug) {
#         message(sprintf(
#           "lon=%.4f  lat=%.4f  projX=%.1f  projY=%.1f  inside=%s  val=%s -> flag=%s",
#           lon, lat, coords[1], coords[2], inside,
#           ifelse(exists("val_dil"), val_dil, "NA"), flag
#         ))
#       }
#     }, silent = TRUE)
#   }
#   
#   return(flag)
# }
# 
# .outside_check_impl <- function(lon, lat, rpath, buffer_km = 50) {
#   r <- terra::rast(rpath)
#   
#   mask <- r > 0
#   mask <- terra::ifel(is.na(mask), 0, mask)
#   
#   poly <- try(terra::as.polygons(mask, dissolve = TRUE), silent = TRUE)
#   if (inherits(poly, "try-error"))             return(list(flag = TRUE, reason = "mask_to_polygon_failed"))
#   poly <- poly[poly[[1]] == 1, , drop = FALSE]
#   if (nrow(poly) == 0)                         return(list(flag = TRUE, reason = "empty_mask_polygon"))
#   
#   ea_crs <- "EPSG:6933"
#   poly_ea  <- try(terra::project(poly, ea_crs), silent = TRUE)
#   if (inherits(poly_ea, "try-error"))          return(list(flag = TRUE, reason = "poly_project_failed"))
#   buf_m    <- buffer_km * 1000
#   poly_buf <- try(terra::buffer(poly_ea, width = buf_m), silent = TRUE)
#   if (inherits(poly_buf, "try-error"))         return(list(flag = TRUE, reason = "buffer_failed"))
#   
#   # 4) point-in-buffer test
#   pt_wgs <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326)
#   pt_ea  <- try(sf::st_transform(pt_wgs, ea_crs), silent = TRUE)
#   if (inherits(pt_ea, "try-error"))            return(list(flag = TRUE, reason = "point_transform_failed"))
#   
#   inside <- as.logical(sf::st_intersects(pt_ea, sf::st_as_sf(poly_buf), sparse = FALSE))
#   if (!inside) return(list(flag = TRUE,  reason = "outside_after_buffer"))
#   list(flag = FALSE, reason = "inside_after_buffer")
# }

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

few_summaries <- few_summaries %>%
  rowwise() %>%
  mutate(
    outside_range_flag = check_outside_flag(sampling_lon, sampling_lat, raster_path_ea),
    outside_reason     = check_outside_reason(sampling_lon, sampling_lat, raster_path_ea)
  ) %>%
  ungroup()

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

all_species <- (df %>% distinct(iucn_name, .keep_all = TRUE)) %>%
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
  )


few_all_species <- df %>%
  left_join(
    few_summaries %>%
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
  )


# write_csv(all_species, "results/summary/all_species_summary.csv")
