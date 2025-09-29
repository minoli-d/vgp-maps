suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

# message("args are", paste(args, collapse = " | "))

species         <- args[1]
shapefile       <- args[2]
metadata_csv    <- args[3]
output_dir      <- args[4]
tmpl_goode_path <- args[5]

df  <- read.csv(metadata_csv, stringsAsFactors = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("---", species, "---")

# sampling location (from metadata)
df_sp <- df %>%
  filter((scientific_name == species | iucn_name == species) & 
           !is.na(lon) & !is.na(lat))

if (nrow(df_sp) == 0) {
  sampling_lon <- NA
  sampling_lat <- NA
} else {
  sampling_lon <- df_sp %>% pull(lon)
  sampling_lat <- df_sp %>% pull(lat)
}
message(". Computed sampling coords")

ea_crs <- "EPSG:6933"
wgs_crs <- "EPSG:4326"
goode_crs <- "+proj=igh +datum=WGS84 +no_defs"
raster_resolution <- 10000

layer_name <- sf::st_layers(shapefile)$name[1]
esc <- function(x) gsub("'", "''", x)
qry <- sprintf("SELECT * FROM \"%s\" WHERE sci_name LIKE '%s'", layer_name, esc(species))
poly_sf <- sf::st_read(shapefile, query = qry, quiet = TRUE)

if (nrow(poly_sf) == 0) {
  message(species, " has no polygons")
  quit(save = "no", status = 1)
}

poly <- terra::vect(poly_sf)
rm(poly_sf); gc()
terraOptions(todisk = TRUE, memfrac = 0.6)

message(". Loaded shapefiles")

poly_diss <- aggregate(poly, by = "sci_name")
poly_simp <- simplifyGeom(poly_diss, tolerance = 0.05)

poly_proj_ea <- project(poly_simp, ea_crs)
range_area_km2 <- sum(expanse(poly_proj_ea, unit = "km"))

# weighted centroids for range
areas     <- expanse(poly_proj_ea, unit = "km")
centroids <- centroids(poly_proj_ea)
coords    <- crds(centroids)
weighted  <- colSums(coords * areas) / sum(areas)
weighted_sf <- st_as_sf(data.frame(x = weighted[1], y = weighted[2]), coords = c("x", "y"), crs = ea_crs)
weighted_wgs <- st_transform(weighted_sf, wgs_crs)
centroid_coords <- st_coordinates(weighted_wgs)
message(". Computed centroid")

# rasters for range
poly_proj_ea$pres <- 1
r_template_ea <- rast(poly_proj_ea, res = raster_resolution)
values(r_template_ea) <- 0
range_rast_ea <- rasterize(poly_proj_ea, r_template_ea,
                           field = "pres", background = 0)

raster_path_ea <- file.path(output_dir, "range_raster.tif")
writeRaster(range_rast_ea, raster_path_ea, overwrite = TRUE)
message(". Saved Equal Area raster")


tmpl_goode <- terra::rast(tmpl_goode_path)
poly_goode <- terra::project(poly_simp, terra::crs(tmpl_goode))
poly_goode$pres <- 1
raster_path_goode <- file.path(output_dir, "range_raster_goode.tif")
range_rast_goode <- terra::rasterize(
  poly_goode, tmpl_goode,
  field      = "pres", background = 0)
writeRaster(range_rast_goode, raster_path_goode, overwrite = TRUE)
message(". Saved Goode Homolosine raster")

# results as a df row
result <- data.frame(
  species            = species,
  centroid_lon       = centroid_coords[1],
  centroid_lat       = centroid_coords[2],
  n_range_polygons   = nrow(poly),
  range_area_km2     = range_area_km2,
  sampling_lon       = sampling_lon,
  sampling_lat       = sampling_lat,
  raster_path_ea     = raster_path_ea,
  raster_path_goode  = raster_path_goode,
  stringsAsFactors   = FALSE
)

row_path <- file.path(output_dir, "species_summary.csv")
write.csv(result, row_path, row.names = FALSE)
message(". Wrote result")
