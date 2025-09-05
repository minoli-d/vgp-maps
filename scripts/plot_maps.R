#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(terra)
  library(dplyr)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(viridis)
  library(raster)
  library(sp)
})

summary_csv <- "results/summary/all_species_summary.csv"
summary_plots_dir <- "results/summary/plots"
dir.create(summary_plots_dir, recursive=TRUE, showWarnings=FALSE)

df <- read.csv(summary_csv, stringsAsFactors=FALSE)
df <- df %>% filter(!is.na(centroid_lon), !is.na(sampling_lon))

# goode homolosine setup
lats <- c(
  90:-90,             # right edge down
  -90:0,  0:-90,      # third cut
  -90:0,  0:-90,      # second cut
  -90:0,  0:-90,      # first cut
  -90:90,             # left edge up
  90:0,   0:90,       # top cut
  90                  # close
)
longs <- c(
  rep(180, 181),                    # right edge
  rep(c(80.01, 79.99), each = 91),  # third cut
  rep(c(-19.99, -20.01), each = 91),# second cut
  rep(c(-99.99, -100.01), each = 91),# first cut
  rep(-180, 181),                   # left edge
  rep(c(-40.01, -39.99), each = 91),# top cut
  180                               # close
)
goode_outline_ll <- st_polygon(list(cbind(longs, lats))) %>%
  st_sfc(crs = "+proj=longlat +datum=WGS84")
crs_goode <- "+proj=igh +datum=WGS84 +no_defs"
goode_outline   <- st_transform(goode_outline_ll, crs_goode)
bb   <- st_bbox(goode_outline)
xlim <- bb[c("xmin","xmax")] * c(1.1,1.1)
ylim <- bb[c("ymin","ymax")] * c(1.1,1.1)
goode_encl_rect <- st_polygon(list(rbind(
  c(xlim[1], ylim[1]),
  c(xlim[2], ylim[1]),
  c(xlim[2], ylim[2]),
  c(xlim[1], ylim[2]),
  c(xlim[1], ylim[1])
))) %>% st_sfc(crs = crs_goode)
goode_without <- st_difference(goode_encl_rect, goode_outline)
world_ll   <- ne_countries(scale="medium", returnclass="sf")
world_goode<- st_transform(world_ll, crs_goode)
ocean_goode<- st_difference(goode_encl_rect, st_union(world_goode))

tmpl_goode <- terra::rast(
  extent     = terra::ext(xlim[1], xlim[2], ylim[1], ylim[2]),
  resolution = 10000,                     # ~20 km in Goode meters; bump if needed
  crs        = crs_goode
)


df <- read.csv("results/summary/all_species_summary.csv", stringsAsFactors=FALSE) %>%
  filter(!is.na(sampling_lon), !is.na(centroid_lon))

crs_goode <- "+proj=igh +datum=WGS84 +no_defs"
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs_goode)

# sampling locations
sampling_pts <- st_as_sf(df, coords = c("sampling_lon", "sampling_lat"), crs = 4326) %>%
  st_transform(crs_goode)

p1 <- ggplot() +
  geom_sf(data = ocean_goode, fill="white", color=NA) +
  geom_sf(data = world_goode, fill="grey85", color="white", size=0.2) +
  geom_sf(data = goode_without, fill="white", color=NA) +
  geom_sf(data = sampling_pts, color="red3", size=1, shape=21, fill=NA) +
  geom_sf(data = goode_outline, fill=NA, color="gray70", size=0.5) +
  theme_minimal() +
  labs(title="Sampling Locations")

ggsave(file.path(summary_plots_dir, "141_sampling_locations_map.png"), p1, width = 10, height = 6)


# range centroids 
centroids <- st_as_sf(df, coords = c("centroid_lon", "centroid_lat"), crs = 4326) %>%
  st_transform(crs_goode)

p2 <- ggplot() +
  geom_sf(data = ocean_goode, fill="white", color=NA) +
  geom_sf(data = world_goode, fill="grey85", color="white", size=0.2) +
  geom_sf(data = goode_without, fill="white", color=NA) +
  geom_sf(data = centroids, color="royalblue3", size=1, shape=21, fill=NA) +
  geom_sf(data = goode_outline, fill=NA, color="gray70", size=0.5) +
  theme_minimal() +
  labs(title="Range Centroids")

ggsave(file.path(summary_plots_dir, "141_range_centroids_map.png"), p2, width = 10, height = 6)


# heatmap
aligned_files <- character(length(df$raster_path))
for (i in seq_along(df$raster_path)) {
  r  <- rast(df$raster_path[i])
  pr <- project(
    r, tmpl_goode, method = "near",
    filename  = file.path(tempdir(), sprintf("al_%03d.tif", i)),
    overwrite = TRUE,
    wopt      = list(datatype = "INT2U", gdal = c("COMPRESS=DEFLATE"))
  )
  # treat NA as 0 without pulling whole raster into RAM
  pr0 <- ifel(is.na(pr), 0, pr)
  aligned_files[i] <- writeRaster(
    pr0, filename = file.path(tempdir(), sprintf("al0_%03d.tif", i)),
    overwrite = TRUE,
    wopt      = list(datatype = "INT2U", gdal = c("COMPRESS=DEFLATE"))
  )
}

# 2) Stack aligned layers and sum across them (streamed to disk)
x <- rast(aligned_files)  # geometries now match
sum_raster_goode <- app(
  x, sum, na.rm = TRUE,
  filename  = file.path(tempdir(), "sum_goode.tif"),
  overwrite = TRUE,
  wopt      = list(datatype = "INT4U", gdal = c("COMPRESS=DEFLATE"))
)

# 3) To a dataframe (if needed)
r_df <- as.data.frame(sum_raster_goode, xy = TRUE, na.rm = FALSE)
names(r_df)[3] <- "count"

p3 <- ggplot() +
  geom_sf(data = ocean_goode,       fill = "white",  color = NA) +
  geom_sf(data = world_goode,       fill = "grey85", color = "white", size = 0.2) +
  geom_raster(data = r_df, aes(x = x, y = y, fill = count), interpolate = TRUE) +
  scale_fill_viridis_c(name = "Count") +
  geom_sf(data = goode_without,     fill = "white",  color = NA) +
  geom_sf(data = goode_outline,     fill = NA,       color = "gray70", size = 0.5) +
  coord_sf(crs = crs_goode, expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(title = "Species Range Density", x=NULL, y=NULL)

ggsave(file.path(summary_plots_dir, "141_range_heatmap.png"), p3, width = 10, height = 6)


# individual species maps

for (i in seq_len(nrow(df))) {
  sp <- df$species[i]
  sp_dir <- file.path("results/species", gsub(" ", "_", sp))
  dir.create(sp_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    r0 <- terra::rast(df$raster_path[i])
    r  <- terra::project(r0, tmpl_goode, method = "near",
                         filename = "", overwrite = TRUE) 
    
    r_df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    names(r_df)[3] <- "presence"
    r_df <- r_df[r_df$presence > 0, , drop = FALSE]
    
    sp_pts <- data.frame(
      name = c("Sampling", "Centroid"),
      lon  = c(df$sampling_lon[i], df$centroid_lon[i]),
      lat  = c(df$sampling_lat[i], df$centroid_lat[i])
    )
    sp_sf <- sf::st_as_sf(sp_pts, coords = c("lon","lat"), crs = 4326) |>
      sf::st_transform(crs_goode)
    
    p <- ggplot() +
      geom_sf(data = ocean_goode,   fill = "white",  color = NA) +
      geom_sf(data = world_goode,   fill = "grey85", color = "white", size = 0.2) +
      geom_sf(data = goode_without, fill = "white",  color = NA) +
      geom_sf(data = goode_outline, fill = NA,       color = "gray70", size = 0.5) +
      geom_tile(data = r_df, aes(x = x, y = y), fill = "chartreuse3", alpha = 0.5) +
      geom_sf(data = sp_sf[sp_sf$name == "Sampling", ], shape = 21, fill = "red3",      size = 2) +
      geom_sf(data = sp_sf[sp_sf$name == "Centroid", ], shape = 21, fill = "royalblue3", size = 2) +
      coord_sf(crs = crs_goode, expand = FALSE) +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), axis.title = element_blank()) +
      labs(title = paste("Range and Sample for", sp))
    
    ggsave(file.path(sp_dir, paste0(gsub(" ", "_", sp), "_map.png")),
           p, width = 8, height = 5, dpi = 200)
  }, error = function(e) {
    message("Plot failed for ", sp, ": ", conditionMessage(e))
  })
}