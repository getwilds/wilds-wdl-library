#!/usr/bin/env Rscript
library(optparse)
library(terra)
library(sf)
library(data.table)

option_list <- list(
  make_option("--points", dest="points", type="character", help="Path to points CSV (patientID, longitude, latitude)"),
  make_option("--raster", dest="raster", type="character", help="Path to raster TIF"),
  make_option("--buffer", dest="buffer", type="numeric", default=0.0, help="Buffer distance in meters"),
  make_option("--fun", dest="fun", type="character", default="mean", help="Summary function"),
  make_option("--out", dest="out", type="character", help="Output CSV path")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Equal-area projection (meters) used for buffering
ALBERS <- "ESRI:102003"

# 1. Load data
df <- fread(opt$points)
r <- rast(opt$raster)

# 2. Create spatial object (input coordinates are assumed to be WGS84, EPSG:4326)
pts <- st_as_sf(df, coords = c("longitude", "latitude"), remove = FALSE, crs = 4326)
pts_proj <- st_transform(pts, crs = ALBERS)

# 3. Handle buffering
if (opt$buffer > 0) {
  spatial_obj <- st_buffer(pts_proj, dist = opt$buffer)
} else {
  spatial_obj <- pts_proj
}

# 4. Project into the raster's CRS (no-op if already matching)
spatial_obj <- st_transform(spatial_obj, crs = terra::crs(r))
spatial_vect <- terra::vect(spatial_obj)

# 5. Extract (summary function only applies to buffered polygons)
if (opt$buffer > 0) {
  res <- terra::extract(r, spatial_vect, fun = match.fun(opt$fun), na.rm = TRUE)
} else {
  res <- terra::extract(r, spatial_vect)
}

# 6. Save, recording which raster and buffer each value came from
output <- data.frame(
  patientID = df$patientID,
  raster = tools::file_path_sans_ext(basename(opt$raster)),
  buffer_m = opt$buffer,
  value = res[, 2]
)
fwrite(output, opt$out)
