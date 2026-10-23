#!/usr/bin/env Rscript
library(optparse)
library(terra)
library(sf)
library(data.table)

option_list <- list(
  make_option(dest="points", type="character", help="Path to points CSV"),
  make_option(dest="raster", type="character", help="Path to raster TIF"),
  make_option(dest="buffer", type="numeric", default=0.0, help="Buffer distance in meters"),
  make_option(dest="fun", type="character", default="mean", help="Summary function"),
  make_option(dest="out", type="character", help="Output CSV path")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Coordinate system config
ALBERS <- "ESRI:102003"

# 1. Load data
df <- fread(opt$points)
r <- rast(opt$raster)

# 2. Create spatial object
# Assuming input is WGS84 (EPSG:4326)
pts <- st_as_sf(df, coords = c("longitude", "latitude"), remove = FALSE, crs = 4326)
pts_proj <- st_transform(pts, crs = ALBERS)

# 3. Handle Buffering
if (opt$buffer > 0) {
  spatial_obj <- st_buffer(pts_proj, dist = opt$buffer)
} else {
  spatial_obj <- pts_proj
}

# 4. Projection alignment
if (!terra::same.crs(spatial_obj, r)) {
  spatial_obj <- st_transform(spatial_obj, crs = terra::crs(r))
}

# 5. Extract
# For buffer extraction, use the summary function. For point extraction, use simple.
exact <- opt$fun %in% c("mean", "min", "max")
res <- terra::extract(r, spatial_obj, fun = opt$fun, na.rm = TRUE, exact = exact)

# 6. Save (PatientID and the extracted value)
output <- data.frame(patientID = df$patientID, value = res[, 2])
fwrite(output, opt$out)
