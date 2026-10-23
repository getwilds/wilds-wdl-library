# ww-geolink Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for geospatial linkage and exposure measure extraction. It provides atomic tasks to extract values from environmental rasters for specific points or buffers, and computes aggregated annual exposure metrics.

## Overview

This module solves the "bottleneck" of sequential raster extraction by providing granular tasks that can be scattered across a compute cluster. It is designed to work with the `ww-exposure` pipeline to generate high-resolution exposure datasets for large patient cohorts.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library):

- **Main WDL file**: `ww-geolink.wdl` - Contains task definitions
- **Worker Scripts**: 
    - `extract_val.R`: Point/buffer extraction logic
    - `compute_wf_metrics.R`: Wildfire annual metrics and lags
    - `join_all.R`: Final dataset aggregation
- **Docker Config**: `Dockerfile_latest` - Specification for the `r-spatial` image

## Available Tasks

### `extract_raster`

Extracts a value or statistic from a single TIF raster for a set of points or buffers.

**Inputs:**
- `points_csv` (File): CSV containing patient coordinates (latitude/longitude)
- `raster_tif` (File): TIF raster file
- `buffer_dist` (Float, default=0.0): Buffer radius in meters; 0.0 performs a simple point extraction
- `stat_fun` (String, default="mean"): Summary function for buffered areas (e.g., "mean", "max", "min")
- `output_prefix` (String, default="extract_"): Prefix for the resulting CSV

**Outputs:**
- `result_csv` (File): CSV containing `patientID` and the extracted `value`

### `compute_wildfire_metrics`

Computes specialized annual wildfire exposure metrics and 3-year lags from daily extracts.

**Inputs:**
- `daily_extracts` (Array[File]): Array of daily exposure CSVs
- `p90_pctile` (Float, default=0.90): Percentile threshold for smoke-wave events
- `min_days` (Int, default=2): Minimum consecutive days to qualify as a smoke-wave

**Outputs:**
- `annual_metrics` (File): CSV containing aggregated annual metrics per patient

### `join_exposures`

Aggregates fragmented extraction CSVs and annual metrics into a final wide-format dataset.

**Inputs:**
- `patient_base` (File): Base patient demographics/coordinates CSV
- `all_extracts` (Array[File]): All extracted daily exposure CSVs
- `annual_wf_metrics` (File): The output from `compute_wildfire_metrics`

**Outputs:**
- `final_dataset` (File): The final joined exposure dataset

## Docker Container

This module requires the `getwilds/r-spatial` image (based on `rocker/geospatial`), which includes:
- **R Spatial Stack**: `sf`, `terra`, `raster`
- **Data Utilities**: `data.table`, `dplyr`, `lubridate`
- **System Dependencies**: GDAL, GEOS, PROJ
