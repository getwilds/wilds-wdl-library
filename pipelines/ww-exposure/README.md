# ww-exposure Pipeline

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline for massively parallel environmental exposure linkage. It transforms a sequential R-based extraction process into a distributed workflow, reducing runtime from days to hours.

## Overview

This pipeline automates the process of linking patient coordinates to a large library of environmental raster data (e.g., PM2.5, temperature, ozone). By using a nested scatter strategy, it extracts values for every combination of raster and buffer distance in parallel across a compute cluster.

## Pipeline Structure

This pipeline integrates the `ww-geolink` module:

1. **Parallel Extraction**: Scatters over all provided rasters and buffer distances.
2. **Metric Calculation**: Computes complex annual wildfire metrics from the extracted daily data.
3. **Final Join**: Aggregates all extracted features into a final wide-format CSV for downstream analysis.

## Module Dependencies

This pipeline depends on:
- **ww-geolink module**: For atomic spatial extractions and aggregation.

## Usage

### Input Configuration

Create an inputs JSON file specifying the patient list and the rasters to be processed:

```json
{
  "exposure_linkage.patient_csv": "/path/to/patients.csv",
  "exposure_linkage.rasters_manifest": [
    "/path/to/pm25_2020_01.tif",
    "/path/to/pm25_2020_02.tif",
    "/path/to/temp_2020.tif"
  ],
  "exposure_linkage.buffer_distances": [0.0, 1000.0, 5000.0]
}
```

### Running the Pipeline

```bash
# Using miniWDL
miniwdl run ww-exposure.wdl -i inputs.json

# Using Sprocket
sprocket run ww-exposure.wdl @inputs.json
```

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `patient_csv` | CSV containing patient IDs and WGS84 coordinates | File | Yes | - |
| `rasters_manifest` | Array of TIF files to process | Array[File] | Yes | - |
| `buffer_distances` | Array of radii in meters (0.0 = point extraction) | Array[Float] | No | `[0.0, 1000.0, 5000.0, 10000.0]` |

## Output Files

| Output | Description |
|--------|-------------|
| `final_dataset` | A wide-format CSV containing the base patient data joined with all extracted exposure values and annual metrics. |

## Resource Considerations

Because this pipeline creates a very large number of small tasks, users should be aware of:
- **Task Overhead**: The number of tasks is $|rasters| \times |buffers|$. For 100 rasters and 4 buffers, 400 tasks are launched.
- **Disk I/O**: Ensure the underlying filesystem can handle high concurrent read volume from S3 or shared storage.
