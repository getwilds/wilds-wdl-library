# ww-geolink
[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for geospatial linkage: extracting environmental raster values for patient locations and computing aggregated annual exposure metrics.

## Overview

This module provides granular tasks that extract a value (or summary statistic) from a single raster for a set of patient coordinates, compute annual wildfire smoke exposure metrics from daily extracts, and join everything into one wide-format dataset. Because each extraction is a separate task, extractions can be scattered across a compute cluster rather than run sequentially. The `ww-exposure` pipeline demonstrates this pattern.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `extract_raster`, `compute_wildfire_metrics`, `join_exposures`
- **Scripts**: `extract_val.R`, `compute_wf_metrics.R`, `join_all.R` (fetched via wget at runtime)
- **Test workflow**: `testrun.wdl` (demonstration workflow using synthetic patients and rasters)
- **Containers**: `getwilds/r-spatial:0.1.0`
- **Dependencies**: Integrates with the `ww-testdata` module (`generate_test_geodata`) for test data generation

## Scripts

All scripts are fetched from this repository at runtime via `wget`. This keeps the Docker image generic and makes it easy to iterate on analysis logic without rebuilding containers.

| Script | Used by | Language | Description |
|--------|---------|----------|-------------|
| [`extract_val.R`](extract_val.R) | `extract_raster` | R | Projects patient points (optionally buffered in meters) to the raster CRS and extracts values with `terra` |
| [`compute_wf_metrics.R`](compute_wf_metrics.R) | `compute_wildfire_metrics` | R | Computes annual metrics per patient and buffer distance from daily extracts |
| [`join_all.R`](join_all.R) | `join_exposures` | R | Pivots extracts and annual metrics wide and joins them to the base patient table |

## Data Conventions

- **Patient CSV**: must contain `patientID`, `longitude`, and `latitude` (WGS84, EPSG:4326). Any additional columns are carried through to the final dataset.
- **Extract CSVs** (output of `extract_raster`): columns `patientID`, `raster` (raster file name without extension), `buffer_m`, and `value`. Files are named `<output_prefix>_<raster>_buf<N>m.csv`.
- **Dates**: `compute_wildfire_metrics` currently parses the date (`YYYY-MM-DD` or `YYYYMMDD`) from the raster file name, so daily rasters must be named accordingly (e.g. `pm25_2020-01-01.tif`). The task fails with an error if a date cannot be found.

## Tasks

### `extract_raster`

Extracts a value or statistic from a single raster for a set of points or buffers.

**Inputs:**
- `points_csv` (File): CSV with `patientID`, `longitude`, `latitude`
- `raster_tif` (File): TIF raster file
- `buffer_dist` (Float, default=0.0): Buffer radius in meters; 0.0 performs a simple point extraction
- `stat_fun` (String, default="mean"): Summary function for buffered areas (e.g. "mean", "max", "min")
- `output_prefix` (String, default="extract"): Prefix for the resulting CSV name
- `cpu_cores` (Int, default=1): Number of CPU cores
- `memory_gb` (Int, default=4): Memory in GB
- `docker_image` (String, default="getwilds/r-spatial:0.1.0"): Docker image

**Outputs:**
- `result_csv` (File): CSV with `patientID`, `raster`, `buffer_m`, and `value`

### `compute_wildfire_metrics`

Computes annual wildfire exposure metrics from daily extracts, separately for each buffer distance: `n_weeks_above_threshold` (weeks with mean above `weekly_threshold`), `n_days_gt0`, `peak_week_mean`, `n_smoke_waves`, `mean_daily` (annual mean), and `mean_lagged` (mean of the previous `lag_years` annual means).

**Inputs:**
- `daily_extracts` (Array[File]): Extract CSVs from `extract_raster`, one per raster and buffer distance
- `p90_pctile` (Float, default=0.90): Percentile threshold for smoke-wave events
- `min_days` (Int, default=2): Minimum consecutive days to qualify as a smoke wave
- `weekly_threshold` (Float, default=5.0): Weekly mean above which a week counts as high-exposure
- `lag_years` (Int, default=3): Number of prior years averaged for `mean_lagged`
- `cpu_cores` (Int, default=1): Number of CPU cores
- `memory_gb` (Int, default=8): Memory in GB
- `docker_image` (String, default="getwilds/r-spatial:0.1.0"): Docker image

**Outputs:**
- `annual_metrics` (File): CSV with one row per patient, year, and buffer distance

### `join_exposures`

Pivots all extracts to one column per raster and buffer distance, pivots the annual metrics to one column per metric, year, and buffer distance, and joins both to the base patient table (one row per patient).

**Inputs:**
- `patient_base` (File): Base patient CSV
- `all_extracts` (Array[File]): All extract CSVs from `extract_raster`
- `annual_wf_metrics` (File): Output of `compute_wildfire_metrics`
- `cpu_cores` (Int, default=1): Number of CPU cores
- `memory_gb` (Int, default=8): Memory in GB
- `docker_image` (String, default="getwilds/r-spatial:0.1.0"): Docker image

**Outputs:**
- `final_dataset` (File): Wide-format dataset with one row per patient

## Usage as a Module

### Importing into Your Workflow

```wdl
version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-geolink/ww-geolink.wdl" as ww_geolink

workflow my_exposure_analysis {
  input {
    File patients
    Array[File] rasters
  }

  scatter (raster in rasters) {
    call ww_geolink.extract_raster { input:
        points_csv = patients,
        raster_tif = raster,
        buffer_dist = 1000.0
    }
  }

  call ww_geolink.compute_wildfire_metrics { input:
      daily_extracts = extract_raster.result_csv
  }

  call ww_geolink.join_exposures { input:
      patient_base = patients,
      all_extracts = extract_raster.result_csv,
      annual_wf_metrics = compute_wildfire_metrics.annual_metrics
  }

  output {
    File final_dataset = join_exposures.final_dataset
  }
}
```

## Testing the Module

The module includes a test workflow that runs with zero configuration, using synthetic patients and daily rasters from `ww-testdata`:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

## Future Considerations

Design ideas that are intentionally not implemented yet, pending feedback on how collaborators use this module:

- **Explicit dates instead of filename parsing**: add an optional `raster_date` input to `extract_raster` that writes a `date` column to each extract CSV, so `compute_wildfire_metrics` no longer depends on dates in raster file names.
- **Chunked extraction**: for very large raster libraries, let each `extract_raster` task handle several rasters (or all buffer distances for one raster) to reduce per-task container startup and repeated reads of the patient CSV.

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker or Apptainer/Singularity support
- Internet access for fetching scripts from GitHub at runtime
- R spatial stack (`terra`, `sf`, `data.table`, `dplyr`, `lubridate`, `readr`, `tidyr`, `purrr`, `optparse`), provided by `getwilds/r-spatial:0.1.0`

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
