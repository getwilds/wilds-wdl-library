# ww-exposure Pipeline
[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline for massively parallel environmental exposure linkage, replacing sequential raster extraction with a scattered workflow.

## Overview

This pipeline links patient coordinates to a large library of environmental rasters (e.g. daily PM2.5 or wildfire smoke). A nested scatter extracts values for every combination of raster and buffer distance in parallel, then computes annual exposure metrics and joins everything into a single wide-format dataset.

**Complexity level:** Basic (2 modules)

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library).

## Pipeline Steps

1. **Optional S3 download** (using `ww-aws-sso` module): downloads the patient CSV and/or rasters from S3.
2. **Parallel extraction** (using `ww-geolink` module): one `extract_raster` task per raster per buffer distance.
3. **Metric calculation** (using `ww-geolink` module): `compute_wildfire_metrics` computes annual metrics per patient and buffer distance.
4. **Final join** (using `ww-geolink` module): `join_exposures` aggregates all extracts and metrics into one row per patient.

## Module Dependencies

This pipeline imports the following modules:
- **ww-geolink**: spatial extraction, metric calculation, and aggregation
- **ww-aws-sso**: optional S3 downloads

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker or Apptainer/Singularity support
- Internet access for fetching scripts from GitHub at runtime
- Patient CSV with `patientID`, `longitude`, `latitude` (WGS84)
- Daily TIF rasters whose file names contain the date (`YYYY-MM-DD` or `YYYYMMDD`)

### Input Configuration

Provide the patient file and the rasters either locally or via S3. If both a local value and an S3 value are given, the local value is used. See [inputs.json](inputs.json) for a template:

```json
{
  "exposure_linkage.patient_csv": "/path/to/patients.csv",
  "exposure_linkage.rasters_local": [
    "/path/to/pm25_2020-01-01.tif",
    "/path/to/pm25_2020-01-02.tif"
  ],
  "exposure_linkage.buffer_distances": [0.0, 1000.0, 5000.0]
}
```

For S3 inputs, replace `patient_csv` and `rasters_local` with `patient_csv_s3` and `rasters_manifest_s3`, and optionally provide `aws_config` and `aws_creds` for private buckets (public buckets are accessed without signing).

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-exposure.wdl -i inputs.json

# Using miniWDL
miniwdl run ww-exposure.wdl -i inputs.json

# Using Sprocket
sprocket run ww-exposure.wdl inputs.json
```

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `patient_csv` | Local patient CSV (`patientID`, `longitude`, `latitude`) | File | One of `patient_csv` or `patient_csv_s3` | - |
| `patient_csv_s3` | S3 URI of the patient CSV | String | One of `patient_csv` or `patient_csv_s3` | - |
| `rasters_local` | Local TIF rasters | Array[File] | One of `rasters_local` or `rasters_manifest_s3` | - |
| `rasters_manifest_s3` | S3 URIs of TIF rasters | Array[String] | One of `rasters_local` or `rasters_manifest_s3` | - |
| `buffer_distances` | Radii in meters (0.0 = point extraction) | Array[Float] | No | `[0.0, 1000.0, 5000.0, 10000.0]` |
| `aws_config` | AWS config file for S3 access | File | No | - |
| `aws_creds` | AWS credentials file for S3 access | File | No | - |

## Output Files

| Output | Description |
|--------|-------------|
| `final_dataset` | Wide-format CSV with the base patient data joined to one column per raster and buffer distance, plus annual exposure metrics per year and buffer distance |

## Resource Considerations

This pipeline creates a large number of small tasks:
- **Task count**: the number of extraction tasks is the number of rasters multiplied by the number of buffer distances. For 100 rasters and 4 buffers, 400 tasks are launched.
- **Disk I/O**: ensure the underlying filesystem can handle many concurrent reads from S3 or shared storage.
- **Final join**: `join_exposures` produces one column per raster and buffer, so very large raster libraries yield very wide tables.

## Pipeline Testing

The pipeline includes a test workflow that runs with zero configuration, using synthetic data from `ww-testdata`:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
