# ww-sjl
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for calculating sunrise and sunset time differences across an array of geographic tiles as part of the Solar Jetlag (SJL) project.

## Overview

Solar jetlag is a high-resolution geospatial model that measures circadian misalignment from geographic variation in light exposure timing across time zones. This module provides a reusable WDL task for calculating sunrise/sunset times and solar time differences across geographic tiles as part of the SJL model pipeline. Sunrise and sunset times are calculated using NOAA solar calculator variables averaged over a full year; each tile point is then matched to its time zone's easternmost boundary at the corresponding latitude to derive solar time difference values (measured in seconds).

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `sjl_tiles`
- **Test workflow**: `testrun.wdl` (demonstration workflow that executes all tasks)
- **Script**: `sjl_tiles.R` (pulled from GitHub at runtime)
- **Container**: `rocker/tidyverse:4`

## Tasks

### `sjl_tiles`
Calculates sunrise/sunset times and sun time differences for each point within a single geographic tile.

**Inputs:**
- `tile_path` (File): Input tile RDS file (any filename is accepted). Input tiles contain points spaced 30 m apart with point location and timezone. 
- `border_points_path` (File): Border points CSV file containing point location, sunrise/ sunset time, and time zone data.
- `year` (Int): Year for solar calculations (e.g. `2022`) point 
- `matched_prefix` (String, default: `"matched_"`): Filename prefix for matched results output. Set to `""` to keep the original input filename.
- `missing_prefix` (String, default: `"missing_"`): Filename prefix for missing results output. Set to `""` to keep the original input filename.
- `cpu_cores` (Int, default: 1): Number of CPU cores to use
- `memory_gb` (Int, default: 8): Memory allocation in GB

**Outputs:**
- `matched_points` (File): RDS file containing points with sunrise/sunset difference values (in seconds) (named `<matched_prefix><input_filename>.rds`)
- `missing_points` (File): RDS file containing points that could not be matched to border points, may be empty (named `<missing_prefix><input_filename>.rds`)

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sjl/ww-sjl.wdl" as ww_sjl

workflow my_sjl_pipeline {
  input {
    File tile_path
    File border_points_path
    Int year
  }

  call ww_sjl.sjl_tiles {
    input:
      tile_path = tile_path,
      border_points_path = border_points_path,
      year = year
  }

  output {
    File matched_points = sjl_tiles.matched_points
    File missing_points = sjl_tiles.missing_points
  }
}
```

## Testing the Module

The module includes a test workflow that demonstrates the `sjl_tiles` task:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl --entrypoint sjl_example

# Using Sprocket
sprocket run testrun.wdl --entrypoint sjl_example
```

## Reproducibility

This module uses the latest version of `sjl_tiles.R` at runtime. This script is not bundled inside the Docker container. For maximum reproducibility, we recommend one of the following approaches:

1. Use a commit hash instead of a branch reference in the `wget` command of your own copy of `ww-sjl.wdl`. You can make your own copy by forking the repo or downloading the script. For example: `https://raw.githubusercontent.com/getwilds/wilds-wdl-library/<COMMIT HASH HERE>/modules/ww-sjl/sjl_tiles.R`
2. Build a custom Docker container that includes the script inside.

## Acknowledgments

This module was developed in collaboration with [@cnondin](https://github.com/cnondin) and the [VoPham Lab](https://www.geoexlab.com/) at Fred Hutch through the [WILDS WDL Development Program](https://sciwiki.fredhutch.org/datascience/wilds_workflow_dev/). For questions about the underlying science, please reach out to the VoPham Lab through their website. Thank you to the VoPham Lab for their contributions and domain expertise!

## Support and Documentation

For questions about this module:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for detailed guidelines
- For questions about the SJL model, please contact Caroline Nondin, MS at [cnondin@fredhutch.org](mailto:cnondin@fredhutch.org), or Trang VoPham, PhD, MS at [trang@fredhutch.org](mailto:trang@fredhutch.org)
