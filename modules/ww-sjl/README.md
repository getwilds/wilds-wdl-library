# ww-sjl
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for Solar Jetlag (SJL) tile processing.

## Overview

This module provides a reusable WDL task for calculating sunrise/sunset times and sun time differences for geographic tiles, as part of the SJL model pipeline. It uses NOAA solar calculator variables to compute solar geometry over a full year, then matches tile points to timezone border points to produce sun time difference values.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `sjl_tiles`
- **Test workflow**: `testrun.wdl` (demonstration workflow that executes all tasks)
- **Script**: `sjl_tiles.R` (pulled from GitHub at runtime)
- **Container**: `rocker/tidyverse:4`

## Tasks

### `sjl_tiles`
Calculates sunrise/sunset times and sun time differences for a single geographic tile.

**Inputs:**
- `tile_path` (File): Input tile RDS file
- `border_points_path` (File): Border points CSV file containing timezone boundary data
- `tile_num` (String): Tile identifier (e.g. `0042`)
- `year` (Int): Year for solar calculations (e.g. `2022`)
- `cpu_cores` (Int, default: 1): Number of CPU cores to use
- `memory_gb` (Int, default: 8): Memory allocation in GB

**Outputs:**
- `matched_points` (File): RDS file containing points with sunrise/sunset difference values
- `missing_points` (File): RDS file containing points that could not be matched to border points (may be empty)

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sjl/ww-sjl.wdl" as ww_sjl

workflow my_sjl_pipeline {
  input {
    File tile_path
    File border_points_path
    String tile_num
    Int year
  }

  call ww_sjl.sjl_tiles {
    input:
      tile_path = tile_path,
      border_points_path = border_points_path,
      tile_num = tile_num,
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

This module fetches the `sjl_tiles.R` script from GitHub at runtime rather than bundling it inside the Docker container. This design simplifies development and customization but introduces a reproducibility consideration: the script and container are tracked separately.

**For production workflows**, we recommend one of the following approaches:

1. **Pin to a specific commit**: Replace the branch reference in the wget URL with a commit hash:
   ```
   # Instead of refs/heads/main, use a specific commit:
   https://raw.githubusercontent.com/getwilds/wilds-wdl-library/<commit-sha>/modules/ww-sjl/sjl_tiles.R
   ```

2. **Fork the repository**: Create your own fork and reference your stable branch or tagged release.

3. **Bundle the script**: For maximum reproducibility, build a custom Docker container that includes the script.

## Support and Documentation

For questions about this module:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines
