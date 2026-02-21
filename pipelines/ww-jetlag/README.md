# ww-jetlag Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline for running Solar Jetlag (SJL) tile processing across an array of geographic tiles.

## Overview

This pipeline runs the `sjl_tiles` task from the `ww-sjl` module in parallel across an array of input tiles, using a shared set of timezone border points. It is designed as a scalable wrapper around the `ww-sjl` module for processing multiple tiles as part of the SJL model pipeline.

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Main WDL file**: `ww-jetlag.wdl` - Workflow definition
- **Test workflow**: `testrun.wdl` - Demonstration workflow using test data
- **Example inputs**: `inputs.json` - Starting point for input configuration

## Module Dependencies

This pipeline imports and uses:
- **ww-sjl module**: For SJL tile processing (`sjl_tiles` task)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access (to pull the `sjl_tiles.R` script from GitHub at runtime)

### Input Configuration

Edit `inputs.json` with your tile files and settings:

```json
{
  "jetlag.tile_paths": ["/path/to/tile_0001.rds", "/path/to/tile_0002.rds"],
  "jetlag.tile_nums": ["0001", "0002"],
  "jetlag.border_points_path": "/path/to/border_points.csv",
  "jetlag.year": 2022,
  "jetlag.cpu_cores": 1,
  "jetlag.memory_gb": 8
}
```

**Note**: `tile_paths` and `tile_nums` must be the same length and in the same order — each `tile_nums[i]` must correspond to `tile_paths[i]`.

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-jetlag.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-jetlag.wdl -i inputs.json

# Using Sprocket
sprocket run ww-jetlag.wdl inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this pipeline directly to the on-premise HPC cluster.

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `tile_paths` | Array of input tile RDS files | Array[File] | Yes | - |
| `tile_nums` | Array of tile identifiers corresponding to each tile_path | Array[String] | Yes | - |
| `border_points_path` | Border points CSV file shared across all tiles | File | Yes | - |
| `year` | Year for solar calculations | Int | Yes | - |
| `cpu_cores` | Number of CPU cores per tile task | Int | No | 1 |
| `memory_gb` | Memory allocation in GB per tile task | Int | No | 8 |

## Output Files

| Output | Description |
|--------|-------------|
| `matched_points` | Array of RDS files containing points with sunrise/sunset difference values, one per input tile |
| `missing_points` | Array of RDS files containing points that could not be matched to border points, one per input tile |

## Testing the Pipeline

The pipeline includes a test workflow that can be run independently:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl --entrypoint jetlag_example

# Using Sprocket
sprocket run testrun.wdl --entrypoint jetlag_example
```

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
