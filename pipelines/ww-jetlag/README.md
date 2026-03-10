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

Tile files can be provided in two ways. Provide **one or the other**, not both.

**Option A: Manifest file** (recommended for large numbers of tiles)

Create a text file listing your tile RDS file paths, one per line:
```
/path/to/tile_0001.rds
/path/to/tile_0002.rds
/path/to/my_custom_tile.rds
```

Then reference it in `inputs.json`:
```json
{
  "jetlag.tile_manifest": "/path/to/tile_manifest.txt",
  "jetlag.border_points_path": "/path/to/border_points.csv",
  "jetlag.year": 2022,
  "jetlag.matched_prefix": "matched_",
  "jetlag.missing_prefix": "missing_",
  "jetlag.cpu_cores": 1,
  "jetlag.memory_gb": 8
}
```

**Tip**: Generate a manifest from a directory of tiles with `ls /path/to/tiles/*.rds > tile_manifest.txt`.

**Option B: Direct array** (convenient for smaller numbers of tiles)

```json
{
  "jetlag.tile_paths": ["/path/to/tile_0001.rds", "/path/to/tile_0002.rds"],
  "jetlag.border_points_path": "/path/to/border_points.csv",
  "jetlag.year": 2022,
  "jetlag.matched_prefix": "matched_",
  "jetlag.missing_prefix": "missing_",
  "jetlag.cpu_cores": 1,
  "jetlag.memory_gb": 8
}
```

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
| `tile_paths` | Array of input tile RDS files (provide this OR `tile_manifest`) | Array[File]? | No | - |
| `tile_manifest` | Text file listing tile paths, one per line (provide this OR `tile_paths`) | File? | No | - |
| `border_points_path` | Border points CSV file shared across all tiles | File | Yes | - |
| `year` | Year for solar calculations | Int | Yes | - |
| `matched_prefix` | Filename prefix for matched results output. Set to `""` to keep the original input filename. | String | No | `"matched_"` |
| `missing_prefix` | Filename prefix for missing results output. Set to `""` to keep the original input filename. | String | No | `"missing_"` |
| `cpu_cores` | Number of CPU cores per tile task | Int | No | 1 |
| `memory_gb` | Memory allocation in GB per tile task | Int | No | 8 |

## Output Files

| Output | Description |
|--------|-------------|
| `matched_points` | Array of RDS files containing points with sunrise/sunset difference values, one per input tile (named `<matched_prefix><input_filename>.rds`) |
| `missing_points` | Array of RDS files containing points that could not be matched to border points, one per input tile (named `<missing_prefix><input_filename>.rds`) |

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

## Acknowledgments

This pipeline was developed in collaboration with [@cnondin](https://github.com/cnondin) and the [VoPham Lab](https://www.geoexlab.com/) at Fred Hutch through the [WILDS WDL Development Program](https://sciwiki.fredhutch.org/datascience/wilds_workflow_dev/). For questions about the underlying science, please reach out to the VoPham Lab through their website. Thank you to the VoPham Lab for their contributions and domain expertise!

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
