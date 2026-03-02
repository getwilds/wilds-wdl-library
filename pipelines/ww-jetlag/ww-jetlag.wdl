## WILDS WDL pipeline for Solar Jetlag (SJL) tile processing.
## Runs the sjl_tiles task across an array of geographic tiles
## using a shared set of timezone border points.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sjl/ww-sjl.wdl" as ww_sjl

workflow jetlag {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL pipeline to calculate sunrise/sunset times and sun time differences across an array of geographic tiles as part of the SJL model pipeline"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-jetlag/ww-jetlag.wdl"
    outputs: {
      matched_points: "array of RDS files containing points with sunrise/sunset difference values, one per input tile",
      missing_points: "array of RDS files containing points that could not be matched to border points, one per input tile"
    }
  }

  parameter_meta {
    tile_paths: "array of input tile RDS files to process (provide this OR tile_manifest, not both)"
    tile_manifest: "text file listing input tile RDS file paths, one per line (provide this OR tile_paths, not both)"
    border_points_path: "border points CSV file containing timezone boundary data, shared across all tiles"
    year: "year for solar calculations (e.g. 2022)"
    matched_prefix: "filename prefix for matched results output (default: 'matched_')"
    missing_prefix: "filename prefix for missing results output (default: 'missing_')"
    cpu_cores: "number of CPU cores to use per tile task"
    memory_gb: "memory allocation in GB per tile task"
  }

  input {
    Array[File]? tile_paths
    File? tile_manifest
    File border_points_path
    Int year
    String matched_prefix = "matched_"
    String missing_prefix = "missing_"
    Int cpu_cores = 1
    Int memory_gb = 8
  }

  # Resolve tile files from whichever input was provided
  Array[File] resolved_tile_paths = if defined(tile_paths) then select_first([tile_paths]) else read_lines(select_first([tile_manifest]))

  scatter (tile_path in resolved_tile_paths) {
    call ww_sjl.sjl_tiles { input:
      tile_path          = tile_path,
      border_points_path = border_points_path,
      year               = year,
      matched_prefix     = matched_prefix,
      missing_prefix     = missing_prefix,
      cpu_cores          = cpu_cores,
      memory_gb          = memory_gb
    }
  }

  output {
    Array[File] matched_points = sjl_tiles.matched_points
    Array[File] missing_points = sjl_tiles.missing_points
  }
}
