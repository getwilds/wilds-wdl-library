## WILDS WDL pipeline for Solar Jetlag (SJL) tile processing.
## Runs the sjl_tiles task across an array of geographic tiles
## using a shared set of timezone border points.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jetlag/modules/ww-sjl/ww-sjl.wdl" as ww_sjl

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
    tile_paths: "array of input tile RDS files to process"
    tile_nums: "array of tile identifiers corresponding to each tile_path (e.g. ['0001', '0002'])"
    border_points_path: "border points CSV file containing timezone boundary data, shared across all tiles"
    year: "year for solar calculations (e.g. 2022)"
    cpu_cores: "number of CPU cores to use per tile task"
    memory_gb: "memory allocation in GB per tile task"
  }

  input {
    Array[File] tile_paths
    Array[String] tile_nums
    File border_points_path
    Int year
    Int cpu_cores = 1
    Int memory_gb = 8
  }

  scatter (i in range(length(tile_paths))) {
    call ww_sjl.sjl_tiles { input:
      tile_path          = tile_paths[i],
      border_points_path = border_points_path,
      tile_num           = tile_nums[i],
      year               = year,
      cpu_cores          = cpu_cores,
      memory_gb          = memory_gb
    }
  }

  output {
    Array[File] matched_points = sjl_tiles.matched_points
    Array[File] missing_points = sjl_tiles.missing_points
  }
}
