version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/jetlag-manifest/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/jetlag-manifest/pipelines/ww-jetlag/ww-jetlag.wdl" as ww_jetlag

workflow jetlag_example {
  # Generate synthetic test tile and border points data
  call ww_testdata.generate_sjl_data { }

  # Create a manifest from the generated tile (the engine localizes the File
  # input so the path written into the manifest is valid for downstream tasks)
  call create_manifest { input:
    tile_file = generate_sjl_data.tile_rds
  }

  # Run the jetlag pipeline over a small array of test tiles
  call ww_jetlag.jetlag { input:
    tile_manifest      = create_manifest.manifest,
    border_points_path = generate_sjl_data.border_points_csv,
    year               = 2022
  }

  output {
    Array[File] matched_points = jetlag.matched_points
    Array[File] missing_points = jetlag.missing_points
  }
}

task create_manifest {
  meta {
    description: "Create a tile manifest file from a single tile for testing"
  }

  input {
    File tile_file
  }

  command <<<
    # Write the localized absolute path of the tile into a manifest
    realpath "~{tile_file}" > tile_manifest.txt
  >>>

  output {
    File manifest = "tile_manifest.txt"
  }

  runtime {
    docker: "ubuntu:latest"
    memory: "1 GB"
    cpu: 1
  }
}
