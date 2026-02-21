version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jetlag/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jetlag/pipelines/ww-jetlag/ww-jetlag.wdl" as ww_jetlag

workflow jetlag_example {
  # Download test tile and border points data
  call ww_testdata.download_sjl_data { }

  # Run the jetlag pipeline over a small array of test tiles
  call ww_jetlag.jetlag { input:
    tile_paths         = download_sjl_data.tile_rds_array,
    tile_nums          = download_sjl_data.tile_num_array,
    border_points_path = download_sjl_data.border_points_csv,
    year               = 2022
  }

  output {
    Array[File] matched_points = jetlag.matched_points
    Array[File] missing_points = jetlag.missing_points
  }
}
