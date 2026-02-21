version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jetlag/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jetlag/pipelines/ww-jetlag/ww-jetlag.wdl" as ww_jetlag

workflow jetlag_example {
  # Generate synthetic test tile and border points data
  call ww_testdata.generate_sjl_data { }

  # Run the jetlag pipeline over a small array of test tiles
  call ww_jetlag.jetlag { input:
    tile_paths         = [generate_sjl_data.tile_rds],
    tile_nums          = ["0001"],
    border_points_path = generate_sjl_data.border_points_csv,
    year               = 2022
  }

  output {
    Array[File] matched_points = jetlag.matched_points
    Array[File] missing_points = jetlag.missing_points
  }
}
