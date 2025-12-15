version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-shapemapper/modules/ww-shapemapper/ww-shapemapper.wdl" as ww_shapemapper
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-shapemapper/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct ShapeMapperSample {
    String name
    File target_fa
    File modified_r1
    File modified_r2
    File untreated_r1
    File untreated_r2
}

workflow shapemapper_example {
  # Auto-download ShapeMapper example test data (TPP riboswitch)
  call ww_testdata.download_shapemapper_data as download_shapemapper_data { }

  # Create samples array using ShapeMapper example data
  # This uses the official TPP riboswitch example data from the ShapeMapper repository
  Array[ShapeMapperSample] final_samples = [
    {
      "name": "TPP_riboswitch",
      "target_fa": download_shapemapper_data.target_fa,
      "modified_r1": download_shapemapper_data.modified_r1,
      "modified_r2": download_shapemapper_data.modified_r2,
      "untreated_r1": download_shapemapper_data.untreated_r1,
      "untreated_r2": download_shapemapper_data.untreated_r2
    }
  ]

  # Process each sample
  scatter (sample in final_samples) {
    call ww_shapemapper.run_shapemapper { input:
        sample_name = sample.name,
        target_fa = sample.target_fa,
        modified_r1 = sample.modified_r1,
        modified_r2 = sample.modified_r2,
        untreated_r1 = sample.untreated_r1,
        untreated_r2 = sample.untreated_r2,
        min_depth = 1000,
        is_amplicon = true
    }
  }

  output {
    Array[File] shape_files = run_shapemapper.shape_file
    Array[File] log_files = run_shapemapper.log_file
  }
}
