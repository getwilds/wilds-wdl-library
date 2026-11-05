version 1.0

import "ww-geolink.wdl" as ww_geolink
import "../ww-testdata/ww-testdata.wdl" as ww_testdata

workflow geolink_example {
  input {
    Int n_patients = 5
    Int n_rasters = 3
  }

  # Buffer distances to test: 0.0 exercises point extraction, 1000.0 exercises buffered extraction
  Array[Float] buffer_distances = [0.0, 1000.0]

  # Generate synthetic test data
  call ww_testdata.generate_test_geodata { input:
      n_patients = n_patients,
      n_rasters = n_rasters
  }

  # Test extraction for each raster and buffer distance
  scatter (raster in generate_test_geodata.test_rasters) {
    scatter (dist in buffer_distances) {
      call ww_geolink.extract_raster { input:
          points_csv = generate_test_geodata.patient_csv,
          raster_tif = raster,
          buffer_dist = dist
      }
    }
  }

  Array[File] all_extracts = flatten(extract_raster.result_csv)

  # Test wildfire metrics computation
  call ww_geolink.compute_wildfire_metrics { input:
      daily_extracts = all_extracts,
      p90_pctile = 0.90,
      min_days = 2
  }

  # Test joining everything together
  call ww_geolink.join_exposures { input:
      patient_base = generate_test_geodata.patient_csv,
      all_extracts = all_extracts,
      annual_wf_metrics = compute_wildfire_metrics.annual_metrics
  }

  output {
    File final_dataset = join_exposures.final_dataset
  }
}
