version 1.0

import "ww-geolink.wdl"
import "../ww-testdata/ww-testdata.wdl"

workflow geolink_example {
  input {
    Int n_patients = 5
    Int n_rasters = 3
  }

  # 1. Generate synthetic test data
  var test_data = generate_test_geodata {
    n_patients = n_patients,
    n_rasters = n_rasters
  }

  # 2. Test point extraction for each raster
  scatter (raster in test_data.test_rasters) {
    var extraction = extract_raster {
      points_csv = test_data.patient_csv,
      raster_tif = raster,
      buffer_dist = 0.0,
      stat_fun = "mean",
      output_prefix = "test_ext_"
    }
  }

  # 3. Test wildfire metrics computation (on a subset of extracts)
  var metrics = compute_wildfire_metrics {
    daily_extracts = extraction.result_csv,
    p90_pctile = 0.90,
    min_days = 2
  }

  # 4. Test joining everything together
  var final = join_exposures {
    patient_base = test_data.patient_csv,
    all_extracts = extraction.result_csv,
    annual_wf_metrics = metrics.annual_metrics
  }

  output {
    File final_dataset = final.final_dataset
  }
}
