version 1.0

import "ww-exposure.wdl"
import "../../modules/ww-testdata/ww-testdata.wdl"

workflow exposure_example {
  input {
    Int n_patients = 5
    Int n_rasters = 2
  }

  # 1. Generate synthetic test data
  var test_data = generate_test_geodata {
    n_patients = n_patients,
    n_rasters = n_rasters
  }

  # 2. Run the exposure linkage pipeline
  var pipeline = exposure_linkage {
    patient_csv = test_data.patient_csv,
    rasters_local = test_data.test_rasters,
    buffer_distances = [0.0, 1000.0]
  }

  output {
    File final_dataset = pipeline.final_dataset
  }
}
