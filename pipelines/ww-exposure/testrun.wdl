version 1.0

import "ww-exposure.wdl" as ww_exposure
import "../../modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow exposure_example {
  input {
    Int n_patients = 5
    Int n_rasters = 3
  }

  # Generate synthetic test data
  call ww_testdata.generate_test_geodata { input:
      n_patients = n_patients,
      n_rasters = n_rasters
  }

  # Run the exposure linkage pipeline
  call ww_exposure.exposure_linkage { input:
      patient_csv = generate_test_geodata.patient_csv,
      rasters_local = generate_test_geodata.test_rasters,
      buffer_distances = [0.0, 1000.0]
  }

  output {
    File final_dataset = exposure_linkage.final_dataset
  }
}
