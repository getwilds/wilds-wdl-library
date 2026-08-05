version 1.0

# Import module in question as well as the testdata module for automatic demo functionality
import "./ww-dropletutils.wdl" as ww_dropletutils
import "../ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####
# Define test workflow to demonstrate module functionality

workflow dropletutils_example {
  # Auto-download test data for testing purposes
  call ww_testdata.download_10x_h5_data as download_filtered_matrix { }
  call ww_testdata.download_10x_raw_h5_data as download_raw_matrix { }

  # Load a filtered feature-barcode matrix directly into a SingleCellExperiment
  call ww_dropletutils.read10x_counts { input:
      h5_matrix = download_filtered_matrix.h5_matrix,
      sample_name = "demo_filtered_sample"
  }

  # Call real cells from a raw feature-barcode matrix using emptyDrops
  call ww_dropletutils.empty_drops_filter { input:
      raw_h5_matrix = download_raw_matrix.raw_h5_matrix,
      sample_name = "demo_raw_sample"
  }

  output {
    File sce_rds = read10x_counts.sce_rds
    File filtered_sce_rds = empty_drops_filter.filtered_sce_rds
    File empty_drops_csv = empty_drops_filter.empty_drops_csv
    File barcode_rank_pdf = empty_drops_filter.barcode_rank_pdf
  }
}
