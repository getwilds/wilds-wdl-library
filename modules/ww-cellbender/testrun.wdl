version 1.0

# Import module under test and testdata module using relative paths for local linting/testing
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/organize-sra-cellranger-outputs/modules/ww-cellbender/ww-cellbender.wdl" as ww_cellbender
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/organize-sra-cellranger-outputs/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####

workflow cellbender_example {
  # Download a 10x raw feature-barcode matrix for testing
  call ww_testdata.download_10x_raw_h5_data as download_raw_h5 { }

  # Run CellBender remove-background on the test sample (CPU-only for CI).
  # total_droplets_included caps the barcode count to keep runtime manageable.
  call ww_cellbender.remove_background { input:
    input_h5 = download_raw_h5.raw_h5_matrix,
    sample_name = "pbmc_10k_v3",
    expected_cells = 500,
    total_droplets_included = 5000,
    epochs = 10,
    gpu_enabled = false,
    cpu_cores = 4,
    memory_gb = 8
  }

  output {
    File output_h5         = remove_background.output_h5
    File filtered_h5       = remove_background.filtered_h5
    File report_html       = remove_background.report_html
    File summary_pdf       = remove_background.summary_pdf
    File log_file          = remove_background.log_file
    File cell_barcodes_csv = remove_background.cell_barcodes_csv
    File metrics_csv       = remove_background.metrics_csv
    File checkpoint_tar    = remove_background.checkpoint_tar
  }
}
