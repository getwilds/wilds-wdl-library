version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/organize-sra-cellranger-outputs/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/organize-sra-cellranger-outputs/pipelines/ww-sra-cellranger/ww-sra-cellranger.wdl" as sra_cellranger_workflow

workflow sra_cellranger_example {
  # Download a small GEX reference for Cell Ranger
  call ww_testdata.download_test_cellranger_ref { }

  # Call the workflow with two test samples:
  #   - SRR7722937: 10x's own "1k PBMCs from a Healthy Donor"
  #     demo dataset (should run successfully).
  #   - SRR1039508: Bulk RNA-seq (should fail and appear in
  #     skipped_sample_list).
  call sra_cellranger_workflow.sra_cellranger { input:
    sra_id_list = ["SRR7722937", "SRR1039508"],
    ref_gex = download_test_cellranger_ref.ref_tar,
    ncpu = 2,
    memory_gb = 6,
    max_reads = 5000000,
    skip_on_chemistry_failure = true,
    organize_results = true,
    cellbender_gpu_enabled = false,
    cellbender_memory_gb = 8
  }

  Int n_results = length(sra_cellranger.cellranger_results)

  call validate_outputs { input:
    single_cell_sample_list = sra_cellranger.single_cell_sample_list,
    skipped_sample_list = sra_cellranger.skipped_sample_list,
    n_results = n_results,
    expected_single_cell_id = "SRR7722937",
    expected_skipped_id = "SRR1039508",
    organized_results = sra_cellranger.organized_results,
    cellbender_output_h5s = sra_cellranger.cellbender_output_h5s,
    cellbender_filtered_h5s = sra_cellranger.cellbender_filtered_h5s
  }

  output {
    File single_cell_sample_list = sra_cellranger.single_cell_sample_list
    File skipped_sample_list = sra_cellranger.skipped_sample_list
    Array[File] cellranger_results = sra_cellranger.cellranger_results
    Array[File] cellranger_web_summaries = sra_cellranger.cellranger_web_summaries
    Array[File] cellranger_metrics = sra_cellranger.cellranger_metrics
    Array[File] cellranger_filtered_h5s = sra_cellranger.cellranger_filtered_h5s
    Array[File] cellranger_raw_h5s = sra_cellranger.cellranger_raw_h5s
    Array[File] cellbender_output_h5s = sra_cellranger.cellbender_output_h5s
    Array[File] cellbender_filtered_h5s = sra_cellranger.cellbender_filtered_h5s
    File? organized_results = sra_cellranger.organized_results
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Assert that the skip_on_chemistry_failure=true scatter correctly partitioned the single-cell and non-single-cell samples, and that the organized results ZIP contains the expected per-sample subdirectory."
    outputs: {
        report: "Validation report"
    }
  }

  parameter_meta {
    single_cell_sample_list: "single_cell_sample_list.txt produced by ww-sra-cellranger"
    skipped_sample_list: "skipped_sample_list.txt produced by ww-sra-cellranger"
    n_results: "Length of the filtered cellranger_results array (should be 1)"
    expected_single_cell_id: "Sample ID expected in single_cell_sample_list"
    expected_skipped_id: "Sample ID expected in skipped_sample_list"
    organized_results: "ZIP archive produced by organize_outputs (optional; checked when present)"
    cellbender_output_h5s: "CellBender cleaned H5 files (one per successful sample)"
    cellbender_filtered_h5s: "CellBender filtered H5 files (one per successful sample)"
  }

  input {
    File single_cell_sample_list
    File skipped_sample_list
    Int n_results
    String expected_single_cell_id
    String expected_skipped_id
    File? organized_results
    Array[File] cellbender_output_h5s
    Array[File] cellbender_filtered_h5s
  }

  command <<<
    set -eo pipefail

    echo "=== ww-sra-cellranger skip-behavior validation ===" > validation_report.txt

    expected_single_cell="~{expected_single_cell_id}"
    expected_skipped="~{expected_skipped_id}"

    actual_single_cell=$(tr -d '[:space:]' < "~{single_cell_sample_list}")
    actual_skipped=$(tr -d '[:space:]' < "~{skipped_sample_list}")

    status=PASSED

    if [ "$actual_single_cell" != "$expected_single_cell" ]; then
      echo "single_cell_sample_list mismatch: expected '$expected_single_cell', got '$actual_single_cell'" >> validation_report.txt
      status=FAILED
    else
      echo "single_cell_sample_list = $expected_single_cell - PASSED" >> validation_report.txt
    fi

    if [ "$actual_skipped" != "$expected_skipped" ]; then
      echo "skipped_sample_list mismatch: expected '$expected_skipped', got '$actual_skipped'" >> validation_report.txt
      status=FAILED
    else
      echo "skipped_sample_list = $expected_skipped - PASSED" >> validation_report.txt
    fi

    if [ "~{n_results}" -ne 1 ]; then
      echo "cellranger_results length: expected 1, got ~{n_results}" >> validation_report.txt
      status=FAILED
    else
      echo "cellranger_results length = 1 - PASSED" >> validation_report.txt
    fi

    n_cellbender_output=~{length(cellbender_output_h5s)}
    if [ "$n_cellbender_output" -ne 1 ]; then
      echo "cellbender_output_h5s length: expected 1, got $n_cellbender_output" >> validation_report.txt
      status=FAILED
    else
      echo "cellbender_output_h5s length = 1 - PASSED" >> validation_report.txt
    fi

    n_cellbender_filtered=~{length(cellbender_filtered_h5s)}
    if [ "$n_cellbender_filtered" -ne 1 ]; then
      echo "cellbender_filtered_h5s length: expected 1, got $n_cellbender_filtered" >> validation_report.txt
      status=FAILED
    else
      echo "cellbender_filtered_h5s length = 1 - PASSED" >> validation_report.txt
    fi

    if [ -n "~{organized_results}" ]; then
      tar_listing=$(tar -tf "~{organized_results}")
      echo "=== organized_results tar listing ===" >> validation_report.txt
      echo "$tar_listing" >> validation_report.txt
      echo "=====================================" >> validation_report.txt
      if echo "$tar_listing" | grep -q "/$expected_single_cell/"; then
        echo "organized_results contains $expected_single_cell subdirectory - PASSED" >> validation_report.txt
      else
        echo "organized_results missing $expected_single_cell subdirectory" >> validation_report.txt
        status=FAILED
      fi
      if echo "$tar_listing" | grep -q "/$expected_skipped/"; then
        echo "organized_results contains $expected_skipped subdirectory (skipped sample should be absent) - FAILED" >> validation_report.txt
        status=FAILED
      else
        echo "organized_results does not contain $expected_skipped subdirectory (skipped sample correctly absent) - PASSED" >> validation_report.txt
      fi
    fi

    echo "" >> validation_report.txt
    echo "Overall Status: $status" >> validation_report.txt
    cat validation_report.txt

    if [ "$status" != "PASSED" ]; then
      exit 1
    fi
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: 1
    memory: "2 GB"
  }
}
