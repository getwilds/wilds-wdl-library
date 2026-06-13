version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/run-selector-input/modules/ww-cellranger/ww-cellranger.wdl" as ww_cellranger

workflow cellranger_example {

  # Download a small GEX reference
  call ww_testdata.download_test_cellranger_ref { }

  # Download input FASTQs with Cell Ranger naming convention
  call ww_testdata.download_fastq_data { input:
    prefix = "testdata",
    gzip_output = true
  }

  # Run cellranger count
  call ww_cellranger.run_count { input:
    r1_fastqs = [download_fastq_data.r1_fastq],
    r2_fastqs = [download_fastq_data.r2_fastq],
    ref_gex = download_test_cellranger_ref.ref_tar,
    sample_id = "testdata",
    create_bam = false,
    cpu_cores = 2,
    memory_gb = 6,
    chemistry = "SC3Pv2"
  }

  # Unwrap the Optional outputs from run_count. With single-cell test
  # data and an explicit chemistry, the success branch always populates
  # these — select_first turns File? back into File.
  File results_tar_unwrapped = select_first([run_count.results_tar])
  File web_summary_unwrapped = select_first([run_count.web_summary])
  File metrics_summary_unwrapped = select_first([run_count.metrics_summary])
  File filtered_h5_unwrapped = select_first([run_count.filtered_h5])

  # Validate outputs
  call validate_outputs { input:
    sample_id = "testdata",
    chemistry_status = run_count.chemistry_status,
    results_tar = results_tar_unwrapped,
    web_summary = web_summary_unwrapped,
    metrics_summary = metrics_summary_unwrapped,
    filtered_h5 = filtered_h5_unwrapped
  }

  output {
    File chemistry_status = run_count.chemistry_status
    File results_tar = results_tar_unwrapped
    File web_summary = web_summary_unwrapped
    File metrics_summary = metrics_summary_unwrapped
    File filtered_h5 = filtered_h5_unwrapped
    File validation_report = validate_outputs.report
  }
}


task validate_outputs {
  meta {
    description: "Validate that all expected output files were generated correctly"
    outputs: {
        report: "Validation report summarizing file checks and statistics"
    }
  }

  parameter_meta {
    sample_id: "Sample ID used for the analysis"
    chemistry_status: "Chemistry-status marker file from run_count"
    results_tar: "Compressed tarball of Cell Ranger count output directory"
    web_summary: "Web summary HTML file"
    metrics_summary: "Metrics summary CSV file"
    filtered_h5: "Filtered feature-barcode matrix HDF5 file"
  }

  input {
    String sample_id
    File chemistry_status
    File results_tar
    File web_summary
    File metrics_summary
    File filtered_h5
  }

  command <<<
    set -eo pipefail

    validate_file() {
      local file_path="$1"
      local file_label="$2"

      if [[ -f "$file_path" && -s "$file_path" ]]; then
        echo "$file_label: $file_path - PASSED" >> validation_report.txt
        return 0
      else
        echo "$file_label: $file_path - MISSING OR EMPTY" >> validation_report.txt
        return 1
      fi
    }

    echo "=== Cell Ranger Count Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt
    echo "Sample ID: ~{sample_id}" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Validate all files using the function
    validate_file "~{chemistry_status}" "Chemistry status" || validation_passed=false
    validate_file "~{results_tar}" "Results tarball" || validation_passed=false
    validate_file "~{web_summary}" "Web summary" || validation_passed=false
    validate_file "~{metrics_summary}" "Metrics summary" || validation_passed=false
    validate_file "~{filtered_h5}" "Filtered h5" || validation_passed=false

    # The test fixture is real single-cell data, so chemistry_status
    # must be "ok" — anything else means run_count silently skipped.
    if [[ "$(cat "~{chemistry_status}")" != "ok" ]]; then
      echo "Chemistry status check: expected 'ok', got '$(cat "~{chemistry_status}")' - FAILED" >> validation_report.txt
      validation_passed=false
    fi

    # Final validation summary
    {
      echo ""
      echo "=== Validation Summary ==="
      echo "Total files validated: 5"
    } >> validation_report.txt

    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> validation_report.txt
    else
      echo "Overall Status: FAILED" >> validation_report.txt
      exit 1
    fi

    # Also output to stdout for immediate feedback
    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    memory: "2 GB"
    cpu: 1
  }
}
