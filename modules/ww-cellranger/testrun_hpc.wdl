version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/split-cicd-hpc-testruns/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/split-cicd-hpc-testruns/modules/ww-cellranger/ww-cellranger.wdl" as ww_cellranger

#### TEST WORKFLOW DEFINITION ####
# HPC variant of the Cell Ranger testrun. Calls run_count_hpc_sprocket
# because the monthly HPC test run executes via Sprocket; the Cromwell
# variant (run_count_hpc_cromwell) is exercised by users running Cromwell
# on their HPC by hand. Uses the same tiny test data as testrun.wdl so
# the only thing this run validates is the module-load path.

workflow cellranger_example {

  call ww_testdata.download_test_cellranger_ref { }

  call ww_testdata.download_fastq_data { input:
    prefix = "testdata",
    gzip_output = true
  }

  call ww_cellranger.run_count_hpc_sprocket { input:
    r1_fastqs = [download_fastq_data.r1_fastq],
    r2_fastqs = [download_fastq_data.r2_fastq],
    ref_gex = download_test_cellranger_ref.ref_tar,
    sample_id = "testdata",
    create_bam = false,
    cpu_cores = 2,
    memory_gb = 6,
    chemistry = "SC3Pv2"
  }

  call validate_outputs { input:
    sample_id = "testdata",
    results_tar = run_count_hpc_sprocket.results_tar,
    web_summary = run_count_hpc_sprocket.web_summary,
    metrics_summary = run_count_hpc_sprocket.metrics_summary
  }

  output {
    File results_tar = run_count_hpc_sprocket.results_tar
    File web_summary = run_count_hpc_sprocket.web_summary
    File metrics_summary = run_count_hpc_sprocket.metrics_summary
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
    results_tar: "Compressed tarball of Cell Ranger count output directory"
    web_summary: "Web summary HTML file"
    metrics_summary: "Metrics summary CSV file"
  }

  input {
    String sample_id
    File results_tar
    File web_summary
    File metrics_summary
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

    echo "=== Cell Ranger Count Validation Report (HPC) ===" > validation_report.txt
    echo "" >> validation_report.txt
    echo "Sample ID: ~{sample_id}" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    validate_file "~{results_tar}" "Results tarball" || validation_passed=false
    validate_file "~{web_summary}" "Web summary" || validation_passed=false
    validate_file "~{metrics_summary}" "Metrics summary" || validation_passed=false

    {
      echo ""
      echo "=== Validation Summary ==="
      echo "Total files validated: 3"
    } >> validation_report.txt

    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> validation_report.txt
    else
      echo "Overall Status: FAILED" >> validation_report.txt
      exit 1
    fi

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
