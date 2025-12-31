version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/cellranger-r1-r2/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/cellranger-r1-r2/modules/ww-sra/ww-sra.wdl" as ww_sra
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/cellranger-r1-r2/modules/ww-cellranger/ww-cellranger.wdl" as ww_cellranger

workflow cellranger_example {

  # Download a small GEX reference
  call ww_testdata.download_test_cellranger_ref { }

  # Download input FASTQs
  call ww_sra.fastqdump { input:
    sra_id = "SRR9134714",
    ncpu = 2,
    max_reads = 50000
  }

  # Prepare FASTQs (rename to Cell Ranger convention for this testrun WDL)
  call ww_cellranger.prepare_fastqs { input:
    r1_fastqs = [fastqdump.r1_end],
    r2_fastqs = [fastqdump.r2_end],
    sample_name = "SRR9134714"
  }

  # Run cellranger count
  call ww_cellranger.run_count { input:
    r1_fastqs = prepare_fastqs.renamed_r1_fastqs,
    r2_fastqs = prepare_fastqs.renamed_r2_fastqs,
    ref_gex = download_test_cellranger_ref.ref_tar,
    sample_id = "SRR9134714",
    create_bam = false,
    cpu_cores = 2,
    memory_gb = 6,
    chemistry = "SC3Pv2"
  }

  # Validate outputs
  call validate_outputs { input:
    sample_id = "SRR9134714",
    results_tar = run_count.results_tar,
    web_summary = run_count.web_summary,
    metrics_summary = run_count.metrics_summary
  }

  output {
    File results_tar = run_count.results_tar
    File web_summary = run_count.web_summary
    File metrics_summary = run_count.metrics_summary
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

    echo "=== Cell Ranger Count Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt
    echo "Sample ID: ~{sample_id}" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Validate all files using the function
    validate_file "~{results_tar}" "Results tarball" || validation_passed=false
    validate_file "~{web_summary}" "Web summary" || validation_passed=false
    validate_file "~{metrics_summary}" "Metrics summary" || validation_passed=false

    # Final validation summary
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
