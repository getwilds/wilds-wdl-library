version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl" as ww_sra
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellranger/ww-cellranger.wdl" as ww_cellranger

workflow cellranger_example {
  # Download test GEX reference
  call download_test_reference { }

  # Download SRA data for SRR9134714 (10x scRNA-seq data)
  # Limiting to 50k reads for faster testing
  call ww_sra.fastqdump { input:
    sra_id = "SRR9134714",
    ncpu = 2,
    max_reads = 50000
  }

  # Run cellranger count
  call ww_cellranger.run_count { input:
    gex_fastqs = [fastqdump.r1_end, fastqdump.r2_end],
    ref_gex = download_test_reference.ref_tar,
    sample_id = "SRR9134714",
    create_bam = false,
    cpu_cores = 4,
    memory_gb = 16,
    expect_cells = 100
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

task download_test_reference {
  meta {
    description: "Download a minimal Cell Ranger reference for testing"
    outputs: {
        ref_tar: "Cell Ranger reference transcriptome tarball"
    }
  }

  command <<<
    set -eo pipefail

    # Download a minimal human reference from 10x Genomics
    # Using the GRCh38 reference (2020-A version)
    # Note: This downloads the full reference (~11GB). For production testing,
    # consider hosting a smaller test reference or using a pre-built minimal reference.

    echo "Downloading Cell Ranger reference..."
    curl -L -o refdata-gex-GRCh38-2020-A.tar.gz \
      "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"

    echo "Reference download complete"
  >>>

  output {
    File ref_tar = "refdata-gex-GRCh38-2020-A.tar.gz"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    memory: "4 GB"
    cpu: 2
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

    echo "=== Cell Ranger Count Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt
    echo "Sample ID: ~{sample_id}" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Check results tarball
    echo "--- Results Tarball ---" >> validation_report.txt
    if [[ -f "~{results_tar}" && -s "~{results_tar}" ]]; then
      tar_size=$(stat -c%s "~{results_tar}" 2>/dev/null || stat -f%z "~{results_tar}")
      echo "Results tarball: ~{results_tar} (${tar_size} bytes)" >> validation_report.txt

      # List contents of tarball
      echo "Tarball contents:" >> validation_report.txt
      tar -tzf "~{results_tar}" | head -20 >> validation_report.txt
      echo "..." >> validation_report.txt
    else
      echo "Results tarball: MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Check web summary
    echo "--- Web Summary ---" >> validation_report.txt
    if [[ -f "~{web_summary}" && -s "~{web_summary}" ]]; then
      summary_size=$(stat -c%s "~{web_summary}" 2>/dev/null || stat -f%z "~{web_summary}")
      echo "Web summary: ~{web_summary} (${summary_size} bytes)" >> validation_report.txt
    else
      echo "Web summary: MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Check and parse metrics summary
    echo "--- Metrics Summary ---" >> validation_report.txt
    if [[ -f "~{metrics_summary}" && -s "~{metrics_summary}" ]]; then
      metrics_size=$(stat -c%s "~{metrics_summary}" 2>/dev/null || stat -f%z "~{metrics_summary}")
      echo "Metrics summary: ~{metrics_summary} (${metrics_size} bytes)" >> validation_report.txt
      echo "" >> validation_report.txt
      echo "Key metrics:" >> validation_report.txt
      cat "~{metrics_summary}" >> validation_report.txt
    else
      echo "Metrics summary: MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Final validation summary
    echo "=== Validation Summary ===" >> validation_report.txt
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
