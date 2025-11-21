version 1.0

# Import module in question as well as the testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-salmon/ww-salmon.wdl" as ww_salmon
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

# Define data structure for sample inputs
struct SalmonSample {
    String name
    File r1_fastq
    File r2_fastq
}

#### TEST WORKFLOW DEFINITION ####
# Define test workflow to demonstrate module functionality

workflow salmon_example {
  # Auto-download test data for testing purposes
  call ww_testdata.download_fastq_data as download_demo_data { }
  call ww_testdata.download_test_transcriptome as download_transcriptome { }

  # Build Salmon index
  call ww_salmon.build_index { input:
      transcriptome_fasta = download_transcriptome.transcriptome_fasta,
      cpu_cores = 2,
      memory_gb = 8
  }

  # Create samples array using test data
  Array[SalmonSample] final_samples = [
    {
      "name": "demo_sample",
      "r1_fastq": download_demo_data.r1_fastq,
      "r2_fastq": download_demo_data.r2_fastq
    }
  ]

  # Quantify each sample
  scatter (sample in final_samples) {
    call ww_salmon.quantify { input:
        salmon_index_dir = build_index.salmon_index,
        sample_name = sample.name,
        fastq_r1 = sample.r1_fastq,
        fastq_r2 = sample.r2_fastq,
        cpu_cores = 2,
        memory_gb = 8
    }
  }

  # Merge results
  call ww_salmon.merge_results { input:
      salmon_quant_dirs = quantify.salmon_quant_dir,
      sample_names = quantify.output_sample_name,
      cpu_cores = 1,
      memory_gb = 4
  }

  # Validate outputs
  call validate_outputs { input:
      salmon_quant_dirs = quantify.salmon_quant_dir,
      tpm_matrix = merge_results.tpm_matrix,
      counts_matrix = merge_results.counts_matrix
  }

  output {
    File salmon_index_tar = build_index.salmon_index
    Array[File] salmon_quant_dirs = quantify.salmon_quant_dir
    File merged_tpm_matrix = merge_results.tpm_matrix
    File merged_counts_matrix = merge_results.counts_matrix
    File sample_list = merge_results.sample_list
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that all expected Salmon output files were generated correctly"
    outputs: {
        report: "Validation report summarizing file checks and quantification statistics"
    }
  }

  parameter_meta {
    salmon_quant_dirs: "Array of Salmon quantification directories to validate"
    tpm_matrix: "Merged TPM matrix file to validate"
    counts_matrix: "Merged counts matrix file to validate"
  }

  input {
    Array[File] salmon_quant_dirs
    File tpm_matrix
    File counts_matrix
  }

  command <<<
    set -eo pipefail

    echo "=== Salmon Quantification Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Check each quantification directory
    echo "--- Quantification Directories ---" >> validation_report.txt
    for quant_dir in ~{sep=" " salmon_quant_dirs}; do
      if [[ -f "$quant_dir" && -s "$quant_dir" ]]; then
        quant_size=$(stat -c%s "$quant_dir" 2>/dev/null || stat -f%z "$quant_dir" 2>/dev/null)
        echo "Quant file: $quant_dir (${quant_size} bytes)" >> validation_report.txt
      else
        echo "Quant file: $quant_dir - MISSING OR EMPTY" >> validation_report.txt
        validation_passed=false
      fi
    done
    echo "" >> validation_report.txt

    # Check TPM matrix
    echo "--- TPM Matrix ---" >> validation_report.txt
    if [[ -f "~{tpm_matrix}" && -s "~{tpm_matrix}" ]]; then
      tpm_size=$(stat -c%s "~{tpm_matrix}" 2>/dev/null || stat -f%z "~{tpm_matrix}" 2>/dev/null)
      tpm_lines=$(wc -l < "~{tpm_matrix}" 2>/dev/null || echo "N/A")
      echo "TPM matrix: ~{tpm_matrix} (${tpm_size} bytes, ${tpm_lines} lines)" >> validation_report.txt

      # Show first few lines
      echo "First 5 lines:" >> validation_report.txt
      head -5 "~{tpm_matrix}" >> validation_report.txt
    else
      echo "TPM matrix: ~{tpm_matrix} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Check counts matrix
    echo "--- Counts Matrix ---" >> validation_report.txt
    if [[ -f "~{counts_matrix}" && -s "~{counts_matrix}" ]]; then
      counts_size=$(stat -c%s "~{counts_matrix}" 2>/dev/null || stat -f%z "~{counts_matrix}" 2>/dev/null)
      counts_lines=$(wc -l < "~{counts_matrix}" 2>/dev/null || echo "N/A")
      echo "Counts matrix: ~{counts_matrix} (${counts_size} bytes, ${counts_lines} lines)" >> validation_report.txt

      # Show first few lines
      echo "First 5 lines:" >> validation_report.txt
      head -5 "~{counts_matrix}" >> validation_report.txt
    else
      echo "Counts matrix: ~{counts_matrix} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Overall summary
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
    docker: "combinelab/salmon:latest"
    memory: "2 GB"
    cpu: 1
  }
}
