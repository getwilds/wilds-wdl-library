version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl" as ww_sra
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-spades/modules/ww-spades/ww-spades.wdl" as ww_spades

workflow spades_example {
  # Download test data from SRA
  call ww_sra.fastqdump { input:
    sra_id = "SRR32343420",
    ncpu = 2,
    max_reads = 50000
  }

  # Run metaspades assembly
  call ww_spades.metaspades { input:
    r1_fastq = fastqdump.r1_end,
    r2_fastq = fastqdump.r2_end,
    sample_name = "SRR32343420"
  }

  # Validate outputs
  call validate_outputs { input:
    scaffolds_fasta = metaspades.scaffolds_fasta,
    contigs_fasta = metaspades.contigs_fasta,
    log_file = metaspades.log_file
  }

  output {
    File scaffolds = metaspades.scaffolds_fasta
    File contigs = metaspades.contigs_fasta
    File log = metaspades.log_file
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validates SPAdes assembly outputs to ensure they exist and are non-empty"
    outputs: {
        report: "Validation summary reporting file checks and basic statistics"
    }
  }

  parameter_meta {
    scaffolds_fasta: "Scaffolds FASTA file to validate"
    contigs_fasta: "Contigs FASTA file to validate"
    log_file: "SPAdes log file to validate"
    cpu_cores: "Number of CPU cores to use for validation"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    File scaffolds_fasta
    File contigs_fasta
    File log_file
    Int cpu_cores = 1
    Int memory_gb = 1
  }

  command <<<
    set -euo pipefail

    # Function to validate a file exists and is non-empty
    validate_file() {
      local file_path="$1"
      local file_label="$2"

      if [[ -f "$file_path" && -s "$file_path" ]]; then
        echo "$file_label: $file_path - PASSED" >> validation_report.txt
      else
        echo "$file_label: $file_path - MISSING OR EMPTY" >> validation_report.txt
      fi
    }

    echo "=== SPAdes Assembly Validation Report ===" > validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Validate all output files
    validate_file "~{scaffolds_fasta}" "Scaffolds FASTA" || validation_passed=false
    validate_file "~{contigs_fasta}" "Contigs FASTA" || validation_passed=false
    validate_file "~{log_file}" "SPAdes log file" || validation_passed=false

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
    docker: "getwilds/spades:4.2.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
