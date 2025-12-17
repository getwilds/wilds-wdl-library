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
    fasta = metaspades.scaffolds_fasta
  }

  output {
    File fasta = metaspades.scaffolds_fasta
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate SPAdes output"
    outputs: {
        report: "Validation report summarizing file checks and statistics"
    }
  }

  parameter_meta {
    fasta: "FASTA file to validate"
  }

  input {
    File fasta
  }

  command <<<
    set -eo pipefail

    echo "=== SPAdes Assembly Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    # Check FASTA file
    echo "--- Scaffolds FASTA File ---" >> validation_report.txt
    if [[ -f "~{fasta}" && -s "~{fasta}" ]]; then
      fasta_size=$(wc -c < "~{fasta}")
      echo "FASTA file: ~{fasta} (${fasta_size} bytes)" >> validation_report.txt
      echo "" >> validation_report.txt
      echo "Overall Status: PASSED" >> validation_report.txt
    else
      echo "FASTA file: ~{fasta} - MISSING OR EMPTY" >> validation_report.txt
      echo "" >> validation_report.txt
      echo "Overall Status: FAILED" >> validation_report.txt
      exit 1
    fi
    echo "" >> validation_report.txt

    # Also output to stdout for immediate feedback
    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "staphb/spades:4.2.0"
    memory: "1 GB"
    cpu: 1
  }
}
