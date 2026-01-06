version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl" as ww_sra
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-megahit/ww-megahit.wdl" as ww_megahit

workflow megahit_example {
  # Download test data from SRA
  call ww_sra.fastqdump { input:
    sra_id = "SRR32343420",
    ncpu = 2,
    max_reads = 50000
  }

  # Interleave the paired-end reads
  call ww_testdata.interleave_fastq { input:
    r1_fq = fastqdump.r1_end,
    r2_fq = fastqdump.r2_end
  }

  # Run MEGAHIT assembly
  call ww_megahit.megahit { input:
    input_fastq = interleave_fastq.inter_fastq,
    sample_name = "SRR32343420"
  }

  # Validate outputs
  call validate_outputs { input:
    contigs_fasta = megahit.contigs
  }

  output {
    File contigs = megahit.contigs
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validates MEGAHIT assembly outputs to ensure they exist and are non-empty"
    outputs: {
        report: "Validation summary reporting file size check"
    }
  }

  parameter_meta {
    contigs_fasta: "Contigs FASTA file to validate"
    cpu_cores: "Number of CPU cores to use for validation"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    File contigs_fasta
    Int cpu_cores = 1
    Int memory_gb = 1
  }

  command <<<
    set -euo pipefail

    echo "=== MEGAHIT Assembly Validation Report ===" > validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    # Check if file exists and is non-empty
    if [[ -f "~{contigs_fasta}" && -s "~{contigs_fasta}" ]]; then
      echo "Contigs FASTA: ~{contigs_fasta} - PASSED" >> validation_report.txt
      echo "" >> validation_report.txt
      echo "Overall Status: PASSED" >> validation_report.txt
    else
      echo "Contigs FASTA: ~{contigs_fasta} - MISSING OR EMPTY" >> validation_report.txt
      echo "" >> validation_report.txt
      echo "Overall Status: FAILED" >> validation_report.txt
      exit 1
    fi

    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/megahit:1.2.9"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
