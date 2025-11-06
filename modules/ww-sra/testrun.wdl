version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-star-deseq/modules/ww-sra/ww-sra.wdl" as ww_sra

workflow sra_example {
  # Pulling down test SRA data
  scatter (id in ["ERR1258306"]) {
    call ww_sra.fastqdump { input:
        sra_id = id,
        ncpu = 2,
        max_reads = 1000
    }
  }

  call validate_outputs { input:
    sra_ids = ["ERR1258306"],
    r1_files = fastqdump.r1_end,
    r2_files = fastqdump.r2_end,
    is_paired_flags = fastqdump.is_paired_end
  }

  output {
    Array[File] r1_fastqs = fastqdump.r1_end
    Array[File] r2_fastqs = fastqdump.r2_end
    Array[Boolean] is_paired_end = fastqdump.is_paired_end
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
    r1_files: "Array of R1 FASTQ files to validate"
    r2_files: "Array of R2 FASTQ files to validate"
    sra_ids: "List of SRA IDs that were processed"
    is_paired_flags: "Array of paired-end flags for each sample"
  }

  input {
    Array[File] r1_files
    Array[File] r2_files
    Array[String] sra_ids
    Array[Boolean] is_paired_flags
  }

  command <<<
    set -eo pipefail
    
    echo "=== SRA Download Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt
    
    # Arrays for bash processing
    sra_ids=~{sep=" " sra_ids}
    r1_files=~{sep=" " r1_files}
    r2_files=~{sep=" " r2_files}
    is_paired=~{sep=" " is_paired_flags}
    
    validation_passed=true
    
    # Check each sample
    for i in "${!sra_ids[@]}"; do
      sra_id="${sra_ids[$i]}"
      r1_file="${r1_files[$i]}"
      r2_file="${r2_files[$i]}"
      paired="${is_paired[$i]}"
      
      echo "--- Sample: $sra_id ---" >> validation_report.txt
      echo "Paired-end: $paired" >> validation_report.txt
      
      # Check R1 file exists and is not empty
      if [[ -f "$r1_file" && -s "$r1_file" ]]; then
        r1_size=$(stat -c%s "$r1_file")
        echo "R1 file: $r1_file (${r1_size} bytes)" >> validation_report.txt
      else
        echo "R1 file: $r1_file - MISSING OR EMPTY" >> validation_report.txt
        validation_passed=false
      fi
      
      # Check R2 file - should exist but may be empty for single-end
      if [[ -f "$r2_file" ]]; then
        r2_size=$(stat -c%s "$r2_file")
        if [[ "$paired" == "true" && "$r2_size" -gt 0 ]]; then
          echo "R2 file: $r2_file (${r2_size} bytes)" >> validation_report.txt
        elif [[ "$paired" == "false" && "$r2_size" -eq 0 ]]; then
          echo "R2 file: $r2_file (empty placeholder)" >> validation_report.txt
        else
          echo "R2 file: $r2_file - SIZE MISMATCH FOR PAIRING" >> validation_report.txt
          validation_passed=false
        fi
      else
        echo "R2 file: $r2_file - MISSING" >> validation_report.txt
        validation_passed=false
      fi
    done
      
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
    docker: "getwilds/sra-tools:3.1.1"
    memory: "1 GB"
    cpu: 1
  }
}
