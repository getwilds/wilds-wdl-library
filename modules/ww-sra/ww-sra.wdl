## Pulls down paired fastq's for SRA ID's provided

version 1.0

# TODO: add metadata pull using efetch
# TODO: add checksum for validation purposes
# TODO: add logging for better visibility

#### WORKFLOW DEFINITION

workflow sra_example {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow to download raw sequencing data from SRA in parallel"
    url: "https://github.com/getwilds/ww-sra"
    outputs: {
        r1_fastqs: "array of R1 fastq files for each sample",
        r2_fastqs: "array of R2 fastq files for each sample",
        is_paired_end: "array of booleans indicating whether each sample used paired-end sequencing",
        validation_report: "validation report confirming all expected outputs were generated"
    }
  }

  scatter (id in ["ERR1258306"]) {
    call fastqdump { input:
        sra_id = id,
        ncpu = 2,
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

#### TASK DEFINITIONS

task fastqdump {
  meta {
    description: "Task for pulling down fastq data from SRA."
    outputs: {
        r1_end: "R1 fastq file downloaded for the sample in question",
        r2_end: "R2 fastq file downloaded for the sample in question (empty file for single-end reads)",
        is_paired_end: "boolean indicating whether the sample used paired-end sequencing"
    }
  }

  parameter_meta {
    sra_id: "SRA ID of the sample to be downloaded via parallel-fastq-dump"
    ncpu: "number of cpus to use during download"
  }

  input {
    String sra_id
    Int ncpu = 8
  }

  command <<<
    set -eo pipefail
    # check if paired ended
    numLines=$(fastq-dump -X 1 -Z --split-spot "~{sra_id}" | wc -l)
    if [ "$numLines" -eq 8 ]; then
      echo true > paired_file
      parallel-fastq-dump \
        --sra-id "~{sra_id}" \
        --threads ~{ncpu} \
        --outdir ./ \
        --split-files \
        --gzip
    else
      echo false > paired_file
      parallel-fastq-dump \
        --sra-id "~{sra_id}" \
        --threads ~{ncpu} \
        --outdir ./ \
        --gzip
      # Rename the file to match the expected output format
      mv "~{sra_id}.fastq.gz" "~{sra_id}_1.fastq.gz"
      # Create an empty placeholder for R2
      touch "~{sra_id}_2.fastq.gz"
    fi
  >>>

  output {
    File r1_end = "~{sra_id}_1.fastq.gz"
    File r2_end = "~{sra_id}_2.fastq.gz"
    Boolean is_paired_end = read_boolean("paired_file")
  }

  runtime {
    memory: 2 * ncpu + " GB"
    docker: "getwilds/sra-tools:3.1.1"
    cpu: ncpu
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
