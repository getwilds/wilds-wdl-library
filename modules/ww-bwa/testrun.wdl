version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bwa/ww-bwa.wdl" as ww_bwa

struct BwaSample {
    String name
    File reads
    File? mates
}

workflow bwa_example {
  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_fastq_data { }
  call ww_testdata.interleave_fastq { 
    input:
      r1_fq = download_fastq_data.r1_fastq,
      r2_fq = download_fastq_data.r2_fastq
  }

  # Create test samples array
  Array[BwaSample] final_samples = [
    object {
      name: "demo_sample",
      reads: download_fastq_data.r1_fastq,
      mates: download_fastq_data.r2_fastq
    },
    object {
      name: "demo_sample_interleaved", 
      reads: interleave_fastq.inter_fastq
    }
  ]

  call ww_bwa.bwa_index { input:
      reference_fasta = download_ref_data.fasta,
      cpu_cores = 2,
      memory_gb = 8
  }

  scatter (sample in final_samples) {
    call ww_bwa.bwa_mem { input:
        bwa_genome_tar = bwa_index.bwa_index_tar,
        reference_fasta = download_ref_data.fasta,
        reads = sample.reads,
        mates = sample.mates,
        name = sample.name,
        paired_end = true,
        cpu_cores = 2,
        memory_gb = 8
    }
  }

  call validate_outputs { input:
    bam_files = bwa_mem.sorted_bam,
    bai_files = bwa_mem.sorted_bai
  }

  output {
    Array[File] bwa_bam = bwa_mem.sorted_bam
    Array[File] bwa_bai = bwa_mem.sorted_bai
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that all expected BWA-MEM output files were generated correctly"
    outputs: {
        report: "Validation report summarizing file checks and alignment statistics"
    }
  }

  parameter_meta {
    bam_files: "Array of BAM files to validate"
    bai_files: "Array of BAM index files to validate"
  }

  input {
    Array[File] bam_files
    Array[File] bai_files
  }

  command <<<
    set -eo pipefail

    echo "=== BWA-MEM Alignment Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    # Arrays for bash processing
    bam_files=(~{sep=" " bam_files})
    bai_files=(~{sep=" " bai_files})

    validation_passed=true
    total_mapped_reads=0

    # Check each sample
    for i in "${!bam_files[@]}"; do
      bam_file="${bam_files[$i]}"
      bai_file="${bai_files[$i]}"

      echo "--- Sample: $bam_file ---" >> validation_report.txt

      # Check BAM file exists and is not empty
      if [[ -f "$bam_file" && -s "$bam_file" ]]; then
        bam_size=$(stat -c%s "$bam_file")
        echo "BAM file: $bam_file (${bam_size} bytes)" >> validation_report.txt

        # Try to get alignment stats from samtools if available
        if command -v samtools &> /dev/null; then
          mapped_reads=$(samtools view -c -F 4 "$bam_file" 2>/dev/null || echo "N/A")
          total_reads=$(samtools view -c "$bam_file" 2>/dev/null || echo "N/A")
          echo "  Total reads: $total_reads" >> validation_report.txt
          echo "  Mapped reads: $mapped_reads" >> validation_report.txt

          if [[ "$mapped_reads" =~ ^[0-9]+$ ]]; then
            total_mapped_reads=$((total_mapped_reads + mapped_reads))
          fi
        fi
      else
        echo "BAM file: $bam_file - MISSING OR EMPTY" >> validation_report.txt
        validation_passed=false
      fi

      # Check BAI file exists
      if [[ -f "$bai_file" ]]; then
        bai_size=$(stat -c%s "$bai_file")
        echo "BAI file: $bai_file (${bai_size} bytes)" >> validation_report.txt
      else
        echo "BAI file: $bai_file - MISSING" >> validation_report.txt
        validation_passed=false
      fi

      echo "" >> validation_report.txt
    done

    # Overall summary
    echo "=== Validation Summary ===" >> validation_report.txt
    echo "Total samples processed: ${#bam_files[@]}" >> validation_report.txt
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
    docker: "getwilds/bwa:0.7.17"
    cpu: 1
    memory: "2 GB"
  }
}
