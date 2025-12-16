version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-varscan/modules/ww-samtools/ww-samtools.wdl" as ww_samtools

workflow samtools_example {
  # Download test data
  call ww_testdata.download_ref_data { }

  # Download two CRAM files to test merging in crams_to_fastq
  call ww_testdata.download_cram_data as download_cram_1 { input:
      ref_fasta = download_ref_data.fasta
  }
  call ww_testdata.download_cram_data as download_cram_2 { input:
      ref_fasta = download_ref_data.fasta
  }

  # Download two BAM files to test merging in merge_bams_to_cram
  call ww_testdata.download_bam_data as download_bam_1 { }
  call ww_testdata.download_bam_data as download_bam_2 { }

  # Convert multiple CRAMs to FASTQ (tests the merge functionality)
  call ww_samtools.crams_to_fastq { input:
      cram_files = [download_cram_1.cram, download_cram_2.cram],
      ref = download_ref_data.fasta,
      name = "test_sample",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Test merge_bams_to_cram task with multiple BAMs
  call ww_samtools.merge_bams_to_cram { input:
      bams_to_merge = [download_bam_1.bam, download_bam_2.bam],
      base_file_name = "test_merged",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Test mpileup task with a BAM file
  call ww_samtools.mpileup { input:
      bamfile = download_bam_1.bam,
      ref_fasta = download_ref_data.fasta,
      sample_name = "test_sample",
      disable_baq = true,
      cpu_cores = 2,
      memory_gb = 8
  }

  # Validate outputs
  call validate_outputs { input:
      r1_fastq = crams_to_fastq.r1_fastq,
      r2_fastq = crams_to_fastq.r2_fastq,
      merged_cram = merge_bams_to_cram.cram,
      merged_crai = merge_bams_to_cram.crai,
      pileup = mpileup.pileup
  }

  output {
    File r1_fastqs = crams_to_fastq.r1_fastq
    File r2_fastqs = crams_to_fastq.r2_fastq
    File merged_cram = merge_bams_to_cram.cram
    File pileup_file = mpileup.pileup
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate all Samtools output files"
    outputs: {
        report: "Validation report summarizing file checks and statistics"
    }
  }

  parameter_meta {
    r1_fastq: "R1 FASTQ file to validate"
    r2_fastq: "R2 FASTQ file to validate"
    merged_cram: "Merged CRAM file to validate"
    merged_crai: "Merged CRAM index file to validate"
    pileup: "Pileup file to validate"
  }

  input {
    File r1_fastq
    File r2_fastq
    File merged_cram
    File merged_crai
    File pileup
  }

  command <<<
    set -eo pipefail

    echo "=== Samtools Output Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Check R1 FASTQ file
    echo "--- R1 FASTQ File ---" >> validation_report.txt
    # Check if file exists and is non-empty
    if [[ -f "~{r1_fastq}" && -s "~{r1_fastq}" ]]; then
      r1_size=$(wc -c < "~{r1_fastq}")
      echo "R1 FASTQ: ~{r1_fastq} (${r1_size} bytes)" >> validation_report.txt
    else
      echo "R1 FASTQ: ~{r1_fastq} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Check R2 FASTQ file
    echo "--- R2 FASTQ File ---" >> validation_report.txt
    if [[ -f "~{r2_fastq}" && -s "~{r2_fastq}" ]]; then
      r2_size=$(wc -c < "~{r2_fastq}")
      echo "R2 FASTQ: ~{r2_fastq} (${r2_size} bytes)" >> validation_report.txt
    else
      echo "R2 FASTQ: ~{r2_fastq} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Check merged CRAM file
    echo "--- Merged CRAM File ---" >> validation_report.txt
    if [[ -f "~{merged_cram}" && -s "~{merged_cram}" ]]; then
      cram_size=$(wc -c < "~{merged_cram}")
      echo "Merged CRAM: ~{merged_cram} (${cram_size} bytes)" >> validation_report.txt
    else
      echo "Merged CRAM: ~{merged_cram} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Check CRAM index file
    echo "--- CRAM Index File ---" >> validation_report.txt
    if [[ -f "~{merged_crai}" && -s "~{merged_crai}" ]]; then
      crai_size=$(wc -c < "~{merged_crai}")
      echo "CRAM index: ~{merged_crai} (${crai_size} bytes)" >> validation_report.txt
    else
      echo "CRAM index: ~{merged_crai} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Check pileup file
    echo "--- Pileup File ---" >> validation_report.txt
    if [[ -f "~{pileup}" && -s "~{pileup}" ]]; then
      pileup_size=$(wc -c < "~{pileup}")
      echo "Pileup file: ~{pileup} (${pileup_size} bytes)" >> validation_report.txt
    else
      echo "Pileup file: ~{pileup} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Overall summary
    echo "=== Validation Summary ===" >> validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> validation_report.txt
      echo "All Samtools outputs were generated successfully." >> validation_report.txt
    else
      echo "Overall Status: FAILED" >> validation_report.txt
      echo "One or more output files are missing or empty." >> validation_report.txt
      exit 1
    fi

    # Also output to stdout for immediate feedback
    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/samtools:1.19"
    memory: "2 GB"
    cpu: 1
  }
}
