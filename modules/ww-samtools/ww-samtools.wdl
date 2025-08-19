## WILDS WDL for processing genomic files with Samtools.
## Intended for use alone or as a modular component in the WILDS ecosystem.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct SamtoolsSample {
  String name
  Array[File] cram_files
}

workflow samtools_example {
  meta {
    author: "Emma Bishop"
    email: "wilds@fredhutch.org"
    description: "WDL workflow for processing genomic files with Samtools"
    url: "https://github.com/getwilds/wilds-wdl-library/modules/ww-samtools"
    outputs: {
        r1_fastqs: "Array of R1 FASTQ files generated from CRAM/BAM/SAM files",
        r2_fastqs: "Array of R2 FASTQ files generated from CRAM/BAM/SAM files",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    samples: "Array of sample objects, each containing name and paths to CRAM/BAM/SAM files"
    reference_fasta: "Reference genome FASTA file"
    cpus: "Number of CPU cores allocated for each non-validation task"
    memory_gb: "Memory in GB allocated for each non-validation task"
  }

  input {
    Array[SamtoolsSample]? samples
    File? reference_fasta
    Int cpus = 2
    Int memory_gb = 8
  }

  # If no reference genome provided, download test data
  if (!defined(reference_fasta)) {
    call ww_testdata.download_ref_data { }
  }

  # Determine which genome file to use
  File genome_fasta = select_first([reference_fasta, download_ref_data.fasta])

  # If no samples provided, download demonstration data
  if (!defined(samples)) {
    call ww_testdata.download_cram_data { input:
        ref_fasta = genome_fasta
    }
  }

  # Create test sample array when no samples provided
  Array[SamtoolsSample]? test_samples = if defined(download_cram_data.cram) then [
    object {
      name: "test_sample",
      cram_files: [
        select_first([download_cram_data.cram]),
        select_first([download_cram_data.cram]),
        select_first([download_cram_data.cram])
      ]
    }
  ] else []

  # Create the samples array - either from input or from test data
  Array[SamtoolsSample]? final_samples = select_first([samples, test_samples])

  scatter (sample in final_samples) {
    call crams_to_fastq { input:
        cram_files = sample.cram_files,
        ref = genome_fasta,
        name = sample.name,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  call validate_outputs { input:
      r1_fastqs = crams_to_fastq.r1_fastq,
      r2_fastqs = crams_to_fastq.r2_fastq,
      sample_names = crams_to_fastq.sample_name
  }

  output {
    Array[File] r1_fastqs = crams_to_fastq.r1_fastq
    Array[File] r2_fastqs = crams_to_fastq.r2_fastq
    File validation_report = validate_outputs.report
  }
}

task crams_to_fastq {
  meta {
    description: "Merge CRAM/BAM/SAM files and convert to FASTQ's using samtools."
    outputs: {
        r1_fastq: "R1 FASTQ file generated from merged CRAM/BAM/SAM file",
        r2_fastq: "R2 FASTQ file generated from merged CRAM/BAM/SAM file",
        sample_name: "Sample name that was processed"
    }
  }

  parameter_meta {
    cram_files: "List of CRAM/BAM/SAM files for a sample"
    ref: "Reference genome FASTA file"
    name: "Name of the sample (used for file naming)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    Array[File] cram_files
    File ref
    String name
    Int cpu_cores = 2
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    # Merge CRAM/BAM/SAM files if more than one, then collate and convert to FASTQ
    samtools merge -@ ~{cpu_cores} --reference "~{ref}" -u ~{sep=" " cram_files} | \
    samtools collate -@ ~{cpu_cores} --reference "~{ref}" -u -O -T "$TMPDIR" | \
    samtools fastq -@ ~{cpu_cores} --reference "~{ref}" -1 "~{name}_R1.fastq.gz" -2 "~{name}_R2.fastq.gz" -0 /dev/null -s /dev/null -

    # Clean up temp directory
    rm -rf "${TMPDIR:?}"/*
  >>>

  output {
    File r1_fastq = "~{name}_R1.fastq.gz"
    File r2_fastq = "~{name}_R2.fastq.gz"
    String sample_name = "~{name}"
  }

  runtime {
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
    docker: "getwilds/samtools:1.19"
  }
}

task validate_outputs {
  meta {
    description: "Validates that FASTQ files exist and are non-empty after CRAM-to-FASTQ conversion."
    outputs: {
        report: "Validation report with pass/fail summary"
    }
  }

  parameter_meta {
    r1_fastqs: "Array of R1 FASTQ files to check"
    r2_fastqs: "Array of R2 FASTQ files to check"
    sample_names: "Array of sample names corresponding to FASTQ files"
  }

  input {
    Array[File] r1_fastqs
    Array[File] r2_fastqs
    Array[String] sample_names
  }

  command <<<
    set -eo pipefail

    echo "=== Samtools FASTQ Validation Report ===" > samtools_validation_report.txt
    echo "" >> samtools_validation_report.txt

    r1_fastqs=(~{sep=" " r1_fastqs})
    r2_fastqs=(~{sep=" " r2_fastqs})
    sample_names=(~{sep=" " sample_names})

    validation_passed=true

    for i in "${!sample_names[@]}"; do
      sample="${sample_names[$i]}"
      r1="${r1_fastqs[$i]}"
      r2="${r2_fastqs[$i]}"

      echo "--- Sample: $sample ---" >> samtools_validation_report.txt

      if [[ -f "$r1" && -s "$r1" ]]; then
        size=$(stat -c%s "$r1")
        lines=$(zcat "$r1" | wc -l)
        echo "  R1 FASTQ: PASS (${size} bytes, ${lines} lines)" >> samtools_validation_report.txt
      else
        echo "  R1 FASTQ: FAIL - MISSING OR EMPTY" >> samtools_validation_report.txt
        validation_passed=false
      fi

      if [[ -f "$r2" && -s "$r2" ]]; then
        size=$(stat -c%s "$r2")
        lines=$(zcat "$r2" | wc -l)
        echo "  R2 FASTQ: PASS (${size} bytes, ${lines} lines)" >> samtools_validation_report.txt
      else
        echo "  R2 FASTQ: FAIL - MISSING OR EMPTY" >> samtools_validation_report.txt
        validation_passed=false
      fi

      echo "" >> samtools_validation_report.txt
    done

    echo "=== Summary ===" >> samtools_validation_report.txt
    echo "Samples processed: ${#sample_names[@]}" >> samtools_validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
      echo "Status: PASSED" >> samtools_validation_report.txt
    else
      echo "Status: FAILED" >> samtools_validation_report.txt
      exit 1
    fi

    cat samtools_validation_report.txt
  >>>

  output {
    File report = "samtools_validation_report.txt"
  }

  runtime {
    docker: "getwilds/samtools:1.19"
    cpu: 1
    memory: "2 GB"
  }
}
