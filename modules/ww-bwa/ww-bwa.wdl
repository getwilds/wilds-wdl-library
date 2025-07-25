## WILDS WDL for performing sequence alignment using BWA-MEM.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct BwaSample {
    String name
    File reads
    File? mates
}

workflow bwa_example {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "WDL workflow for sequence alignment via BWA-MEM"
    url: "https://github.com/getwilds/wilds-wdl-library/modules/ww-bwa"
    outputs: {
        bwa_bam: "Sorted BWA-MEM alignment output BAM files for each sample",
        bwa_bai: "Index files for the sorted BWA-MEM alignment BAM files",
        validation_report: "Validation report summarizing file checks and alignment statistics"
    }
  }

  parameter_meta {
    reference_fasta: "Optional reference genome FASTA file. If not provided, test data will be used."
    samples: "Optional list of BwaSample objects, each containing sample name and FASTQ file(s). If not provided, workflow will download a demonstration sample from SRA"
    paired: "Optional boolean indicating if reads are paired end (default: true)"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    memory_gb: "Memory allocated for each task in the workflow in GB"
  }

  input {
    File? reference_fasta
    Array[BwaSample]? samples
    Boolean paired = true
    Int cpus = 2
    Int memory_gb = 8
  }

   # If no reference genome provided, download test data
  if (!defined(reference_fasta)) {
    call ww_testdata.download_ref_data { }
  }

   # Determine which genome files to use
  File genome_fasta = select_first([reference_fasta, download_ref_data.fasta])

  # Handle samples - always download test data, but only use if no samples provided
  call ww_testdata.download_fastq_data { }
  call ww_testdata.interleave_fastq { 
    input:
      r1_fq = download_fastq_data.r1_fastq,
      r2_fq = download_fastq_data.r2_fastq
  }

  # Create default samples from test data
  Array[BwaSample] default_samples = [
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

  # Use provided samples or default to test samples
  Array[BwaSample] final_samples = select_first([samples, default_samples])

  call bwa_index { input:
      reference_fasta = genome_fasta,
      cpu_cores = cpus,
      memory_gb = memory_gb
  }

  scatter (sample in final_samples) {
    call bwa_mem { input:
        bwa_genome_tar = bwa_index.bwa_index_tar,
        reference_fasta = genome_fasta,
        reads = sample.reads,
        mates = sample.mates,
        name = sample.name,
        paired_end = paired,
        cpu_cores = cpus,
        memory_gb = memory_gb
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

task bwa_index {
  meta {
    description: "Task for building BWA index files from a reference FASTA"
    outputs: {
        bwa_index_tar: "Compressed tarball containing BWA genome index"
    }
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File reference_fasta
    Int cpu_cores = 8
    Int memory_gb = 32
  }

  String ref_name = basename(reference_fasta)  # Name of local copy

  command <<<
    set -eo pipefail && \
    mkdir -p bwa_index && \
    cp "~{reference_fasta}" bwa_index/"~{ref_name}" && \
    echo "Building BWA index..." && \
    bwa index "bwa_index/~{ref_name}" && \
    tar -czf bwa_index.tar.gz bwa_index/*
  >>>

  output {
    File bwa_index_tar = "bwa_index.tar.gz"
  }

  runtime {
    docker: "getwilds/bwa:0.7.17"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task bwa_mem {
  meta {
    description: "Task for aligning sequence reads using BWA-MEM"
    outputs: {
        sorted_bam: "Sorted BWA-MEM alignment output BAM file",
        sorted_bai: "Index files for the sorted BWA-MEM alignment BAM files"
    }
  }

  parameter_meta {
    bwa_genome_tar: "Compressed tarball containing BWA genome index"
    reference_fasta: "Reference genome FASTA file"
    reads: "FASTQ file for forward (R1) reads or interleaved reads"
    name: "Sample name for read group information"
    mates: "Optional FASTQ file for reverse (R2) reads"
    paired_end: "Optional boolean indicating if reads are paired end (default: true)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bwa_genome_tar
    File reference_fasta
    File reads
    String name
    File? mates
    Boolean paired_end = true
    Int cpu_cores = 8
    Int memory_gb = 16
  }

    # Name of reference FASTA file, which should be in bwa_genome_tar
  String ref_name = basename(reference_fasta)

   # Compute cpu_threads as one less than cpu_cores, with minimum of 1
  Int cpu_threads = if cpu_cores > 1 then cpu_cores - 1 else 1

  command <<<
    set -eo pipefail

    echo "Extracting BWA reference..."
    tar -xvf "~{bwa_genome_tar}"

    echo "Starting BWA alignment..."

  if [[ "~{mates}" == "" && "~{paired_end}" == "true" ]]; then
      # Interleaved (paired-end)
      bwa mem -p -v 3 -t ~{cpu_threads} -M -R "@RG\tID:~{name}\tSM:~{name}\tPL:illumina" \
        "bwa_index/~{ref_name}" "~{reads}" > "~{name}.sam"
    elif [[ "~{mates}" == "" && "~{paired_end}" == "false" ]]; then
      # Single-end
      bwa mem -v 3 -t ~{cpu_threads} -M -R "@RG\tID:~{name}\tSM:~{name}\tPL:illumina" \
        "bwa_index/~{ref_name}" "~{reads}" > "~{name}.sam"
    elif [[ "~{mates}" != "" && "~{paired_end}" == "true" ]]; then
      # Paired-end with forward and reverse fastqs
      bwa mem -v 3 -t ~{cpu_threads} -M -R "@RG\tID:~{name}\tSM:~{name}\tPL:illumina" \
        "bwa_index/~{ref_name}" "~{reads}" "~{mates}" > "~{name}.sam"
    else
      echo "Invalid input: Single-end experiments should only have one input FASTQ file."
      exit 1
    fi

    # Converting to BAM, sorting, and indexing
    samtools sort -@ ~{cpu_threads - 1} -o "~{name}.sorted_aligned.bam" "~{name}.sam"
    samtools index "~{name}.sorted_aligned.bam"
    
    # Cleaning up initial SAM file to save space
    rm -f "~{name}.sam"
  >>>

  output {
    File sorted_bam = "~{name}.sorted_aligned.bam"
    File sorted_bai = "~{name}.sorted_aligned.bam.bai"
  }

  runtime {
    docker: "getwilds/bwa:0.7.17"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
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
