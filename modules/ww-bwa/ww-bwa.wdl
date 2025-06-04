## WILDS WDL for performing sequence alignment using BWA-MEM.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

struct SampleInfo {
    String name
    File r1
    File r2
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
    reference_fasta: "Reference genome FASTA file"
    samples: "List of sample objects, each containing sample name and R1/R2 FASTQ files"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    memory_gb: "Memory allocated for each task in the workflow in GB"
  }

  input {
    File reference_fasta
    Array[SampleInfo] samples
    Int cpus = 8
    Int memory_gb = 32
  }

  call bwa_index { input:
      reference_fasta = reference_fasta,
      cpu_cores = cpus,
      memory_gb = memory_gb
  }

  scatter (sample in samples) {
    call bwa_mem { input:
        bwa_genome_tar = bwa_index.bwa_index_tar,
        reference_fasta = reference_fasta,
        sample_data = sample,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

    call validate_outputs { input:
    sample_names = bwa_mem.name,
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
        name: "Sample name that was processed",
        sorted_bam: "Sorted BWA-MEM alignment output BAM file",
        sorted_bai: "Index files for the sorted BWA-MEM alignment BAM files"
    }
  }

  parameter_meta {
    bwa_genome_tar: "Compressed tarball containing BWA genome index"
    reference_fasta: "Reference genome FASTA file"
    sample_data: "Sample object containing sample name and R1/R2 FASTQ files"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bwa_genome_tar
    File reference_fasta
    SampleInfo sample_data
    Int cpu_cores = 8
    Int memory_gb = 16
  }

   # Name of reference FASTA file, which should be in bwa_genome_tar
  String ref_name = basename(reference_fasta)

   # Compute cpu_threads as one less than cpu_cores, with minimum of 1
  Int cpu_threads = if cpu_cores > 1 then cpu_cores - 1 else 1

  command <<<
  set -eo pipefail && \
  echo "Extracting BWA reference..." && \
  tar -xvf "~{bwa_genome_tar}" && \
  echo "Starting BWA alignment..." && \
  bwa mem -v 3 -t ~{cpu_threads} -M -R '@RG\tID:~{sample_data.name}\tSM:~{sample_data.name}\tPL:illumina' "bwa_index/~{ref_name}" "~{sample_data.r1}" "~{sample_data.r2}" > "~{sample_data.name}.sam" && \
  samtools sort -@ ~{cpu_threads-1} -o "~{sample_data.name}".sorted_aligned.bam "~{sample_data.name}".sam && \
  samtools index "~{sample_data.name}".sorted_aligned.bam
  >>>

  output {
    String name = sample_data.name
    File sorted_bam = "~{sample_data.name}.sorted_aligned.bam"
    File sorted_bai = "~{sample_data.name}.sorted_aligned.bam.bai"
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
    sample_names: "Array of sample names that were processed"
  }

  input {
    Array[File] bam_files
    Array[File] bai_files
    Array[String] sample_names
  }

  command <<<
    set -eo pipefail

    echo "=== BWA-MEM Alignment Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    # Arrays for bash processing
    sample_names=("~{sep=" " sample_names}")
    bam_files=("~{sep=" " bam_files}")
    bai_files=("~{sep=" " bai_files}")

    validation_passed=true
    total_mapped_reads=0

    # Check each sample
    for i in "${!sample_names[@]}"; do
      sample_name="${sample_names[$i]}"
      bam_file="${bam_files[$i]}"
      bai_file="${bai_files[$i]}"

      echo "--- Sample: $sample_name ---" >> validation_report.txt

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
    echo "Total samples processed: ${#sample_names[@]}" >> validation_report.txt
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
