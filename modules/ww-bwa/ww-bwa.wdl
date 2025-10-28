## WILDS WDL for performing sequence alignment using BWA-MEM.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task bwa_index {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
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

    # Cleaning up the bwa_index directory to save space, only need the tarball
    rm -rf bwa_index
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
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
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

    # Cleaning up BWA index and initial SAM file to save space
    rm -rf bwa_index  # Remove BWA index directory
    rm -f "~{name}.sam"  # Remove initial SAM file
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
