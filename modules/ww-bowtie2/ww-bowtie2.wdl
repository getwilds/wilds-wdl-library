## WILDS WDL for performing sequence alignment using Bowtie 2.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task bowtie2_build {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Task for building Bowtie 2 index files from a reference FASTA"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bowtie2/ww-bowtie2.wdl"
    outputs: {
        bowtie2_index_tar: "Compressed tarball containing Bowtie 2 genome index files"
    }
    topic: "genomics,transcriptomics,mapping"
    species: "any"
    operation: "indexing"
    in_sample_data: "none"
    in_sample_format: "none"
    in_ref_data: "dna_sequence"
    in_ref_format: "fasta"
    out_sample_data: "none"
    out_sample_format: "none"
    out_ref_data: "data_index"
    out_ref_format: "tar_format"
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file"
    index_prefix: "Prefix for the Bowtie 2 index files"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File reference_fasta
    String index_prefix = "bowtie2_index"
    Int cpu_cores = 4
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    mkdir -p bowtie2_index

    echo "Building Bowtie 2 index..."
    bowtie2-build \
      --threads ~{cpu_cores} \
      "~{reference_fasta}" \
      "bowtie2_index/~{index_prefix}"

    tar -czf bowtie2_index.tar.gz bowtie2_index/

    # Clean up index directory to save space
    rm -rf bowtie2_index
  >>>

  output {
    File bowtie2_index_tar = "bowtie2_index.tar.gz"
  }

  runtime {
    docker: "getwilds/bowtie2:2.5.4"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task bowtie2_align {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Task for aligning sequence reads to a reference genome using Bowtie 2"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bowtie2/ww-bowtie2.wdl"
    outputs: {
        sorted_bam: "Sorted Bowtie 2 alignment output BAM file",
        sorted_bai: "Index file for the sorted Bowtie 2 alignment BAM file"
    }
    topic: "genomics,transcriptomics,mapping"
    species: "any"
    operation: "sequence_alignment"
    in_sample_data: "nucleic_acid_sequence"
    in_sample_format: "fastq"
    in_ref_data: "data_index"
    in_ref_format: "tar_format"
    out_sample_data: "nucleic_acid_sequence_alignment,data_index"
    out_sample_format: "bam,bai"
    out_ref_data: "none"
    out_ref_format: "none"
  }

  parameter_meta {
    bowtie2_index_tar: "Compressed tarball containing Bowtie 2 genome index files"
    index_prefix: "Prefix used when building the Bowtie 2 index"
    reads: "FASTQ file for forward (R1) reads"
    name: "Sample name for output file naming and read group information"
    mates: "Optional FASTQ file for reverse (R2) reads for paired-end alignment"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bowtie2_index_tar
    String index_prefix = "bowtie2_index"
    File reads
    String name
    File? mates
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    echo "Extracting Bowtie 2 index..."
    tar -xzf "~{bowtie2_index_tar}"

    echo "Starting Bowtie 2 alignment..."
    if [[ "~{mates}" == "" ]]; then
      # Single-end alignment
      bowtie2 \
        --threads ~{cpu_cores} \
        --rg-id "~{name}" \
        --rg "SM:~{name}" \
        --rg "PL:illumina" \
        -x "bowtie2_index/~{index_prefix}" \
        -U "~{reads}" \
        -S "~{name}.sam"
    else
      # Paired-end alignment
      bowtie2 \
        --threads ~{cpu_cores} \
        --rg-id "~{name}" \
        --rg "SM:~{name}" \
        --rg "PL:illumina" \
        -x "bowtie2_index/~{index_prefix}" \
        -1 "~{reads}" \
        -2 "~{mates}" \
        -S "~{name}.sam"
    fi

    # Convert to sorted BAM and index
    samtools sort -@ ~{cpu_cores} -o "~{name}.sorted.bam" "~{name}.sam"
    samtools index "~{name}.sorted.bam"

    # Clean up intermediate files
    rm -rf bowtie2_index
    rm -f "~{name}.sam"
  >>>

  output {
    File sorted_bam = "~{name}.sorted.bam"
    File sorted_bai = "~{name}.sorted.bam.bai"
  }

  runtime {
    docker: "getwilds/bowtie2:2.5.4"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
