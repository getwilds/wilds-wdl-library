## WILDS WDL for performing short-read alignment using Bowtie.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task bowtie_build {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Task for building Bowtie index files from a reference FASTA"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bowtie/ww-bowtie.wdl"
    outputs: {
        bowtie_index_tar: "Compressed tarball containing Bowtie genome index files"
    }
    topic: "genomics,transcriptomics,mapping"
    species: "any"
    operation: "indexing"
    in_sample_req: "none"
    in_sample_opt: "none"
    in_ref_req: "reference_fasta:nucleic_acid_sequence:fasta"
    in_ref_opt: "none"
    out_sample: "none"
    out_ref: "bowtie_index_tar:data_index:tar_format"
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file"
    index_prefix: "Prefix for the Bowtie index files"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File reference_fasta
    String index_prefix = "bowtie_index"
    Int cpu_cores = 4
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    mkdir -p bowtie_index

    echo "Building Bowtie index..."
    bowtie-build \
      --threads ~{cpu_cores} \
      "~{reference_fasta}" \
      "bowtie_index/~{index_prefix}"

    tar -czf bowtie_index.tar.gz bowtie_index/

    # Clean up index directory to save space
    rm -rf bowtie_index
  >>>

  output {
    File bowtie_index_tar = "bowtie_index.tar.gz"
  }

  runtime {
    docker: "getwilds/bowtie:1.3.1"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task bowtie_align {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Task for aligning short reads to a reference genome using Bowtie"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bowtie/ww-bowtie.wdl"
    outputs: {
        sorted_bam: "Sorted Bowtie alignment output BAM file",
        sorted_bai: "Index file for the sorted Bowtie alignment BAM file"
    }
    topic: "genomics,transcriptomics,mapping"
    species: "any"
    operation: "sequence_alignment"
    in_sample_req: "reads:nucleic_acid_sequence:fastq"
    in_sample_opt: "mates:nucleic_acid_sequence:fastq"
    in_ref_req: "bowtie_index_tar:data_index:tar_format"
    in_ref_opt: "none"
    out_sample: "sorted_bam:nucleic_acid_sequence_alignment:bam,sorted_bai:data_index:bai"
    out_ref: "none"
  }

  parameter_meta {
    bowtie_index_tar: "Compressed tarball containing Bowtie genome index files"
    index_prefix: "Prefix used when building the Bowtie index"
    reads: "FASTQ file for forward (R1) reads"
    name: "Sample name for output file naming and read group information"
    mates: "Optional FASTQ file for reverse (R2) reads for paired-end alignment"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bowtie_index_tar
    String index_prefix = "bowtie_index"
    File reads
    String name
    File? mates
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    echo "Extracting Bowtie index..."
    tar -xzf "~{bowtie_index_tar}"

    echo "Starting Bowtie alignment..."
    if [[ "~{mates}" == "" ]]; then
      # Single-end alignment
      bowtie \
        --threads ~{cpu_cores} \
        --sam \
        --best \
        -k 1 \
        -q \
        "bowtie_index/~{index_prefix}" \
        "~{reads}" \
        "~{name}.sam"
    else
      # Paired-end alignment
      bowtie \
        --threads ~{cpu_cores} \
        --sam \
        --best \
        -k 1 \
        -q \
        "bowtie_index/~{index_prefix}" \
        -1 "~{reads}" \
        -2 "~{mates}" \
        "~{name}.sam"
    fi

    # Convert to sorted BAM and index
    samtools sort -@ ~{cpu_cores} -o "~{name}.sorted.bam" "~{name}.sam"
    samtools index "~{name}.sorted.bam"

    # Clean up intermediate files
    rm -rf bowtie_index
    rm -f "~{name}.sam"
  >>>

  output {
    File sorted_bam = "~{name}.sorted.bam"
    File sorted_bai = "~{name}.sorted.bam.bai"
  }

  runtime {
    docker: "getwilds/bowtie:1.3.1"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
