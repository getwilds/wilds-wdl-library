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
    description: "Aligns reads to a reference using Bowtie 2, optionally filtering the output using Samtools"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bowtie2/ww-bowtie2.wdl"
    outputs: {
        sorted_bam: "Sorted Bowtie 2 alignment output BAM file",
        sorted_bai: "Index file for the sorted Bowtie 2 alignment BAM file",
        unaligned_se: "FASTQ of unaligned reads (single-end input). Only produced when capture_unaligned is true and mates is unset.",
        unaligned_r1: "FASTQ of R1 reads that failed to align concordantly. Only produced when capture_unaligned is true and mates is set.",
        unaligned_r2: "FASTQ of R2 reads that failed to align concordantly. Only produced when capture_unaligned is true and mates is set."
    }
  }

  parameter_meta {
    bowtie2_index_tar: "Compressed tarball containing Bowtie 2 genome index files"
    index_prefix: "Prefix used when building the Bowtie 2 index"
    reads: "FASTQ file for forward (R1) reads"
    name: "Sample name for output file naming and read group information"
    mates: "Optional FASTQ file for reverse (R2) reads for paired-end alignment"
    preset: "Optional bowtie2 sensitivity preset (fast, sensitive, very-sensitive, fast-local, sensitive-local, very-sensitive-local). Defaults to bowtie2's end-to-end --sensitive when unset."
    capture_unaligned: "If true, write reads that fail to align concordantly to gzipped FASTQ outputs (--un-gz / --un-conc-gz). Used for e.g. rRNA depletion upstream of a second alignment step."
    min_mapq: "Minimum MAPQ score for samtools view -q post-alignment filter. 0 means no filter (default)."
    samtools_filter_flags: "Extra flags passed to samtools view for post-alignment filtering (e.g. '-f 2' to keep only proper pairs). Empty by default."
    extra_bowtie2_args: "Additional arguments forwarded verbatim to bowtie2 (e.g. '--no-mixed --no-discordant'). Empty by default."
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bowtie2_index_tar
    String index_prefix = "bowtie2_index"
    File reads
    String name
    File? mates
    String? preset
    Boolean capture_unaligned = false
    Int min_mapq = 0
    String samtools_filter_flags = ""
    String extra_bowtie2_args = ""
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    echo "Extracting Bowtie 2 index..."
    tar -xzf "~{bowtie2_index_tar}"

    PRESET_FLAG=""
    if [ ! -z "~{preset}" ]; then
      PRESET_FLAG="--~{preset}"
    fi

    UNALIGNED_FLAG=""
    if [ "~{capture_unaligned}" == "true" ]; then
      if [[ "~{mates}" == "" ]]; then
        UNALIGNED_FLAG="--un-gz ~{name}.unaligned.fq.gz"
      else
        UNALIGNED_FLAG="--un-conc-gz ~{name}.unaligned.R%.fq.gz"
      fi
    fi

    echo "Starting Bowtie 2 alignment..."
    if [[ "~{mates}" == "" ]]; then
      # Single-end alignment
      bowtie2 \
        --threads ~{cpu_cores} \
        --rg-id "~{name}" \
        --rg "SM:~{name}" \
        --rg "PL:illumina" \
        ${PRESET_FLAG} \
        ${UNALIGNED_FLAG} \
        ~{extra_bowtie2_args} \
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
        ${PRESET_FLAG} \
        ${UNALIGNED_FLAG} \
        ~{extra_bowtie2_args} \
        -x "bowtie2_index/~{index_prefix}" \
        -1 "~{reads}" \
        -2 "~{mates}" \
        -S "~{name}.sam"
    fi

    # Apply optional MAPQ + SAM-flag filtering, then sort + index
    samtools view -h -b ~{samtools_filter_flags} -q ~{min_mapq} "~{name}.sam" \
      | samtools sort -@ ~{cpu_cores} -o "~{name}.sorted.bam" -
    samtools index "~{name}.sorted.bam"

    # Clean up intermediate files
    rm -rf bowtie2_index
    rm -f "~{name}.sam"
  >>>

  output {
    File sorted_bam = "~{name}.sorted.bam"
    File sorted_bai = "~{name}.sorted.bam.bai"
    # Bowtie 2 emits one of three filenames depending on inputs; the unused two
    # never exist on disk so the engine resolves them to optional None.
    File? unaligned_se = "~{name}.unaligned.fq.gz"
    File? unaligned_r1 = "~{name}.unaligned.R1.fq.gz"
    File? unaligned_r2 = "~{name}.unaligned.R2.fq.gz"
  }

  runtime {
    docker: "getwilds/bowtie2:2.5.4"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
