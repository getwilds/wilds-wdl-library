## WILDS WDL for processing genomic files with Samtools.
## Intended for use alone or as a modular component in the WILDS ecosystem.

version 1.0

task crams_to_fastq {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Merge CRAM/BAM/SAM files and convert to FASTQ's using samtools."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-samtools/ww-samtools.wdl"
    outputs: {
        r1_fastq: "R1 FASTQ file generated from merged CRAM/BAM/SAM file",
        r2_fastq: "R2 FASTQ file generated from merged CRAM/BAM/SAM file",
        sample_name: "Sample name that was processed"
    }
    topic: "genomics,transcriptomics"
    species: "any"
    operation: "data_formatting"
    input_sample_required: "cram_files:nucleic_acid_sequence_alignment:cram|bam|sam"
    input_sample_optional: "none"
    input_reference_required: "ref:nucleic_acid_sequence:fasta"
    input_reference_optional: "none"
    output_sample: "r1_fastq:nucleic_acid_sequence:fastq,r2_fastq:nucleic_acid_sequence:fastq"
    output_reference: "none"
  }

  parameter_meta {
    cram_files: "List of CRAM/BAM/SAM files for a sample"
    ref: "Reference genome FASTA file"
    name: "Name of the sample (used for file naming)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File ref
    Array[File] cram_files
    String name
    Int cpu_cores = 2
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    # Symlink reference FASTA locally so samtools can write the .fai index
    ln -s "~{ref}" "~{basename(ref)}"

    # Symlink CRAM/BAM/SAM files locally to avoid read-only filesystem issues
    cram_files=(~{sep=" " cram_files})
    local_crams=""
    for cram in "${cram_files[@]}"; do
      ln -s "${cram}" "$(basename "${cram}")"
      local_crams="${local_crams} $(basename "${cram}")"
    done

    # Merge CRAM/BAM/SAM files if more than one, then collate and convert to FASTQ
    samtools merge -@ ~{cpu_cores} --reference "~{basename(ref)}" -u - ${local_crams} | \
    samtools collate -@ ~{cpu_cores} --reference "~{basename(ref)}" -O - | \
    samtools fastq -@ ~{cpu_cores} --reference "~{basename(ref)}" -1 "~{name}_R1.fastq.gz" -2 "~{name}_R2.fastq.gz" -0 /dev/null -s /dev/null -
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

task merge_bams_to_cram {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Merge multiple BAM files into a single CRAM file using samtools merge"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-samtools/ww-samtools.wdl"
    outputs: {
        cram: "Merged CRAM file containing all reads from input BAMs",
        crai: "Index file for the merged CRAM"
    }
    topic: "genomics,transcriptomics"
    species: "any"
    operation: "aggregation"
    input_sample_required: "bams_to_merge:nucleic_acid_sequence_alignment:bam"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "cram:nucleic_acid_sequence_alignment:cram,crai:data_index:crai"
    output_reference: "none"
  }

  parameter_meta {
    bams_to_merge: "Array of BAM files to merge into a single CRAM file"
    base_file_name: "Base name for output CRAM file"
    cpu_cores: "Number of CPU cores to use (threads = cpu_cores - 1)"
    memory_gb: "Memory allocation in GB"
  }

  input {
    Array[File] bams_to_merge
    String base_file_name
    Int cpu_cores = 6
    Int memory_gb = 12
  }

  command <<<
    set -eo pipefail

    samtools merge -@ ~{cpu_cores - 1} \
      --write-index --output-fmt CRAM \
      "~{base_file_name}.merged.cram" ~{sep=" " bams_to_merge}
  >>>

  output {
    File cram = "~{base_file_name}.merged.cram"
    File crai = "~{base_file_name}.merged.cram.crai"
  }

  runtime {
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
    docker: "getwilds/samtools:1.19"
  }
}

task mpileup {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Generate samtools mpileup from a BAM file"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-samtools/ww-samtools.wdl"
    outputs: {
        pileup: "Pileup file"
    }
    topic: "genomics,transcriptomics"
    species: "any"
    operation: "quantification"
    input_sample_required: "bamfile:nucleic_acid_sequence_alignment:bam|cram"
    input_sample_optional: "none"
    input_reference_required: "ref_fasta:nucleic_acid_sequence:fasta"
    input_reference_optional: "none"
    output_sample: "pileup:nucleic_acid_sequence_alignment:pileup"
    output_reference: "none"
  }

  parameter_meta {
    bamfile: "Input BAM or CRAM file"
    ref_fasta: "Reference genome FASTA file"
    sample_name: "Name of the sample (used for output file naming)"
    disable_baq: "Whether to disable per-Base Alignment Quality"
    min_mapq: "Minimum mapping quality for alignments to be included"
    min_baseq: "Minimum base quality for bases to be included"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File bamfile
    File ref_fasta
    String sample_name
    Boolean disable_baq = false
    Int min_mapq = 0
    Int min_baseq = 13
    Int cpu_cores = 2
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Symlink reference FASTA locally so samtools can write the .fai index
    ln -s "~{ref_fasta}" "~{basename(ref_fasta)}"

    samtools mpileup \
      -f "~{basename(ref_fasta)}" \
      -q "~{min_mapq}" \
      -Q "~{min_baseq}" \
      ~{if disable_baq then "--no-BAQ" else ""} \
      "~{bamfile}" \
      > "~{sample_name}.pileup"
  >>>

  output {
    File pileup = "~{sample_name}.pileup"
  }

  runtime {
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
    docker: "getwilds/samtools:1.19"
  }
}
