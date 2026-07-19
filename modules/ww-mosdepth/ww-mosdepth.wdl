## WILDS WDL Mosdepth Module
## Fast BAM/CRAM depth calculation.
## Based on mosdepth by Brent Pedersen.

version 1.0

#### TASK DEFINITIONS ####

task calculate_depth {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Calculates sequencing coverage depth from BAM or CRAM alignments."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-mosdepth/ww-mosdepth.wdl"
    outputs: {
        depth_per_base: "Per-base coverage depth in gzipped BED format",
        depth_summary: "Summary of coverage depth statistics (plain text)",
        region_depth: "Coverage depth per window or region (produced when --by is passed)"
    }
    topic: "any"
    species: "human,eukaryote,prokaryote,virus"
    operation: "data_parsing"
    input_sample_required: "input_bam:bam:bam"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "ref_fasta:fasta:fasta"
    output_sample: "depth_per_base:text_data:bed"
    output_reference: "none"
  }

  parameter_meta {
    sample_name: "Name identifier for the sample"
    input_bam: "Input BAM or CRAM file"
    input_bam_index: "Index file for the input BAM/CRAM"
    ref_fasta: "Reference genome FASTA file (required for CRAM inputs)"
    regions_bed: "Optional BED file specifying regions of interest"
    window_size: "Window size for windowed depth output (default: 100bp)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    String sample_name
    File input_bam
    File input_bam_index
    File? ref_fasta
    File? regions_bed
    Int window_size = 100
    Int cpu_cores = 2
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    mosdepth \
      ~{if defined(ref_fasta) then "--fasta " + ref_fasta else ""} \
      ~{if defined(regions_bed) then "--by " + regions_bed else "--by " + window_size} \
      --threads ~{cpu_cores} \
      ~{sample_name} \
      ~{input_bam}
  >>>

  output {
    File depth_per_base = "~{sample_name}.per-base.bed.gz"
    File depth_summary = "~{sample_name}.mosdepth.summary.txt"
    File? region_depth = "~{sample_name}.regions.bed.gz"
  }

  runtime {
    docker: "getwilds/mosdepth:0.3.14"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
