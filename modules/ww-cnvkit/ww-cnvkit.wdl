## WILDS WDL CNVkit Module
## CNVkit module for copy number variation detection from targeted DNA sequencing data
## This module implements CNVkit workflows for germline and somatic CNV calling
## Uses the CNVkit software package for comprehensive CNV analysis

version 1.0

# Import testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### WORKFLOW DEFINITION ####

workflow cnvkit_example {
  meta {
    author: "WILDS Development Team"
    email: "wilds@fredhutch.org"
    description: "CNVkit copy number variation detection workflow for targeted sequencing data"
    url: "https://github.com/getwilds/wilds-wdl-library"
    outputs: {
        cnv_calls: "CNV call files (.cns files)",
        cnv_segments: "CNV segment files (.cnr files)",
        cnv_plot: "CNV visualization plots",
        reference_file: "CNVkit reference file (.cnn)"
    }
  }

  # Auto-download test data for testing purposes (multiple samples for reference)
  call ww_testdata.download_ref_data as download_ref { }
  call ww_testdata.download_bam_data as download_normal_1 {
    input: filename = "normal_sample_1.bam"
  }
  call ww_testdata.download_bam_data as download_normal_2 {
    input: filename = "normal_sample_2.bam"
  }
  call ww_testdata.download_bam_data as download_normal_3 {
    input: filename = "normal_sample_3.bam"
  }

  # Create reference using multiple normal samples
  call create_reference {
    input:
      bam_files = [download_normal_1.bam, download_normal_2.bam, download_normal_3.bam],
      bam_indices = [download_normal_1.bai, download_normal_2.bai, download_normal_3.bai],
      reference_fasta = download_ref.fasta,
      reference_fasta_index = download_ref.fasta_index,
      cpu_cores = 1,
      memory_gb = 4
  }

  # Process single demo sample
  call run_cnvkit {
    input:
      sample_name = "demo_sample",
      tumor_bam = download_normal_1.bam,
      tumor_bai = download_normal_1.bai,
      reference_cnn = create_reference.reference_cnn,
      cpu_cores = 1,
      memory_gb = 4
  }

  output {
    File cnv_calls = run_cnvkit.cnv_calls
    File cnv_segments = run_cnvkit.cnv_segments
    File cnv_plot = run_cnvkit.cnv_plot
    File reference_file = create_reference.reference_cnn
  }
}

#### TASK DEFINITIONS ####

task create_reference {
  meta {
    description: "Create CNVkit reference from normal samples or pooled reference"
    outputs: {
        reference_cnn: "CNVkit reference file (.cnn)"
    }
  }

  parameter_meta {
    bam_files: "Array of BAM files for reference creation"
    bam_indices: "Array of BAM index files corresponding to BAM files"
    target_bed: "Target regions BED file"
    antitarget_bed: "Antitarget regions BED file"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Reference genome FASTA index file"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    Array[File] bam_files
    Array[File] bam_indices
    File? target_bed
    File? antitarget_bed
    File? reference_fasta
    File? reference_fasta_index
    Int cpu_cores = 4
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    # Link FASTA index alongside FASTA file if both are provided
    if [[ -n "~{reference_fasta}" && -n "~{reference_fasta_index}" ]]; then
      # Create symbolic links with consistent naming
      ln -sf ~{reference_fasta} reference.fa
      ln -sf ~{reference_fasta_index} reference.fa.fai
      fasta_arg="--fasta reference.fa"
    else
      fasta_arg=""
    fi

    # Link BAM files and their indices to working directory
    bam_array=(~{sep=' ' bam_files})
    bai_array=(~{sep=' ' bam_indices})
    linked_bams=""

    for i in "${!bam_array[@]}"; do
      bam_file="${bam_array[$i]}"
      bai_file="${bai_array[$i]}"
      bam_basename=$(basename "$bam_file")

      # Create symbolic links for BAM and BAI files
      ln -sf "$bam_file" "$bam_basename"
      ln -sf "$bai_file" "${bam_basename}.bai"

      linked_bams="$linked_bams $bam_basename"
    done

    # Create reference using CNVkit with linked BAM files
    # Use autodetect mode if no target BED is provided
    if [[ -n "~{target_bed}" ]]; then
      # Use provided target regions
      cnvkit.py batch \
        $linked_bams \
        --normal \
        --targets ~{target_bed} \
        ~{"--antitargets " + antitarget_bed} \
        ${fasta_arg} \
        --output-reference reference.cnn \
        --output-dir . \
        --processes ~{cpu_cores}
    else
      # Auto-detect targets from BAM files (WGS/WES mode)
      cnvkit.py batch \
        $linked_bams \
        --normal \
        --method wgs \
        ${fasta_arg} \
        --output-reference reference.cnn \
        --output-dir . \
        --processes ~{cpu_cores}
    fi
  >>>

  output {
    File reference_cnn = "reference.cnn"
  }

  runtime {
    docker: "getwilds/cnvkit:0.9.10"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task run_cnvkit {
  meta {
    description: "Run CNVkit copy number analysis on tumor sample"
    outputs: {
        cnv_calls: "CNV calls file (.cns)",
        cnv_segments: "CNV segments file (.cnr)",
        cnv_plot: "CNV visualization plot"
    }
  }

  parameter_meta {
    sample_name: "Sample identifier"
    tumor_bam: "Tumor BAM file"
    tumor_bai: "Tumor BAM index file"
    normal_bam: "Optional normal/control BAM file"
    normal_bai: "Optional normal/control BAM index file"
    reference_cnn: "CNVkit reference file"
    target_bed: "Optional target regions BED file"
    paired_analysis: "Whether to perform paired tumor/normal analysis"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    String sample_name
    File tumor_bam
    File tumor_bai
    File? normal_bam
    File? normal_bai
    File reference_cnn
    File? target_bed
    Boolean paired_analysis = false
    Int cpu_cores = 4
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    # Run CNVkit analysis
    if [[ "~{paired_analysis}" == "true" && -n "~{normal_bam}" ]]; then
      # Paired tumor/normal analysis
      cnvkit.py batch \
        ~{tumor_bam} \
        --normal ~{normal_bam} \
        --reference ~{reference_cnn} \
        ~{"--targets " + target_bed} \
        --output-dir . \
        --processes ~{cpu_cores}
    else
      # Tumor-only analysis
      cnvkit.py batch \
        ~{tumor_bam} \
        --reference ~{reference_cnn} \
        ~{"--targets " + target_bed} \
        --output-dir . \
        --processes ~{cpu_cores}
    fi

    # Generate segments and calls
    tumor_base=$(basename ~{tumor_bam} .bam)

    # Call copy number segments
    cnvkit.py call "${tumor_base}.cnr" \
      --output "${tumor_base}.cns" \
      --method threshold

    # Generate visualization plot
    cnvkit.py scatter "${tumor_base}.cnr" \
      --segment "${tumor_base}.cns" \
      --output "~{sample_name}_cnv_plot.pdf"

    # Rename outputs with sample name
    mv "${tumor_base}.cnr" "~{sample_name}.cnr"
    mv "${tumor_base}.cns" "~{sample_name}.cns"
  >>>

  output {
    File cnv_segments = "~{sample_name}.cnr"
    File cnv_calls = "~{sample_name}.cns"
    File cnv_plot = "~{sample_name}_cnv_plot.pdf"
  }

  runtime {
    docker: "getwilds/cnvkit:0.9.10"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}