## WILDS WDL CNVkit Module
## CNVkit module for copy number variation detection from targeted DNA sequencing data
## This module implements CNVkit workflows for germline and somatic CNV calling
## Uses the CNVkit software package for comprehensive CNV analysis

version 1.0

task create_reference {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    File reference_fasta
    File reference_fasta_index
    File? target_bed
    File? antitarget_bed
    Int cpu_cores = 4
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    # Link FASTA and index files to working directory
    # If soft links aren't allowed on your HPC system, copy them locally instead
    ln -s "~{reference_fasta}" "~{basename(reference_fasta)}"
    ln -s "~{reference_fasta_index}" "~{basename(reference_fasta_index)}"

    # Link BAM files and their indices to working directory
    # If soft links aren't allowed on your HPC system, copy them locally instead
    bam_array=(~{sep=' ' bam_files})
    bai_array=(~{sep=' ' bam_indices})
    linked_bams=""

    for i in "${!bam_array[@]}"; do
      bam_file="${bam_array[$i]}"
      bai_file="${bai_array[$i]}"
      bam_basename=$(basename "$bam_file")

      # Create symbolic links for BAM and BAI files
      ln -s "$bam_file" "$bam_basename"
      ln -s "$bai_file" "${bam_basename}.bai"

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
        --fasta ~{basename(reference_fasta)} \
        --output-reference reference.cnn \
        --output-dir . \
        --processes ~{cpu_cores}
    else
      # Auto-detect targets from BAM files (WGS/WES mode)
      cnvkit.py batch \
        $linked_bams \
        --normal \
        --method wgs \
        --fasta ~{basename(reference_fasta)} \
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
