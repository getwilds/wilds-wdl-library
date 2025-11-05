## WILDS WDL for RNA-SeQC quality control analysis of aligned RNA-seq data.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task run_rnaseqc {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Run RNA-SeQC quality control metrics on aligned RNA-seq data"
    outputs: {
        rnaseqc_metrics: "Compressed tarball containing RNA-SeQC quality metrics and coverage reports"
    }
  }

  parameter_meta {
    bam_file: "Aligned reads in BAM format"
    bam_index: "Index file for the aligned BAM file"
    ref_gtf: "Reference genome GTF annotation file (collapsed)"
    sample_name: "Sample name for output files"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bam_file
    File bam_index
    File ref_gtf
    String sample_name
    Int cpu_cores = 2
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    echo "Running RNA-SeQC..."
    rnaseqc "~{ref_gtf}" "~{bam_file}" OUTPUT \
      --sample="~{sample_name}" \
      --coverage

    tar -cvzf "~{sample_name}.QC.tar.gz" OUTPUT/*
  >>>

  output {
    File rnaseqc_metrics = "~{sample_name}.QC.tar.gz"
  }

  runtime {
    docker: "getwilds/rnaseqc:2.4.2"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
