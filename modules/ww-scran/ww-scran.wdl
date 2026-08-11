## WILDS WDL Module: ww-scran
## Description: Module for single-cell RNA-seq QC, normalization, and feature
## selection using scran

version 1.0

task run_scran {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Perform per-cell QC filtering, deconvolution-based normalization, and highly variable gene selection on single-cell RNA-seq data with scran"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-scran/ww-scran.wdl"
    outputs: {
        sce_object: "Normalized SingleCellExperiment RDS object with size factors and PCA",
        qc_plot: "Per-cell QC scatter plot (library size vs. detected genes, colored by discard status)",
        size_factor_plot: "Histogram of scran deconvolution size factors",
        mean_variance_plot: "Mean-variance trend plot from highly variable gene modeling",
        hvg_table: "CSV of all genes ranked by biological variance component"
    }
    topic: "transcriptomics,single_cell"
    species: "human,eukaryote"
    operation: "expression_analysis"
    input_sample_required: "sce_rds:gene_expression_matrix:binary_format"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "sce_object:gene_expression_matrix:binary_format,qc_plot:quality_control_report:png,size_factor_plot:plot:png,mean_variance_plot:plot:png,hvg_table:gene_report:csv"
    output_reference: "none"
  }

  parameter_meta {
    sce_rds: "SingleCellExperiment RDS object (e.g. from ww-dropletutils) containing the raw counts to QC and normalize"
    sample_name: "Sample name used for output file prefixes"
    mito_pattern: "Regex pattern identifying mitochondrial gene symbols"
    nmads: "Number of median absolute deviations from the median used to flag low-quality cells"
    min_mean: "Minimum mean expression for a gene to be used in size factor estimation and HVG modeling"
    n_hvgs: "Number of top highly variable genes to select for PCA"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
    docker_image: "Docker image to use for this task"
  }

  input {
    File sce_rds
    String sample_name
    String mito_pattern = "^MT-"
    Float nmads = 3.0
    Float min_mean = 0.1
    Int n_hvgs = 2000
    Int memory_gb = 8
    Int cpu_cores = 2
    String docker_image = "getwilds/scran:1.40.0"
  }

  command <<<
    set -eo pipefail

    curl -so scran_analysis.R \
      "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-scran/modules/ww-scran/scran_analysis.R"

    Rscript scran_analysis.R \
      --sce_rds="~{sce_rds}" \
      --sample_name="~{sample_name}" \
      --mito_pattern="~{mito_pattern}" \
      --nmads=~{nmads} \
      --min_mean=~{min_mean} \
      --n_hvgs=~{n_hvgs} \
      --output_prefix="~{sample_name}_scran"
  >>>

  output {
    File sce_object = "~{sample_name}_scran.rds"
    File qc_plot = "~{sample_name}_scran_qc.png"
    File size_factor_plot = "~{sample_name}_scran_size_factors.png"
    File mean_variance_plot = "~{sample_name}_scran_mean_variance.png"
    File hvg_table = "~{sample_name}_scran_hvgs.csv"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
