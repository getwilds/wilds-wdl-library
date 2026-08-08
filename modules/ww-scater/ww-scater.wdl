## WILDS WDL Module: ww-scater
## Description: Module for single-cell RNA-seq dimensionality reduction (PCA, UMAP)
## and QC/reduced-dimension visualization using scater

version 1.0

task run_scater {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Perform PCA and UMAP dimensionality reduction on a normalized SingleCellExperiment object and generate QC and reduced-dimension diagnostic plots with scater"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-scater/ww-scater.wdl"
    outputs: {
        sce_object: "SingleCellExperiment RDS object with PCA and UMAP embeddings and per-cell QC metrics",
        pca_plot: "Scatter plot of the first two PCA components",
        umap_plot: "Scatter plot of the UMAP embedding",
        qc_plot: "Per-cell QC scatter plot (library size vs. detected genes, colored by total counts)"
    }
    topic: "transcriptomics,single_cell"
    species: "human,eukaryote"
    operation: "dimensionality_reduction"
    input_sample_required: "sce_rds:gene_expression_matrix:rds"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "sce_object:gene_expression_matrix:rds_format,pca_plot:plot:png,umap_plot:plot:png,qc_plot:quality_control_report:png"
    output_reference: "none"
  }

  parameter_meta {
    sce_rds: "Normalized SingleCellExperiment RDS object (e.g. from ww-scran) containing a logcounts assay"
    sample_name: "Sample name used for output file prefixes"
    n_pcs: "Number of principal components to compute"
    n_hvgs: "Number of most variable genes (by log-expression variance) to use for PCA"
    random_seed: "Random seed for UMAP, for reproducibility"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
    docker_image: "Docker image to use for this task"
  }

  input {
    File sce_rds
    String sample_name
    Int n_pcs = 50
    Int n_hvgs = 2000
    Int random_seed = 100
    Int memory_gb = 8
    Int cpu_cores = 2
    String docker_image = "getwilds/scater:1.40.2"
  }

  command <<<
    set -eo pipefail

    curl -so scater_analysis.R \
      "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-scater/modules/ww-scater/scater_analysis.R"

    Rscript scater_analysis.R \
      --sce_rds="~{sce_rds}" \
      --sample_name="~{sample_name}" \
      --n_pcs=~{n_pcs} \
      --n_hvgs=~{n_hvgs} \
      --random_seed=~{random_seed} \
      --output_prefix="~{sample_name}_scater"
  >>>

  output {
    File sce_object = "~{sample_name}_scater.rds"
    File pca_plot = "~{sample_name}_scater_pca.png"
    File umap_plot = "~{sample_name}_scater_umap.png"
    File qc_plot = "~{sample_name}_scater_qc.png"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
