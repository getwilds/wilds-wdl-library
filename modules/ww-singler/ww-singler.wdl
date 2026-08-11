## WILDS WDL Module: ww-singler
## Description: Module for single-cell RNA-seq clustering, marker gene identification,
## and cell-type annotation using scran and SingleR

version 1.0

task run_singler {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Cluster cells, identify per-cluster marker genes, and annotate cell types on a dimensionality-reduced SingleCellExperiment object using scran and SingleR"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-singler/ww-singler.wdl"
    outputs: {
        sce_object: "SingleCellExperiment RDS object with cluster assignments and predicted cell types",
        cluster_table: "CSV mapping each cell barcode to its assigned cluster",
        marker_table: "CSV of top marker genes per cluster",
        prediction_table: "CSV of SingleR cell-type predictions and scores per cluster"
    }
    topic: "transcriptomics,single_cell"
    species: "human,eukaryote"
    operation: "annotation"
    input_sample_required: "sce_rds:gene_expression_matrix:rds"
    input_sample_optional: "reference_rds:gene_expression_matrix:rds"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "sce_object:gene_expression_matrix:rds_format,cluster_table:gene_report:csv,marker_table:gene_report:csv,prediction_table:gene_report:csv"
    output_reference: "none"
  }

  parameter_meta {
    sce_rds: "Dimensionality-reduced SingleCellExperiment RDS object (e.g. from ww-scater) containing a PCA reducedDim"
    sample_name: "Sample name used for output file prefixes"
    reference_rds: "Optional custom reference SingleCellExperiment RDS object with labeled cell types. If omitted, a celldex reference dataset is fetched by name"
    reference_dataset: "Name of the celldex reference dataset function to fetch (e.g. HumanPrimaryCellAtlasData, MouseRNAseqData), used when reference_rds is not provided"
    label_column: "Column name in the reference object's colData containing cell-type labels"
    n_top_markers: "Number of top marker genes to report per cluster"
    random_seed: "Random seed for clustering, for reproducibility"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
    docker_image: "Docker image to use for this task"
  }

  input {
    File sce_rds
    String sample_name
    File? reference_rds
    String reference_dataset = "HumanPrimaryCellAtlasData"
    String label_column = "label.main"
    Int n_top_markers = 10
    Int random_seed = 100
    Int memory_gb = 16
    Int cpu_cores = 2
    String docker_image = "getwilds/singler:2.14.1"
  }

  command <<<
    set -eo pipefail

    curl -so singler_analysis.R \
      "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-singler/modules/ww-singler/singler_analysis.R"

    Rscript singler_analysis.R \
      --sce_rds="~{sce_rds}" \
      --sample_name="~{sample_name}" \
      --reference_rds="~{default="" reference_rds}" \
      --reference_dataset="~{reference_dataset}" \
      --label_column="~{label_column}" \
      --n_top_markers=~{n_top_markers} \
      --random_seed=~{random_seed} \
      --output_prefix="~{sample_name}_singler"
  >>>

  output {
    File sce_object = "~{sample_name}_singler.rds"
    File cluster_table = "~{sample_name}_singler_clusters.csv"
    File marker_table = "~{sample_name}_singler_markers.csv"
    File prediction_table = "~{sample_name}_singler_predictions.csv"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
