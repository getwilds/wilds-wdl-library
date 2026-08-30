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
    input_sample_required: "sce_rds:gene_expression_matrix:binary_format"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "reference_rds:gene_expression_matrix:binary_format"
    output_sample: "sce_object:gene_expression_matrix:binary_format,cluster_table:gene_report:csv,marker_table:gene_report:csv,prediction_table:gene_report:csv"
    output_reference: "none"
  }

  parameter_meta {
    sce_rds: "Dimensionality-reduced SingleCellExperiment RDS object (e.g. from ww-scater) containing a PCA reducedDim"
    sample_name: "Sample name used for output file prefixes"
    reference_rds: "Optional custom reference SingleCellExperiment RDS object with labeled cell types. If omitted, a celldex reference dataset is fetched by name"
    reference_dataset: "Name of the celldex reference dataset function to fetch (e.g. HumanPrimaryCellAtlasData, MouseRNAseqData), used when reference_rds is not provided"
    reference_ensembl: "Fetch the celldex reference with Ensembl gene IDs instead of gene symbols, to match 10x Cell Ranger output"
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
    Boolean reference_ensembl = true
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
      "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/singler-debugging/modules/ww-singler/singler_analysis.R"

    # celldex/AnnotationHub/gypsum all take a write lock in their cache dir even
    # for a cache hit, and the getwilds/singler image's baked-in caches are
    # read-only under Apptainer. Seed a writable copy of each in the task dir
    # and repoint the corresponding env var (skipped for images without them).
    if [ -d /opt/hubcache/gypsum ]; then
      mkdir -p "${PWD}/gypsum-cache"
      cp -r /opt/hubcache/gypsum/. "${PWD}/gypsum-cache/"
      export GYPSUM_CACHE_DIR="${PWD}/gypsum-cache"
    fi
    if [ -d /opt/hubcache/annotationhub ]; then
      mkdir -p "${PWD}/annotationhub-cache"
      cp -r /opt/hubcache/annotationhub/. "${PWD}/annotationhub-cache/"
      export ANNOTATION_HUB_CACHE="${PWD}/annotationhub-cache"
      # Work purely from the local cache; skip the online metadata refresh and
      # its own read-only lock on the metadata SQLite DB.
      export ANNOTATION_HUB_LOCAL=TRUE
    fi
    if [ -d /opt/hubcache/experimenthub ]; then
      mkdir -p "${PWD}/experimenthub-cache"
      cp -r /opt/hubcache/experimenthub/. "${PWD}/experimenthub-cache/"
      export EXPERIMENT_HUB_CACHE="${PWD}/experimenthub-cache"
    fi

    Rscript singler_analysis.R \
      --sce_rds="~{sce_rds}" \
      --sample_name="~{sample_name}" \
      --reference_rds="~{default="" reference_rds}" \
      --reference_dataset="~{reference_dataset}" \
      --reference_ensembl=~{reference_ensembl} \
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
