## WILDS WDL Module: ww-seurat
## Description: Module for single-cell RNA-seq analysis using Seurat

version 1.0

task run_seurat {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Perform basic single-cell RNA-seq QC and clustering with Seurat"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-seurat/ww-seurat.wdl"
    outputs: {
        seurat_object: "Processed Seurat RDS object",
        qc_plot: "QC violin plot (nFeature_RNA, nCount_RNA, percent.mt)",
        umap_plot: "UMAP plot showing Louvain clusters",
        heatmap_plot: "Heatmap of top 8 marker genes per cluster",
        cluster_markers: "CSV of top 30 marker genes per cluster from FindAllMarkers"
    }
    topic: "transcriptomics,single_cell"
    species: "human,eukaryote"
    operation: "expression_analysis"
    input_sample_required: "input_h5:gene_expression_matrix:h5"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "seurat_object:gene_expression_matrix:rds_format,qc_plot:quality_control_report:png,umap_plot:plot:png,heatmap_plot:plot:png,cluster_markers:gene_report:csv"
    output_reference: "none"
  }

  parameter_meta {
    input_h5: "Cell Ranger filtered feature-barcode matrix HDF5 file (filtered_feature_bc_matrix.h5)"
    sample_name: "Sample name used for output file prefixes"
    min_cells: "Minimum number of cells a gene must be detected in to be included"
    min_features: "Minimum number of features (genes) a cell must have to be included"
    max_percent_mt: "Maximum percentage of mitochondrial gene counts allowed per cell"
    resolution: "Louvain clustering resolution (higher = more clusters)"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
  }

  input {
    File input_h5
    String sample_name
    Int min_cells = 3
    Int min_features = 200
    Float max_percent_mt = 10.0
    Float resolution = 0.5
    Int memory_gb = 16
    Int cpu_cores = 4
  }

  command <<<
    set -eo pipefail

    Rscript seurat_analysis.R \
      --input_h5="~{input_h5}" \
      --sample_name="~{sample_name}" \
      --min_cells=~{min_cells} \
      --min_features=~{min_features} \
      --max_percent_mt=~{max_percent_mt} \
      --resolution=~{resolution} \
      --ram_gb=~{memory_gb} \
      --output_prefix="~{sample_name}_seurat"
  >>>

  output {
    File seurat_object = "~{sample_name}_seurat.rds"
    File qc_plot = "~{sample_name}_seurat_qc.png"
    File umap_plot = "~{sample_name}_seurat_umap.png"
    File heatmap_plot = "~{sample_name}_seurat_heatmap.png"
    File cluster_markers = "~{sample_name}_seurat_top30_markers.csv"
  }

  runtime {
    docker: "getwilds/seurat:5.2.1"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}
