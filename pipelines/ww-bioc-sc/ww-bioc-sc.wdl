version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-dropletutils/ww-dropletutils.wdl" as dropletutils_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-scran/ww-scran.wdl" as scran_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-scater/ww-scater.wdl" as scater_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-singler/ww-singler.wdl" as singler_tasks

struct SingleCellSample {
    String name
    File? h5_matrix
    File? raw_h5_matrix
}

workflow bioc_sc {
  meta {
    author: [
        {
            name: "Taylor Firman",
            email: "tfirman@fredhutch.org"
        }
    ]
    description: "WDL pipeline for standard Bioconductor/OSCA single-cell RNA-seq post-processing: load matrix, QC filter, normalize, select highly variable genes, run dimensionality reduction, cluster, identify marker genes, and annotate cell types"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-bioc-sc/ww-bioc-sc.wdl"
    outputs: {
        sce_rds: "Array of loaded SingleCellExperiment RDS objects for each sample, prior to QC/normalization",
        empty_drops_csv: "Array of per-barcode emptyDrops statistics for samples provided as a raw matrix (empty for samples provided pre-filtered)",
        barcode_rank_plot: "Array of barcode-rank QC plots for samples provided as a raw matrix (empty for samples provided pre-filtered)",
        scran_sce_object: "Array of normalized SingleCellExperiment RDS objects with size factors and PCA for each sample",
        scran_qc_plot: "Array of per-cell QC scatter plots from scran for each sample",
        scran_size_factor_plot: "Array of size factor histograms for each sample",
        scran_mean_variance_plot: "Array of mean-variance trend plots for each sample",
        scran_hvg_table: "Array of CSVs of genes ranked by biological variance for each sample",
        scater_sce_object: "Array of SingleCellExperiment RDS objects with PCA and UMAP embeddings for each sample",
        scater_pca_plot: "Array of PCA scatter plots for each sample",
        scater_umap_plot: "Array of UMAP scatter plots for each sample",
        scater_qc_plot: "Array of per-cell QC scatter plots from scater for each sample",
        singler_sce_object: "Array of final SingleCellExperiment RDS objects with cluster assignments and predicted cell types for each sample",
        singler_cluster_table: "Array of CSVs mapping cell barcodes to clusters for each sample",
        singler_marker_table: "Array of CSVs of top marker genes per cluster for each sample",
        singler_prediction_table: "Array of CSVs of SingleR cell-type predictions for each sample"
    }
  }

  parameter_meta {
    samples: "List of single-cell sample objects, each with a name and either a pre-filtered 10x feature-barcode matrix (h5_matrix) or a raw, unfiltered one (raw_h5_matrix). Exactly one of the two must be provided per sample; if raw_h5_matrix is given, real cells are called with emptyDrops before QC"
    mito_pattern: "Regex pattern identifying mitochondrial gene symbols, passed to ww-scran"
    nmads: "Number of median absolute deviations from the median used to flag low-quality cells, passed to ww-scran"
    n_hvgs: "Number of top highly variable genes to select in ww-scran and use for PCA in ww-scater"
    n_pcs: "Number of principal components to compute in ww-scater"
    reference_dataset: "Name of the celldex reference dataset function to fetch for cell-type annotation, passed to ww-singler"
    reference_ensembl: "Fetch the celldex reference with Ensembl gene IDs instead of gene symbols, passed to ww-singler"
    label_column: "Column name in the reference object's colData containing cell-type labels, passed to ww-singler"
    cpu_cores: "Number of CPU cores allocated per task"
    memory_gb: "Memory allocated in GB for ww-dropletutils, ww-scran, and ww-scater tasks"
    singler_memory_gb: "Memory allocated in GB for ww-singler tasks (higher default since celldex reference loading and clustering need more headroom)"
  }

  input {
    Array[SingleCellSample] samples
    String mito_pattern = "^MT-"
    Float nmads = 3.0
    Int n_hvgs = 2000
    Int n_pcs = 50
    String reference_dataset = "HumanPrimaryCellAtlasData"
    Boolean reference_ensembl = true
    String label_column = "label.main"
    Int cpu_cores = 2
    Int memory_gb = 8
    Int singler_memory_gb = 16
  }

  scatter (sample in samples) {
    # Pre-filtered matrix path
    if (defined(sample.h5_matrix)) {
      call dropletutils_tasks.read10x_counts { input:
          h5_matrix = select_first([sample.h5_matrix]),
          sample_name = sample.name,
          cpu_cores = cpu_cores,
          memory_gb = memory_gb
      }
    }

    # Raw (unfiltered) matrix path: call real cells with emptyDrops first
    if (defined(sample.raw_h5_matrix)) {
      call dropletutils_tasks.empty_drops_filter { input:
          raw_h5_matrix = select_first([sample.raw_h5_matrix]),
          sample_name = sample.name,
          cpu_cores = cpu_cores,
          memory_gb = memory_gb
      }
    }

    File sample_sce_rds = select_first([read10x_counts.sce_rds, empty_drops_filter.filtered_sce_rds])

    call scran_tasks.run_scran { input:
        sce_rds = sample_sce_rds,
        sample_name = sample.name,
        mito_pattern = mito_pattern,
        nmads = nmads,
        n_hvgs = n_hvgs,
        cpu_cores = cpu_cores,
        memory_gb = memory_gb
    }

    call scater_tasks.run_scater { input:
        sce_rds = run_scran.sce_object,
        hvg_list = run_scran.hvg_list,
        sample_name = sample.name,
        n_pcs = n_pcs,
        n_hvgs = n_hvgs,
        cpu_cores = cpu_cores,
        memory_gb = memory_gb
    }

    call singler_tasks.run_singler { input:
        sce_rds = run_scater.sce_object,
        sample_name = sample.name,
        reference_dataset = reference_dataset,
        reference_ensembl = reference_ensembl,
        label_column = label_column,
        cpu_cores = cpu_cores,
        memory_gb = singler_memory_gb
    }
  }

  output {
    Array[File] sce_rds = sample_sce_rds
    Array[File?] empty_drops_csv = empty_drops_filter.empty_drops_csv
    Array[File?] barcode_rank_plot = empty_drops_filter.barcode_rank_pdf
    Array[File] scran_sce_object = run_scran.sce_object
    Array[File] scran_qc_plot = run_scran.qc_plot
    Array[File] scran_size_factor_plot = run_scran.size_factor_plot
    Array[File] scran_mean_variance_plot = run_scran.mean_variance_plot
    Array[File] scran_hvg_table = run_scran.hvg_table
    Array[File] scater_sce_object = run_scater.sce_object
    Array[File] scater_pca_plot = run_scater.pca_plot
    Array[File] scater_umap_plot = run_scater.umap_plot
    Array[File] scater_qc_plot = run_scater.qc_plot
    Array[File] singler_sce_object = run_singler.sce_object
    Array[File] singler_cluster_table = run_singler.cluster_table
    Array[File] singler_marker_table = run_singler.marker_table
    Array[File] singler_prediction_table = run_singler.prediction_table
  }
}
