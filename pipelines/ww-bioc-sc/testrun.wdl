version 1.0

import "../../modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "./ww-bioc-sc.wdl" as bioc_sc_workflow

struct SingleCellSample {
    String name
    File? h5_matrix
    File? raw_h5_matrix
}

workflow bioc_sc_example {
  # Download a raw (unfiltered) 10x Genomics H5 matrix of human PBMCs
  call ww_testdata.download_10x_raw_h5_data { }

  SingleCellSample sample1 = {
    "name": "pbmc_10k_v3",
    "raw_h5_matrix": download_10x_raw_h5_data.raw_h5_matrix
  }

  # Run the full pipeline. Human PBMC data lets emptyDrops call real cells
  # from the raw matrix and lets ww-singler annotate against the default
  # HumanPrimaryCellAtlasData reference with biologically meaningful results.
  call bioc_sc_workflow.bioc_sc { input:
      samples = [sample1],
      mito_pattern = "^MT-",
      cpu_cores = 2,
      memory_gb = 4,
      singler_memory_gb = 8
  }

  output {
    Array[File] sce_rds = bioc_sc.sce_rds
    Array[File?] empty_drops_csv = bioc_sc.empty_drops_csv
    Array[File?] barcode_rank_plot = bioc_sc.barcode_rank_plot
    Array[File] scran_sce_object = bioc_sc.scran_sce_object
    Array[File] scran_qc_plot = bioc_sc.scran_qc_plot
    Array[File] scran_size_factor_plot = bioc_sc.scran_size_factor_plot
    Array[File] scran_mean_variance_plot = bioc_sc.scran_mean_variance_plot
    Array[File] scran_hvg_table = bioc_sc.scran_hvg_table
    Array[File] scater_sce_object = bioc_sc.scater_sce_object
    Array[File] scater_pca_plot = bioc_sc.scater_pca_plot
    Array[File] scater_umap_plot = bioc_sc.scater_umap_plot
    Array[File] scater_qc_plot = bioc_sc.scater_qc_plot
    Array[File] singler_sce_object = bioc_sc.singler_sce_object
    Array[File] singler_cluster_table = bioc_sc.singler_cluster_table
    Array[File] singler_marker_table = bioc_sc.singler_marker_table
    Array[File] singler_prediction_table = bioc_sc.singler_prediction_table
    File validation_report = validate_outputs.report
  }

  call validate_outputs { input:
      singler_sce_object = bioc_sc.singler_sce_object,
      singler_cluster_table = bioc_sc.singler_cluster_table,
      singler_marker_table = bioc_sc.singler_marker_table,
      singler_prediction_table = bioc_sc.singler_prediction_table
  }
}

task validate_outputs {
  meta {
    description: "Validate that all expected ww-bioc-sc output files were generated correctly"
    outputs: {
        report: "Validation report summarizing file checks"
    }
  }

  parameter_meta {
    singler_sce_object: "Array of final SingleCellExperiment RDS objects"
    singler_cluster_table: "Array of per-cell cluster assignment CSVs"
    singler_marker_table: "Array of per-cluster marker gene CSVs"
    singler_prediction_table: "Array of per-cluster SingleR prediction CSVs"
  }

  input {
    Array[File] singler_sce_object
    Array[File] singler_cluster_table
    Array[File] singler_marker_table
    Array[File] singler_prediction_table
  }

  command <<<
    set -eo pipefail

    echo "=== ww-bioc-sc Pipeline Validation Report ===" > validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    for file_path in ~{sep=" " singler_sce_object} ~{sep=" " singler_cluster_table} ~{sep=" " singler_marker_table} ~{sep=" " singler_prediction_table}; do
      if [[ -f "$file_path" && -s "$file_path" ]]; then
        file_size=$(stat -c%s "$file_path")
        echo "$(basename $file_path): PASSED (${file_size} bytes)" >> validation_report.txt
      else
        echo "$(basename $file_path): MISSING OR EMPTY" >> validation_report.txt
        validation_passed=false
      fi
    done

    echo "" >> validation_report.txt
    echo "=== Validation Summary ===" >> validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> validation_report.txt
    else
      echo "Overall Status: FAILED" >> validation_report.txt
      exit 1
    fi

    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "ubuntu:22.04"
    memory: "2 GB"
    cpu: 1
  }
}
