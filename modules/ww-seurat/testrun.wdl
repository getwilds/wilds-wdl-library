version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-seurat/ww-seurat.wdl" as ww_seurat

workflow seurat_example {
  # Download 10x Genomics H5 test data
  call ww_testdata.download_10x_h5_data { }

  call ww_seurat.run_seurat { input:
      input_h5 = download_10x_h5_data.h5_matrix,
      sample_name = "2500_Wistar_Rat_PBMCs_Singleplex"
  }

  call validate_outputs { input:
      seurat_object = run_seurat.seurat_object,
      qc_plot = run_seurat.qc_plot,
      umap_plot = run_seurat.umap_plot,
      heatmap_plot = run_seurat.heatmap_plot,
      cluster_markers = run_seurat.cluster_markers
  }

  output {
    File seurat_object = run_seurat.seurat_object
    File qc_plot = run_seurat.qc_plot
    File umap_plot = run_seurat.umap_plot
    File heatmap_plot = run_seurat.heatmap_plot
    File cluster_markers = run_seurat.cluster_markers
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that all expected Seurat output files were generated correctly"
    outputs: {
        report: "Validation report summarizing file checks"
    }
  }

  parameter_meta {
    seurat_object: "Processed Seurat RDS object"
    qc_plot: "QC violin plot PNG"
    umap_plot: "UMAP plot PNG"
    heatmap_plot: "Heatmap plot PNG"
    cluster_markers: "CSV of top marker genes per cluster"
  }

  input {
    File seurat_object
    File qc_plot
    File umap_plot
    File heatmap_plot
    File cluster_markers
  }

  command <<<
    set -eo pipefail

    echo "=== Seurat Analysis Validation Report ===" > validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    for file_path in "~{seurat_object}" "~{qc_plot}" "~{umap_plot}" "~{heatmap_plot}" "~{cluster_markers}"; do
      if [[ -f "$file_path" && -s "$file_path" ]]; then
        file_size=$(stat -c%s "$file_path")
        echo "$(basename $file_path): PASSED (${file_size} bytes)" >> validation_report.txt
      else
        echo "$(basename $file_path): MISSING OR EMPTY" >> validation_report.txt
        validation_passed=false
      fi
    done

    # Seurat object should be at least 2MB for test data, probably more
    min_rds_size=2097152
    rds_size=$(stat -c%s "~{seurat_object}")
    if [[ "$rds_size" -lt "$min_rds_size" ]]; then
      echo "seurat_object size check: FAILED (${rds_size} bytes, expected >= ${min_rds_size})" >> validation_report.txt
      validation_passed=false
    else
      echo "seurat_object size check: PASSED" >> validation_report.txt
    fi

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
