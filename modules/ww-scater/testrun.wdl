version 1.0

import "../ww-testdata/ww-testdata.wdl" as ww_testdata
import "../ww-dropletutils/ww-dropletutils.wdl" as ww_dropletutils
import "../ww-scran/ww-scran.wdl" as ww_scran
import "./ww-scater.wdl" as ww_scater

workflow scater_example {
  # Download 10x Genomics H5 test data
  call ww_testdata.download_10x_h5_data { }

  # Load the filtered matrix into a SingleCellExperiment (upstream ww-dropletutils step)
  call ww_dropletutils.read10x_counts { input:
      h5_matrix = download_10x_h5_data.h5_matrix,
      sample_name = "2500_Wistar_Rat_PBMCs_Singleplex"
  }

  # Normalize with scran (upstream ww-scran step)
  call ww_scran.run_scran { input:
      sce_rds = read10x_counts.sce_rds,
      sample_name = "2500_Wistar_Rat_PBMCs_Singleplex",
      mito_pattern = "^Mt-"
  }

  # Reuse ww-scran's HVG selection instead of recomputing it
  call ww_scater.run_scater { input:
      sce_rds = run_scran.sce_object,
      hvg_list = run_scran.hvg_list,
      sample_name = "2500_Wistar_Rat_PBMCs_Singleplex"
  }

  call validate_outputs { input:
      sce_object = run_scater.sce_object,
      pca_plot = run_scater.pca_plot,
      umap_plot = run_scater.umap_plot,
      qc_plot = run_scater.qc_plot
  }

  output {
    File sce_object = run_scater.sce_object
    File pca_plot = run_scater.pca_plot
    File umap_plot = run_scater.umap_plot
    File qc_plot = run_scater.qc_plot
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that all expected scater output files were generated correctly"
    outputs: {
        report: "Validation report summarizing file checks"
    }
  }

  parameter_meta {
    sce_object: "SingleCellExperiment RDS object with PCA and UMAP embeddings"
    pca_plot: "PCA scatter plot PNG"
    umap_plot: "UMAP scatter plot PNG"
    qc_plot: "Per-cell QC scatter plot PNG"
  }

  input {
    File sce_object
    File pca_plot
    File umap_plot
    File qc_plot
  }

  command <<<
    set -eo pipefail

    echo "=== scater Analysis Validation Report ===" > validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    for file_path in "~{sce_object}" "~{pca_plot}" "~{umap_plot}" "~{qc_plot}"; do
      if [[ -f "$file_path" && -s "$file_path" ]]; then
        file_size=$(stat -c%s "$file_path")
        echo "$(basename $file_path): PASSED (${file_size} bytes)" >> validation_report.txt
      else
        echo "$(basename $file_path): MISSING OR EMPTY" >> validation_report.txt
        validation_passed=false
      fi
    done

    # SCE RDS object should be at least 1MB for test data, probably more
    min_rds_size=1048576
    rds_size=$(stat -c%s "~{sce_object}")
    if [[ "$rds_size" -lt "$min_rds_size" ]]; then
      echo "sce_object size check: FAILED (${rds_size} bytes, expected >= ${min_rds_size})" >> validation_report.txt
      validation_passed=false
    else
      echo "sce_object size check: PASSED" >> validation_report.txt
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
