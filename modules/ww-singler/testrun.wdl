version 1.0

import "../ww-testdata/ww-testdata.wdl" as ww_testdata
import "../ww-dropletutils/ww-dropletutils.wdl" as ww_dropletutils
import "../ww-scran/ww-scran.wdl" as ww_scran
import "../ww-scater/ww-scater.wdl" as ww_scater
import "./ww-singler.wdl" as ww_singler

workflow singler_example {
  # Download a raw (unfiltered) 10x Genomics H5 matrix of human PBMCs
  call ww_testdata.download_10x_raw_h5_data { }

  # Call real cells from the raw matrix (upstream ww-dropletutils step)
  call ww_dropletutils.empty_drops_filter { input:
      raw_h5_matrix = download_10x_raw_h5_data.raw_h5_matrix,
      sample_name = "pbmc_10k_v3"
  }

  # Normalize with scran (upstream ww-scran step)
  call ww_scran.run_scran { input:
      sce_rds = empty_drops_filter.filtered_sce_rds,
      sample_name = "pbmc_10k_v3",
      mito_pattern = "^MT-"
  }

  # Dimensionality reduction with scater (upstream ww-scater step)
  call ww_scater.run_scater { input:
      sce_rds = run_scran.sce_object,
      sample_name = "pbmc_10k_v3"
  }

  call ww_singler.run_singler { input:
      sce_rds = run_scater.sce_object,
      sample_name = "pbmc_10k_v3"
  }

  call validate_outputs { input:
      sce_object = run_singler.sce_object,
      cluster_table = run_singler.cluster_table,
      marker_table = run_singler.marker_table,
      prediction_table = run_singler.prediction_table
  }

  output {
    File sce_object = run_singler.sce_object
    File cluster_table = run_singler.cluster_table
    File marker_table = run_singler.marker_table
    File prediction_table = run_singler.prediction_table
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that all expected singler output files were generated correctly"
    outputs: {
        report: "Validation report summarizing file checks"
    }
  }

  parameter_meta {
    sce_object: "SingleCellExperiment RDS object with cluster assignments and predicted cell types"
    cluster_table: "Per-cell cluster assignment CSV"
    marker_table: "Per-cluster marker gene CSV"
    prediction_table: "Per-cluster SingleR prediction CSV"
  }

  input {
    File sce_object
    File cluster_table
    File marker_table
    File prediction_table
  }

  command <<<
    set -eo pipefail

    echo "=== singler Analysis Validation Report ===" > validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    for file_path in "~{sce_object}" "~{cluster_table}" "~{marker_table}" "~{prediction_table}"; do
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
