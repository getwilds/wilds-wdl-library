version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sjl/ww-sjl.wdl" as ww_sjl

workflow sjl_example {
  # Generate synthetic test tile and border points data
  call ww_testdata.generate_sjl_data { }

  # Run SJL tile processing
  call ww_sjl.sjl_tiles { input:
      tile_path = generate_sjl_data.tile_rds,
      border_points_path = generate_sjl_data.border_points_csv,
      year = 2022
  }

  call validate_outputs { input:
      matched_points = sjl_tiles.matched_points,
      missing_points = sjl_tiles.missing_points
  }

  output {
    File matched_points = sjl_tiles.matched_points
    File missing_points = sjl_tiles.missing_points
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate SJL tile processing outputs and generate summary statistics"
    outputs: {
        report: "Validation summary with tile processing statistics"
    }
  }

  parameter_meta {
    matched_points: "RDS file containing points with sunrise/sunset difference values"
    missing_points: "RDS file containing points that could not be matched to border points"
  }

  input {
    File matched_points
    File missing_points
  }

  command <<<
    set -eo pipefail

    echo "SJL Tile Processing Validation Report" > validation_report.txt
    echo "======================================" >> validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    Rscript - <<'RSCRIPT'
      matched <- readRDS("~{matched_points}")
      missing <- readRDS("~{missing_points}")

      sink("validation_report.txt", append = TRUE)
      cat("Matched Points File:", basename("~{matched_points}"), "\n")
      cat("  Total matched points:", nrow(matched), "\n")
      cat("  Columns:", paste(names(matched), collapse = ", "), "\n")
      cat("\n")
      cat("Missing Points File:", basename("~{missing_points}"), "\n")
      cat("  Total unmatched points:", nrow(missing), "\n")
      cat("\n")
      cat("Validation completed successfully.\n")
      sink()
    RSCRIPT
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "rocker/tidyverse:4"
    memory: "4 GB"
    cpu: 1
  }
}
