## WILDS WDL module for Geospatial Linkage.
## provides tasks for extracting values from raster data for points/buffers
## and computing aggregated exposure metrics.

version 1.0

task extract_raster {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Extract a value or statistic from a single raster for a set of points or buffers"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-geolink/ww-geolink.wdl"
    outputs: {
      result_csv: "CSV file containing extracted values"
    }
  }

  parameter_meta {
    points_csv: "CSV file containing patient coordinates (columns: patientID, longitude, latitude)"
    raster_tif: "TIF raster file for exposure data"
    buffer_dist: "Buffer distance in meters (0 for point extraction)"
    stat_fun: "Summary function to use (e.g., mean, min, max)"
    output_prefix: "Prefix for the output CSV filename (the raster name and buffer distance are appended)"
  }

  input {
    File points_csv
    File raster_tif
    Float buffer_dist = 0.0
    String stat_fun = "mean"
    String output_prefix = "extract"
  }

  String output_name = "~{output_prefix}_~{basename(raster_tif, '.tif')}_buf~{round(buffer_dist)}m.csv"

  command <<<
    set -eo pipefail
    wget -O extract_val.R https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-exposure/modules/ww-geolink/extract_val.R
    Rscript extract_val.R \
      --points "~{points_csv}" \
      --raster "~{raster_tif}" \
      --buffer ~{buffer_dist} \
      --fun "~{stat_fun}" \
      --out "~{output_name}"
  >>>

  output {
    File result_csv = output_name
  }

  runtime {
    docker: "getwilds/r-spatial:0.1.0"
  }
}

task compute_wildfire_metrics {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Compute annual wildfire exposure metrics and lags from daily extracts"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-geolink/ww-geolink.wdl"
    outputs: {
      annual_metrics: "CSV file containing annual wildfire metrics"
    }
  }

  parameter_meta {
    daily_extracts: "Array of daily exposure CSVs"
    p90_pctile: "Percentile for smoke-wave threshold (default 0.90)"
    min_days: "Minimum days to qualify as a smoke-wave (default 2)"
  }

  input {
    Array[File] daily_extracts
    Float p90_pctile = 0.90
    Int min_days = 2
  }

  command <<<
    set -eo pipefail
    wget -O compute_wf_metrics.R https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-exposure/modules/ww-geolink/compute_wf_metrics.R
    Rscript compute_wf_metrics.R \
      --inputs "~{sep=',' daily_extracts}" \
      --pctile ~{p90_pctile} \
      --min_days ~{min_days} \
      --out "annual_wildfire_metrics.csv"
  >>>

  output {
    File annual_metrics = "annual_wildfire_metrics.csv"
  }

  runtime {
    docker: "getwilds/r-spatial:0.1.0"
  }
}

task join_exposures {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Join fragmented exposure extracts and annual metrics into a final dataset"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-geolink/ww-geolink.wdl"
    outputs: {
      final_dataset: "Final joined exposure dataset CSV"
    }
  }

  parameter_meta {
    patient_base: "Base patient demographics CSV"
    all_extracts: "Array of all extracted exposure CSVs"
    annual_wf_metrics: "Annual wildfire metrics CSV"
  }

  input {
    File patient_base
    Array[File] all_extracts
    File annual_wf_metrics
  }

  command <<<
    set -eo pipefail
    wget -O join_all.R https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-exposure/modules/ww-geolink/join_all.R
    Rscript join_all.R \
      --base "~{patient_base}" \
      --extracts "~{sep=',' all_extracts}" \
      --annuals "~{annual_wf_metrics}" \
      --out "final_dataset.csv"
  >>>

  output {
    File final_dataset = "final_dataset.csv"
  }

  runtime {
    docker: "getwilds/r-spatial:0.1.0"
  }
}
