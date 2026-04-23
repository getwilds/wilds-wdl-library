## WILDS WDL module for SJL (Solar Jetlag) tile processing.
## Calculates sunrise/sunset times and sun time differences for geographic tiles
## using NOAA solar calculator variables, as part of the SJL model pipeline.
##
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task sjl_tiles {
  meta {
    author: "Caroline Nondin"
    email: "cnondin@fredhutch.org"
    description: "Calculate sunrise/sunset times and sun time differences for a geographic tile as part of the SJL model pipeline"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sjl/ww-sjl.wdl"
    outputs: {
      matched_points: "RDS file containing points with sunrise/sunset difference values",
      missing_points: "RDS file containing points that could not be matched to border points (may be empty)"
    }
    topic: "public_health_and_epidemiology"
    species: "any"
    operation: "statistical_calculation"
    input_sample_required: "tile_path:accession:binary_format,border_points_path:accession:csv"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "matched_points:report:binary_format,missing_points:report:binary_format"
    output_reference: "none"
  }

  parameter_meta {
    tile_path: "Path to input tile .rds file"
    border_points_path: "Path to border points .csv file containing timezone boundary data"
    year: "Year for solar calculations (e.g. 2022)"
    matched_prefix: "Filename prefix for matched results output (default: 'matched_')"
    missing_prefix: "Filename prefix for missing results output (default: 'missing_')"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File tile_path
    File border_points_path
    Int year
    String matched_prefix = "matched_"
    String missing_prefix = "missing_"
    Int cpu_cores = 1
    Int memory_gb = 8
  }

  String tile_basename = basename(tile_path, ".rds")

  command <<<
    set -eo pipefail

    # Pull sjl_tiles script from GitHub
    # NOTE: For reproducibility in production workflows, replace the branch reference
    # (e.g., "refs/heads/main") with a specific commit hash (e.g., "abc1234...")
    wget -q "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sjl/sjl_tiles.R" \
      -O sjl_tiles.R

    Rscript sjl_tiles.R \
      --tile_path "~{tile_path}" \
      --border_points_path "~{border_points_path}" \
      --year ~{year} \
      --matched_prefix "~{matched_prefix}" \
      --missing_prefix "~{missing_prefix}"
  >>>

  output {
    File matched_points = "~{matched_prefix}~{tile_basename}.rds"
    File missing_points = "~{missing_prefix}~{tile_basename}.rds"
  }

  runtime {
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
    docker: "getwilds/r-utils:0.1.0"
  }
}
