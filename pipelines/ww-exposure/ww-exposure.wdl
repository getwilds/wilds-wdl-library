## WILDS WDL pipeline for Environmental Exposure Linkage.
## Parallelizes the extraction of exposure values from raster data
## and computes aggregated health metrics.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-exposure/modules/ww-geolink/ww-geolink.wdl" as ww_geolink
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-exposure/modules/ww-aws-sso/ww-aws-sso.wdl" as aws_sso

workflow exposure_linkage {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Parallel pipeline to link patient coordinates to environmental exposure rasters and compute annual metrics"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-exposure/ww-exposure.wdl"
    outputs: {
      final_dataset: "The final joined dataset containingall exposure measures and metrics"
    }
  }

  parameter_meta {
    patient_csv: "Local path to patient coordinates CSV (optional if S3 URI provided)"
    patient_csv_s3: "S3 URI of the patient coordinates CSV (optional if local file provided)"
    rasters_local: "Local array of TIF files (optional if S3 URIs provided)"
    rasters_manifest_s3: "Array of S3 URIs for the TIF files (optional if local files provided)"
    buffer_distances: "Array of buffer radii in meters (e.g., [0.0, 1000.0, 5000.0])"
    aws_config: "AWS config file for S3 access"
    aws_creds: "AWS credentials file for S3 access"
  }

  input {
    File? patient_csv,
    String? patient_csv_s3,
    Array[File]? rasters_local,
    Array[String]? rasters_manifest_s3,
    Array[Float] buffer_distances = [0.0, 1000.0, 5000.0, 10000.0],
    File? aws_config,
    File? aws_creds
  }

  # 1. Optional S3 Download for base patient file
  call aws_sso.s3_download_file {
    input:
      s3_uri = patient_csv_s3,
      aws_config_file = aws_config,
      aws_credentials_file = aws_creds
  }

  # Normalize patient file input
  File final_patient_file = select_first([patient_csv, s3_download_file.downloaded_file])

  # 2. Optional S3 Download for rasters
  scatter (s3_uri in select_first([rasters_manifest_s3, []])) {
    call aws_sso.s3_download_file {
      input:
        s3_uri = s3_uri,
        aws_config_file = aws_config,
        aws_credentials_file = aws_creds
    }
  }

  # Normalize raster manifest input
  Array[File] final_rasters = select_first([rasters_local, s3_download_file.downloaded_file])

  # 3. Parallel Extraction: one task per raster per buffer distance
  scatter (raster in final_rasters) {
    scatter (dist in buffer_distances) {
      call ww_geolink.extract_raster {
        input:
          points_csv = final_patient_file,
          raster_tif = raster,
          buffer_dist = dist
      }
    }
  }

  # 4. Compute complex wildfire metrics
  call ww_geolink.compute_wildfire_metrics {
    input:
      daily_extracts = extract_raster.result_csv
  }

  # 5. Aggregate all fragmented outputs into the final dataset
  call ww_geolink.join_exposures {
    input:
      patient_base = final_patient_file,
      all_extracts = extract_raster.result_csv,
      annual_wf_metrics = compute_wildfire_metrics.annual_metrics
  }

  output {
    File final_dataset = join_exposures.final_dataset
  }
}
