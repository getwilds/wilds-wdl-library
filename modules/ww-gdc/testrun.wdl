version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gdc/ww-gdc.wdl" as ww_gdc
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow gdc_client_example {
  meta {
    description: "Test workflow demonstrating GDC Data Transfer Tool module functionality with open-access data"
    outputs: {
        uuid_downloaded_files: "Array of files downloaded using UUIDs",
        uuid_download_log: "Log file for UUID download",
        manifest_downloaded_files: "Array of files downloaded using manifest",
        manifest_download_log: "Log file for manifest download"
    }
  }

  # Example file UUIDs for small open-access TCGA files
  # These are small annotation/metadata files that don't require authentication
  Array[String] test_file_uuids = [
    "6e811713-17b0-4413-a756-af178269824f",
    "4e89ba70-022e-48a3-a8f8-04f5720fb2d0",
    "52fd584b-9ca3-4f7e-bdd5-fd9dce3d630b"
  ]

  # Create test manifest for manifest-based download
  call ww_testdata.create_gdc_manifest

  # Test 1: Download files by UUID (no token needed for open-access data)
  call ww_gdc.download_by_uuids as download_by_uuid {
    input:
      file_uuids = test_file_uuids,
      n_processes = 2,
      cpu_cores = 2,
      memory_gb = 4
  }

  # Test 2: Download using manifest (no token needed for open-access data)
  call ww_gdc.download_by_manifest as download_by_manifest {
    input:
      manifest_file = create_gdc_manifest.manifest,
      n_processes = 4,
      cpu_cores = 2,
      memory_gb = 4
  }

  output {
    Array[File] uuid_downloaded_files = download_by_uuid.downloaded_files
    File uuid_download_log = download_by_uuid.download_log
    Array[File] manifest_downloaded_files = download_by_manifest.downloaded_files
    File manifest_download_log = download_by_manifest.download_log
  }
}
