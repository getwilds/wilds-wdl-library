version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-fastq-to-cram/vignettes/ww-fastq-to-cram/ww-fastq-to-cram.wdl" as ww_fastq_to_cram

workflow fastq_to_cram_example {
  # Download test FASTQ data
  call ww_testdata.download_fastq_data as download_fastq_1 { }
  call ww_testdata.download_fastq_data as download_fastq_2 { }

  # Create test batch with two FASTQ groups to test merging functionality
  Array[ww_fastq_to_cram.FastqGroup] test_groups = [
    {
      "group_name": "test_group_1",
      "fastq_r1_locations": [download_fastq_1.r1_fastq],
      "fastq_r2_locations": [download_fastq_1.r2_fastq]
    },
    {
      "group_name": "test_group_2",
      "fastq_r1_locations": [download_fastq_2.r1_fastq],
      "fastq_r2_locations": [download_fastq_2.r2_fastq]
    }
  ]

  Array[ww_fastq_to_cram.SampleData] test_batch = [
    {
      "sample_name": "NA12878_test",
      "library_name": "test_library",
      "sequencing_center": "test_center",
      "fastq_groups": test_groups
    }
  ]

  # Run the fastq-to-cram workflow
  call ww_fastq_to_cram.fastq_to_cram { input:
      batch_info = test_batch,
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    Array[File] unmapped_crams = fastq_to_cram.unmapped_crams
    Array[File] unmapped_cram_indexes = fastq_to_cram.unmapped_cram_indexes
    Array[File] validation_reports = fastq_to_cram.validation_reports
  }
}
