version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/vignettes/ww-fastq-to-cram/ww-fastq-to-cram.wdl" as fastq_to_cram_wf

# Import structs from the workflow
struct FastqGroup {
  String group_name
  Array[File] fastq_r1_locations
  Array[File] fastq_r2_locations
}

struct SampleData {
  String sample_name
  String? library_name
  String? sequencing_center
  Array[FastqGroup] fastq_groups
}

workflow fastq_to_cram_example {
  # Download test FASTQ data
  call ww_testdata.download_fastq_data { }

  # Create test batch with one FASTQ group
  # Note: We use a single group here because the test data doesn't have multiple
  # distinct FASTQ files to merge. In real usage, multiple groups would represent
  # different flowcells, lanes, or sequencing runs for the same sample.
  FastqGroup group1 = object {
    group_name: "test_group_1",
    fastq_r1_locations: [download_fastq_data.r1_fastq],
    fastq_r2_locations: [download_fastq_data.r2_fastq]
  }

  SampleData sample = object {
    sample_name: "NA12878_test",
    library_name: "test_library",
    sequencing_center: "test_center",
    fastq_groups: [group1]
  }

  # Run the fastq-to-cram workflow
  call fastq_to_cram_wf.fastq_to_cram { input:
      batch_info = [sample],
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    Array[File] unmapped_crams = fastq_to_cram.unmapped_crams
    Array[File] validation_reports = fastq_to_cram.validation_reports
  }
}
