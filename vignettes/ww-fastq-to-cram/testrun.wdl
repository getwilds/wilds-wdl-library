version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-fastq-to-cram/vignettes/ww-fastq-to-cram/ww-fastq-to-cram.wdl" as fastq_to_cram_wf

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
  call ww_testdata.download_fastq_data as download_fastq_1 { }
  call ww_testdata.download_fastq_data as download_fastq_2 { }

  # Create test batch with two FASTQ groups to test merging functionality
  FastqGroup group1 = object {
    group_name: "test_group_1",
    fastq_r1_locations: [download_fastq_1.r1_fastq],
    fastq_r2_locations: [download_fastq_1.r2_fastq]
  }

  FastqGroup group2 = object {
    group_name: "test_group_2",
    fastq_r1_locations: [download_fastq_2.r1_fastq],
    fastq_r2_locations: [download_fastq_2.r2_fastq]
  }

  SampleData sample = object {
    sample_name: "NA12878_test",
    library_name: "test_library",
    sequencing_center: "test_center",
    fastq_groups: [group1, group2]
  }

  # Run the fastq-to-cram workflow
  call fastq_to_cram_wf.fastq_to_cram { input:
      batch_info = [sample],
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    Array[File] unmapped_crams = fastq_to_cram.unmapped_crams
    Array[File] unmapped_cram_indexes = fastq_to_cram.unmapped_cram_indexes
    Array[File] validation_reports = fastq_to_cram.validation_reports
  }
}
