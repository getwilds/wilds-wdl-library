version 1.0

# Import module in question as well as the testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-trimgalore/ww-trimgalore.wdl" as ww_trimgalore
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

# Define data structure for paired-end sample inputs
struct TrimGaloreSample {
    String name
    File r1_fastq
    File r2_fastq
}

#### TEST WORKFLOW DEFINITION ####
# Define test workflow to demonstrate Trim Galore module functionality

workflow trimgalore_example {
  # Auto-download test data for testing purposes
  call ww_testdata.download_fastq_data as download_demo_data { }

  # Create samples array using test data
  Array[TrimGaloreSample] final_samples = [
    {
      "name": "demo_sample",
      "r1_fastq": download_demo_data.r1_fastq,
      "r2_fastq": download_demo_data.r2_fastq
    }
  ]

  # Process each sample with Trim Galore paired-end and single-end trimming
  scatter (sample in final_samples) {
    call ww_trimgalore.trimgalore_paired { input:
        sample_name = sample.name,
        r1_fastq = sample.r1_fastq,
        r2_fastq = sample.r2_fastq,
        cpu_cores = 1,
        memory_gb = 4
    }
    call ww_trimgalore.trimgalore_single { input:
        sample_name = sample.name,
        fastq = sample.r1_fastq,
        cpu_cores = 1,
        memory_gb = 4
    }
  }

  output {
    Array[File] r1_trimmed = trimgalore_paired.r1_trimmed
    Array[File] r2_trimmed = trimgalore_paired.r2_trimmed
    Array[File] paired_r1_reports = trimgalore_paired.r1_report
    Array[File] paired_r2_reports = trimgalore_paired.r2_report
    Array[File] single_trimmed = trimgalore_single.trimmed_fastq
    Array[File] single_reports = trimgalore_single.trimming_report
  }
}
