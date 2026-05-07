version 1.0

# Import module in question as well as the testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-proseq/modules/ww-fastp/ww-fastp.wdl" as ww_fastp
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

# Define data structure for paired-end sample inputs
struct FastpSample {
    String name
    File r1_fastq
    File r2_fastq
}

#### TEST WORKFLOW DEFINITION ####
# Define test workflow to demonstrate fastp module functionality

workflow fastp_example {
  # Auto-download test data for testing purposes
  call ww_testdata.download_fastq_data as download_demo_data { }

  # Create samples array using test data
  Array[FastpSample] final_samples = [
    {
      "name": "demo_sample",
      "r1_fastq": download_demo_data.r1_fastq,
      "r2_fastq": download_demo_data.r2_fastq
    }
  ]

  # Process each sample with fastp paired-end and single-end trimming
  scatter (sample in final_samples) {
    call ww_fastp.fastp_paired { input:
        sample_name = sample.name,
        r1_fastq = sample.r1_fastq,
        r2_fastq = sample.r2_fastq,
        cpu_cores = 1,
        memory_gb = 4
    }
    # Exercise the UMI path on the same sample (e.g. PRO-seq with UMIs on both ends).
    # fastp moves the 6-mer UMI from each mate into the read name as `:UMI1_UMI2`.
    call ww_fastp.fastp_paired as fastp_paired_umi { input:
        sample_name = sample.name + "_umi",
        r1_fastq = sample.r1_fastq,
        r2_fastq = sample.r2_fastq,
        umi_loc = "per_read",
        umi_len = 6,
        cpu_cores = 1,
        memory_gb = 4
    }
    call ww_fastp.fastp_single { input:
        sample_name = sample.name,
        fastq = sample.r1_fastq,
        cpu_cores = 1,
        memory_gb = 4
    }
  }

  output {
    Array[File] r1_trimmed = fastp_paired.r1_trimmed
    Array[File] r2_trimmed = fastp_paired.r2_trimmed
    Array[File] paired_html_reports = fastp_paired.html_report
    Array[File] paired_json_reports = fastp_paired.json_report
    Array[File] r1_trimmed_umi = fastp_paired_umi.r1_trimmed
    Array[File] r2_trimmed_umi = fastp_paired_umi.r2_trimmed
    Array[File] paired_umi_html_reports = fastp_paired_umi.html_report
    Array[File] paired_umi_json_reports = fastp_paired_umi.json_report
    Array[File] single_trimmed = fastp_single.trimmed_fastq
    Array[File] single_html_reports = fastp_single.html_report
    Array[File] single_json_reports = fastp_single.json_report
  }
}
