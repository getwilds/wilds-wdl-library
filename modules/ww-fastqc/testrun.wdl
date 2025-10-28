version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-fastqc/ww-fastqc.wdl" as ww_fastqc

workflow fastqc_example {
  # Auto-download test data for testing purposes
  call ww_testdata.download_fastq_data as download_demo_data { }

  # Process paired-end sample
  call ww_fastqc.run_fastqc as run_fastqc_paired { input:
      r1_fastq = download_demo_data.r1_fastq,
      r2_fastq = download_demo_data.r2_fastq,
      cpu_cores = 2,
      memory_gb = 4
  }

  # Process single-end sample (R1 only)
  call ww_fastqc.run_fastqc as run_fastqc_single { input:
      r1_fastq = download_demo_data.r1_fastq,
      cpu_cores = 2,
      memory_gb = 4
  }

  output {
    Array[File] fastqc_html_reports = flatten([run_fastqc_paired.html_reports, run_fastqc_single.html_reports])
    Array[File] fastqc_zip_reports = flatten([run_fastqc_paired.zip_reports, run_fastqc_single.zip_reports])
  }
}
