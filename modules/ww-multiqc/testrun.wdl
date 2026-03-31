version 1.0

# Import module under test and dependencies
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-multiqc/ww-multiqc.wdl" as ww_multiqc
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-fastqc/ww-fastqc.wdl" as ww_fastqc

#### TEST WORKFLOW DEFINITION ####

workflow multiqc_example {
  # Auto-download test data for testing purposes
  call ww_testdata.download_fastq_data as download_demo_data { }

  # Run FastQC on test data to generate QC reports for MultiQC to aggregate
  call ww_fastqc.run_fastqc as run_fastqc { input:
      r1_fastq = download_demo_data.r1_fastq,
      r2_fastq = download_demo_data.r2_fastq,
      cpu_cores = 2,
      memory_gb = 4
  }

  # Run MultiQC to aggregate FastQC results
  call ww_multiqc.run_multiqc { input:
      input_files = run_fastqc.zip_reports,
      report_title = "WILDS Test Data QC Report",
      output_prefix = "wilds_test",
      cpu_cores = 1,
      memory_gb = 4
  }

  output {
    File multiqc_html_report = run_multiqc.html_report
    File multiqc_data = run_multiqc.data_dir
  }
}
