version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-star-deseq/modules/ww-rnaseqc/ww-rnaseqc.wdl" as ww_rnaseqc
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow rnaseqc_example {
  # Auto-download test reference data for testing purposes
  call ww_testdata.download_ref_data as download_ref {
    input:
      chromo = "chr1",
      version = "hg38"
  }

  # Auto-download test BAM data
  call ww_testdata.download_bam_data as download_bam {
    input:
      filename = "test_sample.bam"
  }

  # Collapse GTF for RNA-SeQC compatibility
  call ww_rnaseqc.collapse_gtf {
    input:
      reference_gtf = download_ref.gtf
  }

  # Run RNA-SeQC on the test BAM file
  call ww_rnaseqc.run_rnaseqc {
    input:
      bam_file = download_bam.bam,
      bam_index = download_bam.bai,
      ref_gtf = collapse_gtf.collapsed_gtf,
      sample_name = "test_sample",
      cpu_cores = 2,
      memory_gb = 4
  }

  output {
    File qc_metrics = run_rnaseqc.rnaseqc_metrics
  }
}
