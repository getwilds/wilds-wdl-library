version 1.0

# Import module in question as well as the testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-star-deseq/modules/ww-rseqc/ww-rseqc.wdl" as ww_rseqc
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bedops/ww-bedops.wdl" as ww_bedops
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####
# Define test workflow to demonstrate module functionality

workflow rseqc_example {
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

  # Convert GTF to BED12 format for RSeQC using ww-bedops
  call ww_bedops.gtf_to_bed {
    input:
      gtf_file = download_ref.gtf
  }

  # Run RSeQC on the test BAM file
  call ww_rseqc.run_rseqc {
    input:
      bam_file = download_bam.bam,
      bam_index = download_bam.bai,
      ref_bed = gtf_to_bed.bed_file,
      sample_name = "test_sample",
      cpu_cores = 2,
      memory_gb = 4
  }

  output {
    File bed_annotation = gtf_to_bed.bed_file
    File read_distribution = run_rseqc.read_distribution
    File gene_body_coverage = run_rseqc.gene_body_coverage
    File infer_experiment = run_rseqc.infer_experiment
    File qc_summary = run_rseqc.rseqc_summary
  }
}
