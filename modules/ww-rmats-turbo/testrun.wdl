version 1.0

## Test workflow for ww-rmats-turbo module
## Demonstrates the module functionality using test data from ww-testdata.
## NOTE: The test BAM files are WGS data, not RNA-seq, so the splicing analysis
## results will not be biologically meaningful but serve to validate the module
## executes correctly.

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jcast/modules/ww-rmats-turbo/ww-rmats-turbo.wdl" as ww_rmats_turbo
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jcast/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow rmats_turbo_example {
  # Download reference GTF annotation
  call ww_testdata.download_ref_data {
    input:
      chromo = "chr1",
      version = "hg38",
      region = "1-10000000"
  }

  # Download test BAM files for sample group 1
  call ww_testdata.download_bam_data as download_bam_1 {
    input:
      filename = "sample1_rep1.bam"
  }

  call ww_testdata.download_bam_data as download_bam_2 {
    input:
      filename = "sample1_rep2.bam"
  }

  # Download test BAM files for sample group 2
  call ww_testdata.download_bam_data as download_bam_3 {
    input:
      filename = "sample2_rep1.bam"
  }

  call ww_testdata.download_bam_data as download_bam_4 {
    input:
      filename = "sample2_rep2.bam"
  }

  # Run main rMATS-turbo analysis (both prep and post in one step)
  call ww_rmats_turbo.rmats {
    input:
      gtf_file = download_ref_data.gtf,
      sample1_bams = [download_bam_1.bam, download_bam_2.bam],
      sample2_bams = [download_bam_3.bam, download_bam_4.bam],
      read_length = 150,
      read_type = "paired",
      output_name = "test_rmats",
      stat_off = true,
      cpu_cores = 2,
      memory_gb = 8
  }

  # Test prep step separately
  call ww_rmats_turbo.rmats_prep {
    input:
      gtf_file = download_ref_data.gtf,
      sample_bams = [download_bam_1.bam, download_bam_2.bam],
      read_length = 150,
      read_type = "paired",
      output_name = "test_prep",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Test stat step on existing output
  call ww_rmats_turbo.rmats_stat {
    input:
      gtf_file = download_ref_data.gtf,
      existing_output = rmats.output_directory,
      read_length = 150,
      read_type = "paired",
      output_name = "test_stat",
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    File rmats_results = rmats.output_directory
    File prep_results = rmats_prep.prep_output
    File stat_results = rmats_stat.output_directory
  }
}
