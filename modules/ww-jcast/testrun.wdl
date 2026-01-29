version 1.0

## Test workflow for the ww-jcast module
## Downloads test data and demonstrates JCAST functionality for alternative splicing proteomics

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jcast/modules/ww-jcast/ww-jcast.wdl" as ww_jcast
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-jcast/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####

workflow jcast_example {
  # Download reference data (GTF and FASTA)
  call ww_testdata.download_ref_data {
    input:
      chromo = "chr1",
      version = "hg38",
      region = "1-5000000"
  }

  # Download rMATS test data from the JCAST repository
  call ww_testdata.download_jcast_test_data { }

  # Run JCAST with test data
  call ww_jcast.jcast {
    input:
      rmats_directory = download_jcast_test_data.rmats_output,
      gtf_file = download_ref_data.gtf,
      genome_fasta = download_ref_data.fasta,
      output_name = "test_jcast",
      min_read_count = 1,
      qvalue_min = 0,
      qvalue_max = 1,
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    File protein_fasta = jcast.output_fasta
    File all_results = jcast.output_directory
  }
}
