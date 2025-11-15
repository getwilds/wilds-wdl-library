version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-star-deseq/modules/ww-bedparse/ww-bedparse.wdl" as ww_bedparse
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow bedparse_example {
  # Auto-download test reference data for testing purposes
  call ww_testdata.download_ref_data as download_ref {
    input:
      chromo = "chr1",
      version = "hg38"
  }

  # Convert GTF to BED12 format using bedparse
  call ww_bedparse.gtf2bed {
    input:
      gtf_file = download_ref.gtf
  }

  output {
    File gtf_annotation = download_ref.gtf
    File bed_annotation = gtf2bed.bed_file
  }
}
