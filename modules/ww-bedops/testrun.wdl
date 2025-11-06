version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-star-deseq/modules/ww-bedops/ww-bedops.wdl" as ww_bedops
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow bedops_example {
  # Auto-download test reference data for testing purposes
  call ww_testdata.download_ref_data as download_ref {
    input:
      chromo = "chr1",
      version = "hg38"
  }

  # Convert GTF to BED12 format
  call ww_bedops.gtf_to_bed {
    input:
      gtf_file = download_ref.gtf
  }

  output {
    File gtf_annotation = download_ref.gtf
    File bed_annotation = gtf_to_bed.bed_file
  }
}
