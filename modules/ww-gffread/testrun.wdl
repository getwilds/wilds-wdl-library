version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/rnaseq-gtf-processing/modules/ww-gffread/ww-gffread.wdl" as ww_gffread
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/rnaseq-gtf-processing/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow gffread_example {
  # Case 1: Bacterial NCBI GTF (the primary use case for this module)
  call ww_testdata.download_pao1_ref { }

  call ww_gffread.normalize_gtf as normalize_bacterial { input:
      input_gtf = download_pao1_ref.gtf,
      output_prefix = "pao1_normalized"
  }

  # Case 2: Eukaryotic Ensembl GTF (pass-through sanity check)
  call ww_testdata.download_jcast_test_data { }

  call ww_gffread.normalize_gtf as normalize_eukaryotic { input:
      input_gtf = download_jcast_test_data.gtf_file,
      output_prefix = "ensembl_chr15_normalized"
  }

  output {
    File bacterial_normalized_gtf = normalize_bacterial.normalized_gtf
    File bacterial_feature_counts = normalize_bacterial.feature_counts
    File eukaryotic_normalized_gtf = normalize_eukaryotic.normalized_gtf
    File eukaryotic_feature_counts = normalize_eukaryotic.feature_counts
  }
}
