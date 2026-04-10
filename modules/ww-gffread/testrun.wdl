version 1.0

# Import the module under test and ww-testdata for auto-provisioned test inputs.
# TODO: switch both imports to refs/heads/main before merging
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/rnaseq-gtf-processing/modules/ww-gffread/ww-gffread.wdl" as ww_gffread
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/rnaseq-gtf-processing/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####
# Exercises normalize_gtf on two very different inputs:
#   1. The NCBI PAO1 GTF (bacterial layout: mostly CDS, ~100 exon rows).
#      Expected behavior: normalize_gtf synthesizes ~5500 exon features.
#   2. The Ensembl human chromosome 15 GTF (well-formed eukaryotic layout).
#      Expected behavior: pass-through — exon features already exist and
#      the feature count is largely unchanged.

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
