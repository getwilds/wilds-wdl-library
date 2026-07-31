version 1.0

import "./ww-gffread.wdl" as ww_gffread
import "../ww-testdata/ww-testdata.wdl" as ww_testdata

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

  # Case 3: GFF3-to-GTF conversion, reusing the PAO1 GFF3 fetched above
  call ww_gffread.gff3_to_gtf { input:
      input_gff3 = download_pao1_ref.gff3,
      output_prefix = "pao1_converted"
  }

  output {
    File bacterial_normalized_gtf = normalize_bacterial.normalized_gtf
    File eukaryotic_normalized_gtf = normalize_eukaryotic.normalized_gtf
    File pao1_converted_gtf = gff3_to_gtf.gtf_file
  }
}
