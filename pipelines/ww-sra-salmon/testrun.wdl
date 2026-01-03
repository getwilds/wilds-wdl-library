version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/consolidate-levels/pipelines/ww-sra-salmon/ww-sra-salmon.wdl" as sra_salmon_workflow

workflow sra_salmon_example {
  # Call testdata workflow to get test transcriptome
  call ww_testdata.download_test_transcriptome { }

  # Call the actual sra_salmon workflow with test data outputs
  # Using SRR3589956: Human (HEK293) RNA-seq for compatibility with GENCODE transcriptome
  # Limiting to 100k reads for fast testing
  call sra_salmon_workflow.sra_salmon { input:
    sra_id_list = ["SRR3589956"],
    transcriptome_fasta = download_test_transcriptome.transcriptome_fasta,
    ncpu = 2,
    memory_gb = 8,
    max_reads = 100000
  }

  output {
    Array[File] salmon_quant_dirs = sra_salmon.salmon_quant_dirs
    File tpm_matrix = sra_salmon.tpm_matrix
    File counts_matrix = sra_salmon.counts_matrix
  }
}
