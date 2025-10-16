version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/vignettes/ww-sra-star/ww-sra-star.wdl" as sra_star_workflow

workflow sra_star_example {
  # Call testdata workflow to get test data
  call ww_testdata.download_ref_data { }

  # Call the actual sra_star workflow with test data outputs
  call sra_star_workflow.sra_star { input:
    sra_id_list = ["SRR13008264"],
    ref_genome = {
      "name": "chr1",
      "fasta": download_ref_data.fasta,
      "gtf": download_ref_data.gtf
    },
    ncpu = 2,
    memory_gb = 8
  }

  output {
    Array[File] star_bam = sra_star.star_bam
    Array[File] star_bai = sra_star.star_bai
    Array[File] star_gene_counts = sra_star.star_gene_counts
    Array[File] star_log_final = sra_star.star_log_final
    Array[File] star_log_progress = sra_star.star_log_progress
    Array[File] star_log = sra_star.star_log
    Array[File] star_sj = sra_star.star_sj
    File validation_report = sra_star.validation_report
  }
}
