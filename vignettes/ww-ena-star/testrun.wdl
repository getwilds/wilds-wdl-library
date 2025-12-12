version 1.0

# import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-ww-ena/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
# import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-ww-ena/vignettes/ww-ena-star/ww-ena-star.wdl" as ena_star_workflow
import "../../modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "ww-ena-star.wdl" as ena_star_workflow

workflow ena_star_example {
  # Call testdata workflow to get test data
  call ww_testdata.download_ref_data { }

  # Call the actual ena_star workflow with test data outputs
  # Using ERR000001: Small test dataset from ENA
  # Note: This is a very small test file to ensure fast execution
  call ena_star_workflow.ena_star { input:
    ena_accession_list = ["ERR10825982"],
    ref_genome = {
      "name": "chr1",
      "fasta": download_ref_data.fasta,
      "gtf": download_ref_data.gtf
    },
    ncpu = 2,
    memory_gb = 8,
    ena_file_format = "READS_FASTQ",
    ena_protocol = "FTP"
  }

  output {
    Array[File] star_bam = ena_star.star_bam
    Array[File] star_bai = ena_star.star_bai
    Array[File] star_gene_counts = ena_star.star_gene_counts
    Array[File] star_log_final = ena_star.star_log_final
    Array[File] star_log_progress = ena_star.star_log_progress
    Array[File] star_log = ena_star.star_log
    Array[File] star_sj = ena_star.star_sj
  }
}
