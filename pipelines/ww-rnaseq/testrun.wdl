version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl" as ww_sra
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-rnaseq/ww-rnaseq.wdl" as rnaseq_workflow

struct RefGenome {
    String name
    File fasta
    File gtf
}

workflow rnaseq_example {
  # Download reference data: 50Mbp region of chr1 to keep STAR index size manageable
  call ww_testdata.download_ref_data {
    input:
      region = "1-50000000",
      output_name = "chr1_50M"
  }

  # Download RNA-seq data from SRA (DESeq2 vignette dataset - airway study)
  # Using 2 treated + 2 untreated samples (DESeq2 requires >=2 replicates per condition)
  # Using 250K reads per sample to keep runtime and resource usage low for testing
  call ww_sra.fastqdump as untreated1 { input: sra_id = "SRR1039509", ncpu = 2, max_reads = 250000 }
  call ww_sra.fastqdump as untreated2 { input: sra_id = "SRR1039513", ncpu = 2, max_reads = 250000 }
  call ww_sra.fastqdump as treated1 { input: sra_id = "SRR1039508", ncpu = 2, max_reads = 250000 }
  call ww_sra.fastqdump as treated2 { input: sra_id = "SRR1039512", ncpu = 2, max_reads = 250000 }

  # Construct RefGenome struct
  RefGenome reference = {
    "name": "chr1_50M",
    "fasta": download_ref_data.fasta,
    "gtf": download_ref_data.gtf
  }

  # Run the full RNA-seq pipeline with reduced resources for testing
  call rnaseq_workflow.rnaseq {
    input:
      sample_names = ["untreated_1", "untreated_2", "treated_1", "treated_2"],
      r1_fastqs = [untreated1.r1_end, untreated2.r1_end, treated1.r1_end, treated2.r1_end],
      r2_fastqs = [untreated1.r2_end, untreated2.r2_end, treated1.r2_end, treated2.r2_end],
      conditions = ["untreated", "untreated", "treated", "treated"],
      reference_genome = reference,
      reference_level = "untreated",
      contrast = "condition,treated,untreated",
      trim_quality = 20,
      trim_length = 20,
      star_cpu = 2,
      star_memory_gb = 6,
      genome_sa_index_nbases = 10
  }

  output {
    Array[Array[File]] pretrim_fastqc_html = rnaseq.pretrim_fastqc_html
    Array[Array[File]] pretrim_fastqc_zip = rnaseq.pretrim_fastqc_zip
    Array[File] trimgalore_r1_trimmed = rnaseq.trimgalore_r1_trimmed
    Array[File] trimgalore_r2_trimmed = rnaseq.trimgalore_r2_trimmed
    Array[File] trimgalore_r1_report = rnaseq.trimgalore_r1_report
    Array[File] trimgalore_r2_report = rnaseq.trimgalore_r2_report
    Array[Array[File]] posttrim_fastqc_html = rnaseq.posttrim_fastqc_html
    Array[Array[File]] posttrim_fastqc_zip = rnaseq.posttrim_fastqc_zip
    Array[File] star_bam = rnaseq.star_bam
    Array[File] star_bai = rnaseq.star_bai
    Array[File] star_gene_counts = rnaseq.star_gene_counts
    Array[File] star_log_final = rnaseq.star_log_final
    Array[File] star_log_progress = rnaseq.star_log_progress
    Array[File] star_log = rnaseq.star_log
    Array[File] star_sj = rnaseq.star_sj
    Array[File] rseqc_qc_summary = rnaseq.rseqc_qc_summary
    File combined_counts_matrix = rnaseq.combined_counts_matrix
    File sample_metadata = rnaseq.sample_metadata
    File deseq2_all_results = rnaseq.deseq2_all_results
    File deseq2_significant_results = rnaseq.deseq2_significant_results
    File deseq2_normalized_counts = rnaseq.deseq2_normalized_counts
    File deseq2_pca_plot = rnaseq.deseq2_pca_plot
    File deseq2_volcano_plot = rnaseq.deseq2_volcano_plot
    File deseq2_heatmap = rnaseq.deseq2_heatmap
    File multiqc_report = rnaseq.multiqc_report
    File multiqc_data = rnaseq.multiqc_data
  }
}
