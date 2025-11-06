version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-star-deseq/modules/ww-sra/ww-sra.wdl" as ww_sra
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-star-deseq/vignettes/ww-star-deseq2/ww-star-deseq2.wdl" as star_deseq2_workflow

struct SampleInfo {
    String name
    File r1
    File r2
    String condition
}

struct RefGenome {
    String name
    File fasta
    File gtf
}

workflow star_deseq2_example {
  # Call ww_testdata tasks to get test reference data
  call ww_testdata.download_ref_data { }

  # Download RNA-seq data from SRA (DESeq2 vignette dataset - airway study)
  # Using 2 treated + 2 untreated samples for testing
  # Limiting to 1M reads per sample for memory efficiency on CI runners
  call ww_sra.fastqdump as untreated1 { input: sra_id = "SRR1039509", ncpu = 2, max_reads = 1000000 }
  call ww_sra.fastqdump as untreated2 { input: sra_id = "SRR1039513", ncpu = 2, max_reads = 1000000 }
  call ww_sra.fastqdump as treated1 { input: sra_id = "SRR1039508", ncpu = 2, max_reads = 1000000 }
  call ww_sra.fastqdump as treated2 { input: sra_id = "SRR1039512", ncpu = 2, max_reads = 1000000 }

  # Construct SampleInfo structs from SRA downloads
  SampleInfo sample1 = {
    "name": "untreated_1",
    "r1": untreated1.r1_end,
    "r2": untreated1.r2_end,
    "condition": "untreated"
  }

  SampleInfo sample2 = {
    "name": "untreated_2",
    "r1": untreated2.r1_end,
    "r2": untreated2.r2_end,
    "condition": "untreated"
  }

  SampleInfo sample3 = {
    "name": "treated_1",
    "r1": treated1.r1_end,
    "r2": treated1.r2_end,
    "condition": "treated"
  }

  SampleInfo sample4 = {
    "name": "treated_2",
    "r1": treated2.r1_end,
    "r2": treated2.r2_end,
    "condition": "treated"
  }

  # Construct RefGenome struct
  RefGenome reference = {
    "name": "hg38",
    "fasta": download_ref_data.fasta,
    "gtf": download_ref_data.gtf
  }

  # Run actual STAR-DESeq2 workflow
  call star_deseq2_workflow.star_deseq2 {
    input:
      samples = [sample1, sample2, sample3, sample4],
      reference_genome = reference,
      reference_level = "untreated",
      contrast = "condition,treated,untreated",
      star_cpu = 2,
      star_memory_gb = 6
  }

  output {
    Array[File] star_bam = star_deseq2.star_bam
    Array[File] star_bai = star_deseq2.star_bai
    Array[File] star_gene_counts = star_deseq2.star_gene_counts
    Array[File] star_log_final = star_deseq2.star_log_final
    Array[File] star_log_progress = star_deseq2.star_log_progress
    Array[File] star_log = star_deseq2.star_log
    Array[File] star_sj = star_deseq2.star_sj
    Array[File] rseqc_qc_summary = star_deseq2.rseqc_qc_summary
    File combined_counts_matrix = star_deseq2.combined_counts_matrix
    File sample_metadata = star_deseq2.sample_metadata
    File deseq2_all_results = star_deseq2.deseq2_all_results
    File deseq2_significant_results = star_deseq2.deseq2_significant_results
    File deseq2_normalized_counts = star_deseq2.deseq2_normalized_counts
    File deseq2_pca_plot = star_deseq2.deseq2_pca_plot
    File deseq2_volcano_plot = star_deseq2.deseq2_volcano_plot
    File deseq2_heatmap = star_deseq2.deseq2_heatmap
  }
}
