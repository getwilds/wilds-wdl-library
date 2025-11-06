version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
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
  # Call ww_testdata tasks to get test data
  call ww_testdata.download_ref_data { }

  call ww_testdata.download_fastq_data { }

  # Construct SampleInfo structs from task outputs
  # Creating two samples with different conditions for DESeq2 analysis
  SampleInfo sample1 = {
    "name": "Sample1",
    "r1": download_fastq_data.r1_fastq,
    "r2": download_fastq_data.r2_fastq,
    "condition": "control"
  }

  SampleInfo sample2 = {
    "name": "Sample2",
    "r1": download_fastq_data.r1_fastq,
    "r2": download_fastq_data.r2_fastq,
    "condition": "treatment"
  }

  # Construct RefGenome struct
  RefGenome reference = {
    "name": "chr1",
    "fasta": download_ref_data.fasta,
    "gtf": download_ref_data.gtf
  }

  # Run actual STAR-DESeq2 workflow
  call star_deseq2_workflow.star_deseq2 {
    input:
      samples = [sample1, sample2],
      reference_genome = reference,
      reference_level = "control",
      contrast = "condition,treatment,control",
      star_cpu = 2,
      star_memory_gb = 8
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
