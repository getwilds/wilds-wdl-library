## WILDS WDL Module: ww-deseq2
## Description: Module for differential expression analysis using DESeq2
## Author: WILDS Team
## Contact: wilds@fredhutch.org

version 1.0

task combine_count_matrices {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Combine STAR gene count files from multiple samples into a single count matrix"
    outputs: {
        counts_matrix: "Combined matrix of gene-level counts from all samples",
        sample_metadata: "Metadata file containing sample names and conditions"
    }
  }

  parameter_meta {
    gene_count_files: "Array of STAR gene count files (ReadsPerGene.out.tab) from each sample"
    sample_names: "Array of sample names corresponding to the gene_count_files"
    sample_conditions: "Array of experimental conditions corresponding to each sample"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
    count_column: "Column number in STAR count files to extract (2=unstranded, 3=forward strand, 4=reverse strand)"
  }

  input {
    Array[File] gene_count_files
    Array[String] sample_names
    Array[String] sample_conditions
    Int memory_gb = 4
    Int cpu_cores = 1
    Int count_column = 2
        # Column to extract from ReadsPerGene.out.tab files:
        # 2 = unstranded counts
        # 3 = stranded counts, first read forward
        # 4 = stranded counts, first read reverse
  }

  command <<<
    set -eo pipefail

    combine_star_counts.py \
      --input ~{sep=" " gene_count_files} \
      --names ~{sep=" " sample_names} \
      --conditions ~{sep=" " sample_conditions} \
      --output combined_counts_matrix.txt \
      --count_column ~{count_column}
  >>>

  output {
    File counts_matrix = "combined_counts_matrix.txt"
    File sample_metadata = "sample_metadata.txt"
  }

  runtime {
    docker: "getwilds/combine-counts:0.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task run_deseq2 {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Perform differential expression analysis using DESeq2"
    outputs: {
        deseq2_results: "Complete results of DESeq2 differential expression analysis with statistics",
        deseq2_significant: "Filtered results containing only statistically significant differentially expressed genes",
        deseq2_normalized_counts: "Normalized count values produced by DESeq2 for all samples",
        deseq2_pca_plot: "Principal Component Analysis plot showing sample clustering based on expression patterns",
        deseq2_volcano_plot: "Volcano plot showing log fold change vs. statistical significance",
        deseq2_heatmap: "Heatmap visualization of differentially expressed genes across samples"
    }
  }

  parameter_meta {
    counts_matrix: "Combined matrix of gene-level counts from all samples"
    sample_metadata: "Metadata file containing sample information including conditions"
    condition_column: "Column name in the metadata file that contains the experimental condition"
    reference_level: "Reference level for the contrast (typically the control condition)"
    contrast: "DESeq2 contrast string in the format 'condition,treatment,control'"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
  }

  input {
    File counts_matrix
    File sample_metadata
    String condition_column = "condition"
    String reference_level = ""
    String contrast = ""
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail

    deseq2_analysis.R \
      --counts_file="~{counts_matrix}" \
      --metadata_file="~{sample_metadata}" \
      --condition_column="~{condition_column}" \
      --reference_level="~{reference_level}" \
      --contrast="~{contrast}" \
      --output_prefix="deseq2_results"
  >>>

  output {
    File deseq2_results = "deseq2_results_all_genes.csv"
    File deseq2_significant = "deseq2_results_significant.csv"
    File deseq2_normalized_counts = "deseq2_results_normalized_counts.csv"
    File deseq2_pca_plot = "deseq2_results_pca.pdf"
    File deseq2_volcano_plot = "deseq2_results_volcano.pdf"
    File deseq2_heatmap = "deseq2_results_heatmap.pdf"
  }

  runtime {
    docker: "getwilds/deseq2:1.40.2"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}
