## WILDS WDL Module: ww-edger
## Description: Module for differential expression analysis using edgeR
## Author: Taylor Firman
## Contact: tfirman@fredhutch.org

version 1.0

task run_edger {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Perform differential expression analysis using edgeR's quasi-likelihood pipeline"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-edger/ww-edger.wdl"
    outputs: {
        edger_results: "Complete edgeR differential expression results with statistics for all genes",
        edger_significant: "Filtered results containing only statistically significant differentially expressed genes",
        edger_normalized_counts: "TMM-normalized CPM values for all samples",
        edger_md_plot: "MD plot showing log fold change vs. mean log CPM",
        edger_volcano_plot: "Volcano plot showing log fold change vs. statistical significance",
        edger_heatmap: "Heatmap visualization of differentially expressed genes across samples",
        edger_bcv_plot: "Biological coefficient of variation plot showing dispersion estimates"
    }
    topic: "transcriptomics,gene_expression"
    species: "human,eukaryote,prokaryote,virus"
    operation: "statistical_calculation"
    input_sample_required: "counts_matrix:gene_expression_matrix:matrix,sample_metadata:report:textual_format"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "edger_results:report:csv,edger_significant:report:csv,edger_normalized_counts:report:csv,edger_md_plot:plot:pdf,edger_volcano_plot:plot:pdf,edger_heatmap:plot:pdf"
    output_reference: "none"
  }

  parameter_meta {
    counts_matrix: "Combined matrix of gene-level counts (genes as rows, samples as columns, tab-delimited with header)"
    sample_metadata: "Metadata file containing sample information including experimental conditions (tab-delimited with header)"
    condition_column: "Column name in the metadata file that contains the experimental condition"
    reference_level: "Reference level for the contrast (typically the control condition; empty = first alphabetically)"
    contrast: "Contrast string in the format 'treatment,control' to specify comparison (empty = infer from condition_column)"
    fdr_threshold: "FDR threshold for calling significant genes"
    lfc_threshold: "Minimum absolute log2 fold change threshold for calling significant genes"
    min_count: "Minimum count threshold used by filterByExpr for low-count gene filtering"
    min_total_count: "Minimum total count across all samples used by filterByExpr"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
    docker_image: "Docker image to use for this task"
  }

  input {
    File counts_matrix
    File sample_metadata
    String condition_column = "condition"
    String reference_level = ""
    String contrast = ""
    Float fdr_threshold = 0.05
    Float lfc_threshold = 1.0
    Int min_count = 10
    Int min_total_count = 15
    Int memory_gb = 8
    Int cpu_cores = 2
    String docker_image = "getwilds/edger:4.10.0"
  }

  command <<<
    set -eo pipefail

    curl -so edger_analysis.R \
      "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-edger-module/modules/ww-edger/edger_analysis.R"

    Rscript edger_analysis.R \
      --counts_file="~{counts_matrix}" \
      --metadata_file="~{sample_metadata}" \
      --condition_column="~{condition_column}" \
      --reference_level="~{reference_level}" \
      --contrast="~{contrast}" \
      --fdr_threshold=~{fdr_threshold} \
      --lfc_threshold=~{lfc_threshold} \
      --min_count=~{min_count} \
      --min_total_count=~{min_total_count} \
      --output_prefix="edger_results"
  >>>

  output {
    File edger_results = "edger_results_all_genes.csv"
    File edger_significant = "edger_results_significant.csv"
    File edger_normalized_counts = "edger_results_normalized_cpm.csv"
    File edger_md_plot = "edger_results_md_plot.pdf"
    File edger_volcano_plot = "edger_results_volcano.pdf"
    File edger_heatmap = "edger_results_heatmap.pdf"
    File edger_bcv_plot = "edger_results_bcv_plot.pdf"
  }

  runtime {
    docker: docker_image
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}
