## WILDS WDL Module: ww-deseq2
## Description: Module for differential expression analysis using DESeq2
## Author: Taylor Firman
## Contact: tfirman@fredhutch.org

version 1.0

task combine_count_matrices {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Combine STAR gene count files from multiple samples into a single count matrix"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deseq2/ww-deseq2.wdl"
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

    pip install pandas==2.2.3

    curl -so combine_star_counts.py \
      "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/rnaseq-feedback/modules/ww-deseq2/combine_star_counts.py"

    python combine_star_counts.py \
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
    docker: "python:3.12"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task run_deseq2 {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Perform differential expression analysis using DESeq2"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deseq2/ww-deseq2.wdl"
    outputs: {
        deseq2_results: "Complete results of DESeq2 differential expression analysis with statistics",
        deseq2_significant: "Filtered results containing only statistically significant differentially expressed genes",
        deseq2_normalized_counts: "Normalized count values produced by DESeq2 for all samples",
        deseq2_pca_plot: "Principal Component Analysis plot showing sample clustering based on expression patterns",
        deseq2_volcano_plot: "Volcano plot showing log fold change vs. statistical significance",
        deseq2_heatmap: "Heatmap visualization of differentially expressed genes across samples",
        deseq2_ma_plot: "MA plot showing log fold change vs. mean expression (unshrunken)",
        deseq2_ma_plot_shrunk: "MA plot with shrunken log fold changes (empty if shrinkage not applied)",
        deseq2_results_shrunk: "DESeq2 results with shrunken log fold changes (empty if shrinkage not applied)"
    }
  }

  parameter_meta {
    counts_matrix: "Combined matrix of gene-level counts from all samples"
    sample_metadata: "Metadata file containing sample information including conditions"
    condition_column: "Column name in the metadata file that contains the experimental condition"
    reference_level: "Reference level for the contrast (typically the control condition)"
    contrast: "DESeq2 contrast string in the format 'condition,treatment,control'"
    min_counts: "Minimum number of counts a gene must have to pass filtering"
    min_samples: "Minimum number of samples that must meet min_counts threshold (0 = use total counts instead)"
    shrinkage_method: "LFC shrinkage method: apeglm, ashr, or normal (empty = no shrinkage)"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
  }

  input {
    File counts_matrix
    File sample_metadata
    String condition_column = "condition"
    String reference_level = ""
    String contrast = ""
    Int min_counts = 10
    Int min_samples = 0
    String shrinkage_method = ""
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail

    curl -so deseq2_analysis.R \
      "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/rnaseq-feedback/modules/ww-deseq2/deseq2_analysis.R"

    Rscript deseq2_analysis.R \
      --counts_file="~{counts_matrix}" \
      --metadata_file="~{sample_metadata}" \
      --condition_column="~{condition_column}" \
      --reference_level="~{reference_level}" \
      --contrast="~{contrast}" \
      --min_counts=~{min_counts} \
      --min_samples=~{min_samples} \
      --shrinkage_method="~{shrinkage_method}" \
      --output_prefix="deseq2_results"
  >>>

  output {
    File deseq2_results = "deseq2_results_all_genes.csv"
    File deseq2_significant = "deseq2_results_significant.csv"
    File deseq2_normalized_counts = "deseq2_results_normalized_counts.csv"
    File deseq2_pca_plot = "deseq2_results_pca.pdf"
    File deseq2_volcano_plot = "deseq2_results_volcano.pdf"
    File deseq2_heatmap = "deseq2_results_heatmap.pdf"
    File deseq2_ma_plot = "deseq2_results_ma_plot.pdf"
    File deseq2_ma_plot_shrunk = "deseq2_results_ma_plot_shrunk.pdf"
    File deseq2_results_shrunk = "deseq2_results_all_genes_shrunk.csv"
  }

  runtime {
    docker: "getwilds/deseq2:1.40.2"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task compile_deseq2_results {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Merge DESeq2 results with normalized counts and GTF gene annotations into a single comprehensive output file"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deseq2/ww-deseq2.wdl"
    outputs: {
        compiled_results: "Combined CSV with DESeq2 statistics, normalized counts, and gene annotations from the GTF"
    }
  }

  parameter_meta {
    deseq2_results: "DESeq2 all-genes results CSV (output of run_deseq2)"
    normalized_counts: "DESeq2 normalized counts CSV (output of run_deseq2)"
    gtf_file: "GTF annotation file used for alignment, for extracting gene descriptions"
    output_name: "Name for the output CSV file"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
  }

  input {
    File deseq2_results
    File normalized_counts
    File gtf_file
    String output_name = "deseq2_compiled_results.csv"
    Int memory_gb = 4
    Int cpu_cores = 1
  }

  command <<<
    set -eo pipefail

    pip install pandas==2.2.3

    curl -so compile_deseq2_results.py \
      "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/rnaseq-feedback/modules/ww-deseq2/compile_deseq2_results.py"

    python compile_deseq2_results.py \
      --results "~{deseq2_results}" \
      --counts "~{normalized_counts}" \
      --gtf "~{gtf_file}" \
      --output "~{output_name}"
  >>>

  output {
    File compiled_results = "~{output_name}"
  }

  runtime {
    docker: "python:3.12"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

