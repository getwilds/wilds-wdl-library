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
    topic: "transcriptomics,gene_expression"
    species: "human,eukaryote,prokaryote,virus"
    operation: "aggregation"
    input_sample_required: "gene_count_files:gene_expression_matrix:tsv"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "counts_matrix:gene_expression_matrix:matrix,sample_metadata:report:textual_format"
    output_reference: "none"
  }

  parameter_meta {
    gene_count_files: "Array of STAR gene count files (ReadsPerGene.out.tab) from each sample"
    sample_names: "Array of sample names corresponding to the gene_count_files"
    sample_conditions: "Array of experimental conditions corresponding to each sample"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
    count_column: "Column number in STAR count files to extract (2=unstranded, 3=forward strand, 4=reverse strand)"
    docker_image: "Docker image to use for this task"
  }

  input {
    Array[File] gene_count_files
    Array[String] sample_names
    Array[String] sample_conditions
    Int memory_gb = 4
    Int cpu_cores = 1
    String docker_image = "getwilds/python-utils:0.1.0"
    Int count_column = 2
        # Column to extract from ReadsPerGene.out.tab files:
        # 2 = unstranded counts
        # 3 = stranded counts, first read forward
        # 4 = stranded counts, first read reverse
  }

  command <<<
    set -eo pipefail

    curl -so combine_star_counts.py \
      "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deseq2/combine_star_counts.py"

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
    docker: docker_image
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task run_deseq2 {
  meta {
    author: [
        {
            name: "Taylor Firman",
            email: "tfirman@fredhutch.org"
        },
        {
            name: "Elaine Glenny",
            email: "eglenny@fredhutch.org"
        }
    ]
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
    topic: "transcriptomics,gene_expression"
    species: "human,eukaryote,prokaryote,virus"
    operation: "statistical_calculation"
    input_sample_required: "counts_matrix:gene_expression_matrix:matrix,sample_metadata:report:textual_format"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "deseq2_results:report:csv,deseq2_significant:report:csv,deseq2_normalized_counts:report:csv,deseq2_pca_plot:plot:pdf,deseq2_volcano_plot:plot:pdf,deseq2_heatmap:plot:pdf"
    output_reference: "none"
  }

  parameter_meta {
    counts_matrix: "Combined matrix of gene-level counts from all samples"
    sample_metadata: "Metadata file containing sample information including conditions"
    condition_column: "Column name in the metadata file that contains the experimental condition"
    reference_level: "Reference level for the contrast (typically the control condition)"
    contrast: "DESeq2 contrast string in the format 'condition,treatment,control'"
    min_counts: "Minimum number of counts a gene must have to pass filtering"
    min_samples: "A gene must meet the `min_counts` threshold in this many samples to be kept. 0 = Gene is kept if its count across all samples meets the `min_count` threshold (default: 0)"
    shrinkage_method: "LFC shrinkage method: apeglm, ashr, or normal (empty = no shrinkage)"
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
    Int min_counts = 10
    Int min_samples = 0
    String shrinkage_method = ""
    Int memory_gb = 8
    Int cpu_cores = 2
    String docker_image = "getwilds/deseq2:1.40.2"
  }

  command <<<
    set -eo pipefail

    curl -so deseq2_analysis.R \
      "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deseq2/deseq2_analysis.R"

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
    docker: docker_image
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task run_deseq2_biocjobs {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Run the canonical DESeq2 differential expression workflow on a raw count matrix using the biocjobs package"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deseq2/ww-deseq2.wdl"
    outputs: {
        results: "Per-gene DESeq2 results table with log2 fold changes and adjusted p-values",
        normalized_counts: "Median-of-ratios normalized count values for all samples",
        ma_plot: "MA plot showing log fold change vs. mean expression"
    }
    topic: "transcriptomics,gene_expression"
    species: "human,eukaryote,prokaryote,virus"
    operation: "statistical_calculation"
    input_sample_required: "counts:gene_expression_matrix:tsv,coldata:report:tsv"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "results:report:tsv,normalized_counts:report:tsv,ma_plot:plot:pdf"
    output_reference: "none"
    biocjobs_package: "DESeq2"
    biocjobs_job: "deseq2-differential-expression"
    biocjobs_job_version: "1.0.0"
  }
  
  parameter_meta {
    counts: "Raw count matrix. Tab-separated matrix of raw (un-normalized) integer read counts: first column gene identifiers, one further column per sample. Do not supply TPM/FPKM or otherwise normalized values.. Format: tsv."
    coldata: "Sample table. Tab-separated sample annotation: first column sample identifiers matching the count matrix column names, one further column per covariate (e.g. condition, batch). Every variable used in the design formula must be a column of this table.. Format: tsv."
    design: "Design formula. R formula over sample table columns. Control for covariates by putting them first; the tested variable comes last (e.g. \"~ batch + condition\")."
    contrast_factor: "Factor to test. Sample table column containing the tested condition."
    contrast_numerator: "Treatment level (numerator). Level of the tested factor whose effect is reported."
    contrast_denominator: "Reference level (denominator). Baseline level the treatment is compared against."
    alpha: "FDR threshold (alpha). Significance cutoff used to optimize independent filtering. Set it to the adjusted p-value threshold you will actually use.. Range: 0 to 1."
    shrinkage: "Log2 fold change shrinkage. Shrinkage estimator applied to the reported log2 fold changes; \"apeglm\" is the recommended default, \"none\" reports maximum likelihood estimates.. One of: apeglm, ashr, normal, none."
    prefilter: "Pre-filter low-count genes. Remove genes with insufficient counts before testing."
    prefilter_min_count: "Pre-filter minimum count. Range: 1 to Inf."
    prefilter_min_samples: "Pre-filter minimum samples. Keep genes with at least the minimum count in at least this many samples. 0 (default) means automatic: the size of the smallest group of the tested factor.. Range: 0 to Inf."
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
    docker_image: "Docker image to use for this task"
  }

  input {
    File counts
    File coldata
    String design = "~ condition"
    String contrast_factor
    String contrast_numerator
    String contrast_denominator
    String alpha = "0.1"
    String shrinkage = "apeglm"
    Boolean prefilter = true
    Int prefilter_min_count = 10
    Int prefilter_min_samples = 0
    Int memory_gb = 8
    Int cpu_cores = 2
    String docker_image = "ghcr.io/almahmoud/deseq2:devel"
  }

  command <<<
    set -euo pipefail
    Rscript -e 'BiocJobs::execJob("DESeq2", "deseq2-differential-expression")' \
      --counts '~{sub(counts, "'", "'\\''")}' \
      --coldata '~{sub(coldata, "'", "'\\''")}' \
      --design '~{sub(design, "'", "'\\''")}' \
      --contrast_factor '~{sub(contrast_factor, "'", "'\\''")}' \
      --contrast_numerator '~{sub(contrast_numerator, "'", "'\\''")}' \
      --contrast_denominator '~{sub(contrast_denominator, "'", "'\\''")}' \
      --alpha '~{sub(alpha, "'", "'\\''")}' \
      --shrinkage '~{sub(shrinkage, "'", "'\\''")}' \
      --prefilter ~{prefilter} \
      --prefilter_min_count ~{prefilter_min_count} \
      --prefilter_min_samples ~{prefilter_min_samples} \
      --results 'results.tsv' \
      --normalized_counts 'normalized_counts.tsv' \
      --ma_plot 'ma_plot.pdf'
  >>>

  output {
    File results = "results.tsv"
    File normalized_counts = "normalized_counts.tsv"
    File ma_plot = "ma_plot.pdf"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task compile_deseq2_results {
  meta {
    author: [
        {
            name: "Taylor Firman",
            email: "tfirman@fredhutch.org"
        },
        {
            name: "Elaine Glenny",
            email: "eglenny@fredhutch.org"
        }
    ]
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
    docker_image: "Docker image to use for this task"
  }

  input {
    File deseq2_results
    File normalized_counts
    File gtf_file
    String output_name = "deseq2_compiled_results.csv"
    Int memory_gb = 8
    Int cpu_cores = 2
    String docker_image = "getwilds/python-utils:0.1.0"
  }

  command <<<
    set -eo pipefail

    curl -so compile_deseq2_results.py \
      "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deseq2/compile_deseq2_results.py"

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
    docker: docker_image
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

