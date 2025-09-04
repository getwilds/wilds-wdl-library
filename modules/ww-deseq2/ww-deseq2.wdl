## WILDS WDL Module: ww-deseq2
## Description: Module for differential expression analysis using DESeq2
## Author: WILDS Team
## Contact: wilds@fredhutch.org

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-deseq2/modules/ww-testdata/ww-testdata.wdl" as testdata

workflow deseq2_example {
  meta {
    description: "Example workflow for ww-deseq2 module tasks with optional test data integration"
    outputs: {
        combined_counts_matrix: "Combined count matrix from STAR gene counts or generated test data",
        final_metadata: "Sample metadata file with experimental conditions",
        deseq2_all_results: "Complete DESeq2 differential expression results",
        deseq2_significant_results: "Filtered significant DESeq2 results",
        deseq2_normalized_counts: "DESeq2 normalized count values",
        deseq2_pca_plot: "PCA plot of samples",
        deseq2_volcano_plot: "Volcano plot of differential expression",
        deseq2_heatmap: "Heatmap of differentially expressed genes",
        validation_report: "Validation report summarizing analysis outputs"
    }
  }

  parameter_meta {
    gene_count_files: "Optional array of STAR gene count files. If not provided, test data will be generated"
    sample_names: "Optional array of sample names. If not provided, defaults will be used with test data"
    sample_conditions: "Optional array of experimental conditions. If not provided, defaults will be used with test data"
    counts_matrix: "Optional pre-combined count matrix. If provided, gene_count_files will be ignored"
    sample_metadata: "Optional sample metadata file. If provided, sample metadata generation will be skipped"
    reference_level: "Reference level for DESeq2 contrast (typically control condition)"
    contrast: "DESeq2 contrast string in format 'condition,treatment,control'"
    condition_column: "Column name in metadata containing experimental conditions"
  }

  input {
    # Optional user-provided inputs - if not provided, test data will be used
    Array[File]? gene_count_files
    Array[String]? sample_names
    Array[String]? sample_conditions
    File? counts_matrix
    File? sample_metadata
    # DESeq2 analysis parameters
    String condition_column = "condition"
    String reference_level = ""
    String contrast = ""
  }

  # Generate test data if user inputs are not provided
  if (!defined(counts_matrix) && !defined(gene_count_files)) {
    call testdata.generate_pasilla_counts as generate_test_data { input:
        n_samples = 7,
        n_genes = 10000,
        condition_name = condition_column
    }
  }

  # Determine which count files and metadata to use for combination
  Array[File] files_to_combine = select_first([gene_count_files, generate_test_data.individual_count_files])
  Array[String] names_to_use = select_first([sample_names, generate_test_data.sample_names])
  Array[String] conditions_to_use = select_first([sample_conditions, generate_test_data.sample_conditions])

  # If we have individual count files (either user-provided or generated), combine them
  # Skip this step only if user provided a pre-combined matrix
  if (!defined(counts_matrix)) {
    call combine_count_matrices { input:
        gene_count_files = files_to_combine,
        sample_names = names_to_use,
        sample_conditions = conditions_to_use
    }
  }

  # Determine which count matrix and metadata to use
  File final_counts_matrix = select_first([
    counts_matrix,  # User provided pre-combined matrix
    combine_count_matrices.counts_matrix  # Combined from individual files
  ])

  File final_sample_metadata = select_first([
    sample_metadata,  # User provided metadata
    combine_count_matrices.sample_metadata  # Generated during combination
  ])

  call run_deseq2 { input:
      counts_matrix = final_counts_matrix,
      sample_metadata = final_sample_metadata,
      condition_column = condition_column,
      reference_level = reference_level,
      contrast = contrast
  }

  call validate_deseq2_outputs { input:
      deseq2_results = run_deseq2.deseq2_results,
      deseq2_significant = run_deseq2.deseq2_significant,
      normalized_counts = run_deseq2.deseq2_normalized_counts,
      expected_samples = length(names_to_use)
  }

  output {
    File combined_counts_matrix = final_counts_matrix
    File final_metadata = final_sample_metadata
    File deseq2_all_results = run_deseq2.deseq2_results
    File deseq2_significant_results = run_deseq2.deseq2_significant
    File deseq2_normalized_counts = run_deseq2.deseq2_normalized_counts
    File deseq2_pca_plot = run_deseq2.deseq2_pca_plot
    File deseq2_volcano_plot = run_deseq2.deseq2_volcano_plot
    File deseq2_heatmap = run_deseq2.deseq2_heatmap
    File validation_report = validate_deseq2_outputs.validation_report
  }
}

task combine_count_matrices {
  meta {
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

task validate_deseq2_outputs {
  meta {
    description: "Validate DESeq2 analysis outputs for correctness and completeness"
    outputs: {
        validation_report: "Report summarizing validation results and any issues found"
    }
  }

  parameter_meta {
    deseq2_results: "DESeq2 results file to validate"
    deseq2_significant: "Significant results file to validate"
    normalized_counts: "Normalized counts file to validate"
    expected_samples: "Expected number of samples in the analysis"
    expected_genes_min: "Minimum expected number of genes in results"
  }

  input {
    File deseq2_results
    File deseq2_significant
    File normalized_counts
    Int expected_samples
    Int expected_genes_min = 1000
  }

  command <<<
    set -eo pipefail

    python3 << 'EOF'
import pandas as pd
import sys

validation_issues = []

# Validate main results file
try:
    results_df = pd.read_csv("~{deseq2_results}")
    print(f"Results file loaded: {len(results_df)} genes")
    
    required_cols = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    missing_cols = [col for col in required_cols if col not in results_df.columns]
    if missing_cols:
        validation_issues.append(f"Missing required columns in results: {missing_cols}")
    else:
        print("All required columns present in results")
        
    if len(results_df) < ~{expected_genes_min}:
        validation_issues.append(f"Too few genes in results: {len(results_df)} < ~{expected_genes_min}")
        
except Exception as e:
    validation_issues.append(f"Failed to validate results file: {str(e)}")

# Validate significant results  
try:
    sig_df = pd.read_csv("~{deseq2_significant}")
    print(f"Significant results file loaded: {len(sig_df)} genes")
    
    if len(sig_df) > len(results_df):
        validation_issues.append("More significant genes than total genes")
        
except Exception as e:
    validation_issues.append(f"Failed to validate significant results file: {str(e)}")

# Validate normalized counts
try:
    counts_df = pd.read_csv("~{normalized_counts}", index_col=0)
    print(f"Normalized counts loaded: {counts_df.shape[0]} genes, {counts_df.shape[1]} samples")
    
    if counts_df.shape[1] != ~{expected_samples}:
        validation_issues.append(f"Expected ~{expected_samples} samples, got {counts_df.shape[1]}")
        
    if (counts_df < 0).any().any():
        validation_issues.append("Negative values found in normalized counts")
        
except Exception as e:
    validation_issues.append(f"Failed to validate normalized counts: {str(e)}")

# Write validation report
with open("validation_report.txt", "w") as f:
    f.write("DESeq2 Output Validation Report\n")
    f.write("=" * 35 + "\n\n")
    
    if validation_issues:
        f.write("VALIDATION FAILED\n")
        f.write("Issues found:\n")
        for issue in validation_issues:
            f.write(f"- {issue}\n")
        sys.exit(1)
    else:
        f.write("VALIDATION PASSED\n")
        f.write("All outputs validated successfully.\n")

EOF

    echo "Validation completed successfully"
  >>>

  output {
    File validation_report = "validation_report.txt"
  }

  runtime {
    docker: "python:3.9-slim"
    memory: "2 GB"
    cpu: 1
  }
}
