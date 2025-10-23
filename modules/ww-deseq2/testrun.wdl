version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deseq2/ww-deseq2.wdl" as ww_deseq2

workflow deseq2_example {
  # Generate test data
  call ww_testdata.generate_pasilla_counts as generate_test_data { input:
      n_samples = 7,
      n_genes = 10000,
      condition_name = "condition"
  }

  # Combine the test count files
  call ww_deseq2.combine_count_matrices { input:
      gene_count_files = generate_test_data.individual_count_files,
      sample_names = generate_test_data.sample_names,
      sample_conditions = generate_test_data.sample_conditions
  }

  call ww_deseq2.run_deseq2 { input:
      counts_matrix = combine_count_matrices.counts_matrix,
      sample_metadata = combine_count_matrices.sample_metadata,
      condition_column = "condition",
      reference_level = "",
      contrast = ""
  }

  call validate_outputs { input:
      deseq2_results = run_deseq2.deseq2_results,
      deseq2_significant = run_deseq2.deseq2_significant,
      normalized_counts = run_deseq2.deseq2_normalized_counts,
      expected_samples = length(generate_test_data.sample_names)
  }

  output {
    File combined_counts_matrix = combine_count_matrices.counts_matrix
    File final_metadata = combine_count_matrices.sample_metadata
    File deseq2_all_results = run_deseq2.deseq2_results
    File deseq2_significant_results = run_deseq2.deseq2_significant
    File deseq2_normalized_counts = run_deseq2.deseq2_normalized_counts
    File deseq2_pca_plot = run_deseq2.deseq2_pca_plot
    File deseq2_volcano_plot = run_deseq2.deseq2_volcano_plot
    File deseq2_heatmap = run_deseq2.deseq2_heatmap
    File validation_report = validate_outputs.validation_report
  }
}

task validate_outputs {
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

    pip3 install pandas
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
