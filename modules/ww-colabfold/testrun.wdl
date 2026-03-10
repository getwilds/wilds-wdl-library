version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-colabfold/ww-colabfold.wdl" as ww_colabfold
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####
# Tests ColabFold prediction with a tiny protein sequence on CPU.
# Uses single_sequence MSA mode and minimal model settings for fast execution.
# Note: Production use requires GPU; this test demonstrates the workflow on CPU.

workflow colabfold_example {
  # Create a minimal test FASTA with a short peptide sequence
  call ww_testdata.create_test_protein_fasta { }

  # Download AlphaFold2 model weights (download once, reuse across predictions)
  call ww_colabfold.download_weights { input:
      cpu_cores = 2,
      memory_gb = 8
  }

  # Run ColabFold prediction with minimal settings for CI testing
  call ww_colabfold.colabfold_predict { input:
      fasta_file = create_test_protein_fasta.test_fasta,
      weights_tarball = download_weights.weights_tarball,
      output_prefix = "test_protein",
      num_recycle = 1,
      num_models = 1,
      num_seeds = 1,
      msa_mode = "single_sequence",
      use_amber = false,
      stop_at_score = 50,
      use_gpu_relax = false,
      model_type = "auto",
      cpu_cores = 2,
      memory_gb = 8,
      gpu_enabled = false
  }

  # Validate that prediction outputs were generated
  call validate_outputs { input:
      results_tarball = colabfold_predict.results_tarball
  }

  output {
    File colabfold_weights = download_weights.weights_tarball
    File colabfold_results = colabfold_predict.results_tarball
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that ColabFold prediction generated expected output files"
    outputs: {
        report: "Validation report summarizing file checks and prediction results"
    }
  }

  parameter_meta {
    results_tarball: "Compressed tarball of ColabFold prediction results to validate"
  }

  input {
    File results_tarball
  }

  command <<<
    set -eo pipefail

    echo "=== ColabFold Prediction Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    # Extract results
    tar -xzf ~{results_tarball}

    validation_passed=true

    # Check that the results directory exists and has files
    if [ -d "results" ]; then
      file_count=$(find results/ -type f | wc -l)
      echo "Results directory: found ${file_count} files" >> validation_report.txt
    else
      echo "Results directory: MISSING" >> validation_report.txt
      validation_passed=false
    fi

    # Check for PDB structure files
    pdb_count=$(find results/ -name "*.pdb" 2>/dev/null | wc -l)
    echo "PDB structure files: ${pdb_count}" >> validation_report.txt
    if [ "${pdb_count}" -eq 0 ]; then
      echo "WARNING: No PDB files found" >> validation_report.txt
      validation_passed=false
    fi

    # List all output files
    echo "" >> validation_report.txt
    echo "--- Output Files ---" >> validation_report.txt
    find results/ -type f -exec ls -lh {} \; >> validation_report.txt 2>/dev/null || true

    # Summary
    echo "" >> validation_report.txt
    echo "=== Validation Summary ===" >> validation_report.txt
    if [ "${validation_passed}" = "true" ]; then
      echo "Overall Status: PASSED" >> validation_report.txt
    else
      echo "Overall Status: FAILED" >> validation_report.txt
      exit 1
    fi

    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: 1
    memory: "2 GB"
  }
}
