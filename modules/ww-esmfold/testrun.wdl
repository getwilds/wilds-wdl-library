version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-esmfold/ww-esmfold.wdl" as ww_esmfold
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####
# Tests ESMFold structure prediction with a tiny protein sequence on CPU.
# Uses the Trp-cage miniprotein (20 residues) for fast execution.
# Note: Production use requires GPU; this test demonstrates the workflow on CPU.

workflow esmfold_example {
  # Create a minimal test FASTA with a short peptide sequence
  call ww_testdata.create_test_protein_fasta { }

  # Run ESMFold prediction with minimal settings for CI testing
  call ww_esmfold.esmfold_predict { input:
      fasta_file = create_test_protein_fasta.test_fasta,
      output_prefix = "test_protein",
      num_recycles = 1,
      chunk_size = 0,
      cpu_only = true,
      cpu_offload = false,
      max_tokens_per_batch = 0,
      cpu_cores = 2,
      memory_gb = 24,
      gpu_enabled = false
  }

  # Validate that prediction outputs were generated
  call validate_outputs { input:
      results_tarball = esmfold_predict.pdb_output
  }

  output {
    File esmfold_results = esmfold_predict.pdb_output
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that ESMFold prediction generated expected output files"
    outputs: {
        report: "Validation report summarizing file checks and prediction results"
    }
  }

  parameter_meta {
    results_tarball: "Compressed tarball of ESMFold prediction results to validate"
  }

  input {
    File results_tarball
  }

  command <<<
    set -eo pipefail

    echo "=== ESMFold Prediction Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    # Extract results
    tar -xzf ~{results_tarball}

    validation_passed=true

    # Check that the output directory exists and has files
    if [ -d "pdb_output" ]; then
      file_count=$(find pdb_output/ -type f | wc -l)
      echo "Output directory: found ${file_count} files" >> validation_report.txt
    else
      echo "Output directory: MISSING" >> validation_report.txt
      validation_passed=false
    fi

    # Check for PDB structure files
    pdb_count=$(find pdb_output/ -name "*.pdb" 2>/dev/null | wc -l)
    echo "PDB structure files: ${pdb_count}" >> validation_report.txt
    if [ "${pdb_count}" -eq 0 ]; then
      echo "WARNING: No PDB files found" >> validation_report.txt
      validation_passed=false
    fi

    # List all output files
    echo "" >> validation_report.txt
    echo "--- Output Files ---" >> validation_report.txt
    find pdb_output/ -type f -exec ls -lh {} \; >> validation_report.txt 2>/dev/null || true

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
