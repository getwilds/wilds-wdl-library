version 1.2

# HPC test workflow for ww-colabfold. Mirrors testrun.wdl coverage but uses
# realistic protein input and full GPU + AMBER settings — intended to validate
# end-to-end ColabFold behavior on Fred Hutch HPC infrastructure.

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/split-cicd-hpc-testruns/modules/ww-colabfold/ww-colabfold.wdl" as ww_colabfold
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/split-cicd-hpc-testruns/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####

workflow colabfold_example {
  # Use a realistic-size protein (human ubiquitin, 76 residues) instead of a tiny peptide
  call ww_testdata.create_realistic_protein_fasta { }

  # Download AlphaFold2 model weights
  call ww_colabfold.download_weights { input:
      cpu_cores = 2,
      memory_gb = 8
  }

  # Production-style ColabFold prediction on GPU with AMBER relaxation
  call ww_colabfold.colabfold_predict { input:
      fasta_file = create_realistic_protein_fasta.test_fasta,
      weights_tarball = download_weights.weights_tarball,
      output_prefix = "test_ubiquitin",
      num_recycle = 3,
      num_models = 5,
      num_seeds = 1,
      msa_mode = "mmseqs2_uniref_env",
      use_amber = true,
      stop_at_score = 85,
      use_gpu_relax = true,
      model_type = "auto",
      gpu_enabled = true,
      cpu_cores = 8,
      memory_gb = 48
  }

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

    echo "=== ColabFold HPC Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    tar -xzf ~{results_tarball}

    validation_passed=true

    if [ -d "results" ]; then
      file_count=$(find results/ -type f | wc -l)
      echo "Results directory: found ${file_count} files" >> validation_report.txt
    else
      echo "Results directory: MISSING" >> validation_report.txt
      validation_passed=false
    fi

    pdb_count=$(find results/ -name "*.pdb" 2>/dev/null | wc -l)
    echo "PDB structure files: ${pdb_count}" >> validation_report.txt
    if [ "${pdb_count}" -eq 0 ]; then
      echo "WARNING: No PDB files found" >> validation_report.txt
      validation_passed=false
    fi

    echo "" >> validation_report.txt
    echo "--- Output Files ---" >> validation_report.txt
    find results/ -type f -exec ls -lh {} \; >> validation_report.txt 2>/dev/null || true

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
