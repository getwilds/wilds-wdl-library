version 1.2

# HPC test workflow for ww-starling. Mirrors testrun.wdl coverage with GPU
# enabled on the inference tasks and the module's default 400 conformations
# — intended to validate end-to-end behavior on Fred Hutch HPC infrastructure.

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/wdl-v1.2-support/modules/ww-starling/ww-starling.wdl" as ww_starling
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/hpc-testruns/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####

workflow starling_example {
  # p53 N-terminal transactivation domain (residues 1-39) — a well-studied IDP
  String test_sequence = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQA"

  # Single-sequence ensemble generation on GPU
  call ww_starling.generate_ensemble { input:
    sequence = test_sequence,
    sample_name = "test_p53_ntad",
    num_conformations = 400,
    gpu_enabled = true,
    cpu_cores = 4,
    memory_gb = 8
  }

  # Ensemble metadata query (CPU-only task)
  call ww_starling.ensemble_info { input:
    starling_file = generate_ensemble.starling_file,
    sample_name = "test_p53_ntad"
  }

  # Multi-sequence batch ensemble generation on GPU
  call ww_testdata.create_test_idp_fasta { }

  call ww_starling.generate_ensemble_batch { input:
    fasta_file = create_test_idp_fasta.test_fasta,
    num_conformations = 400,
    gpu_enabled = true,
    cpu_cores = 4,
    memory_gb = 8
  }

  # FASTA splitting (CPU-only task)
  call ww_starling.split_fasta { input:
    fasta_file = create_test_idp_fasta.test_fasta,
    sequences_per_batch = 2
  }

  call validate_outputs { input:
      starling_file = generate_ensemble.starling_file,
      pdb_file = generate_ensemble.pdb_file,
      xtc_file = generate_ensemble.xtc_file,
      info_file = ensemble_info.info_file,
      batch_starling_files = generate_ensemble_batch.starling_files,
      split_batch_files = split_fasta.batch_files
  }

  output {
    File starling_file = generate_ensemble.starling_file
    File pdb_file = generate_ensemble.pdb_file
    File xtc_file = generate_ensemble.xtc_file
    File info_file = ensemble_info.info_file
    Array[File] batch_starling_files = generate_ensemble_batch.starling_files
    Array[File] batch_pdb_files = generate_ensemble_batch.pdb_files
    Array[File] batch_xtc_files = generate_ensemble_batch.xtc_files
    Array[File] split_batch_files = split_fasta.batch_files
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that STARLING produced expected ensemble, conversion, info, and split outputs"
    outputs: {
        report: "Validation report summarizing file checks"
    }
  }

  parameter_meta {
    starling_file: "STARLING ensemble file from single-sequence generate_ensemble"
    pdb_file: "PDB topology file from single-sequence generate_ensemble"
    xtc_file: "XTC trajectory file from single-sequence generate_ensemble"
    info_file: "Metadata file from ensemble_info"
    batch_starling_files: "Array of STARLING ensemble files from generate_ensemble_batch"
    split_batch_files: "Array of FASTA batch files from split_fasta"
  }

  input {
    File starling_file
    File pdb_file
    File xtc_file
    File info_file
    Array[File] batch_starling_files
    Array[File] split_batch_files
  }

  command <<<
    set -eo pipefail

    echo "=== STARLING HPC Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    for f in ~{starling_file} ~{pdb_file} ~{xtc_file} ~{info_file}; do
      if [ -s "${f}" ]; then
        echo "OK: ${f} ($(wc -c < "${f}") bytes)" >> validation_report.txt
      else
        echo "MISSING/EMPTY: ${f}" >> validation_report.txt
        validation_passed=false
      fi
    done

    batch_count=$(echo "~{sep=' ' batch_starling_files}" | wc -w)
    echo "Batch ensemble files produced: ${batch_count}" >> validation_report.txt
    if [ "${batch_count}" -eq 0 ]; then
      echo "WARNING: no batch ensemble files produced" >> validation_report.txt
      validation_passed=false
    fi

    split_count=$(echo "~{sep=' ' split_batch_files}" | wc -w)
    echo "FASTA batch files produced: ${split_count}" >> validation_report.txt
    if [ "${split_count}" -eq 0 ]; then
      echo "WARNING: no FASTA batch files produced" >> validation_report.txt
      validation_passed=false
    fi

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
