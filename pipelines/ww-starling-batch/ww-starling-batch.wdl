version 1.0

import "../../modules/ww-starling/ww-starling.wdl" as starling_tasks

workflow starling_batch {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Batch ensemble generation pipeline that splits a large protein FASTA into chunks and processes them in parallel using STARLING"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-starling-batch/ww-starling-batch.wdl"
    outputs: {
        starling_files: "All STARLING ensemble files across all batches",
        pdb_files: "All PDB topology files across all batches",
        xtc_files: "All XTC trajectory files across all batches"
    }
  }

  parameter_meta {
    fasta_file: "Input FASTA file containing protein sequences for ensemble generation"
    sequences_per_batch: "Number of sequences to include in each batch for parallel processing"
    num_conformations: "Number of conformations to generate per sequence"
    gpu_enabled: "Enable GPU for STARLING inference in each batch task"
    cpu_cores: "Number of CPU cores allocated per batch task"
    memory_gb: "Memory allocated per batch task in GB"
  }

  input {
    File fasta_file
    Int sequences_per_batch = 10
    Int num_conformations = 400
    Boolean gpu_enabled = true
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  # Split the input FASTA into batches
  call starling_tasks.split_fasta { input:
    fasta_file = fasta_file,
    sequences_per_batch = sequences_per_batch
  }

  # Scatter ensemble generation across batches
  scatter (batch_fasta in split_fasta.batch_files) {
    call starling_tasks.generate_ensemble_batch { input:
      fasta_file = batch_fasta,
      num_conformations = num_conformations,
      gpu_enabled = gpu_enabled,
      cpu_cores = cpu_cores,
      memory_gb = memory_gb
    }
  }

  output {
    Array[Array[File]] starling_files = generate_ensemble_batch.starling_files
    Array[Array[File]] pdb_files = generate_ensemble_batch.pdb_files
    Array[Array[File]] xtc_files = generate_ensemble_batch.xtc_files
  }
}
