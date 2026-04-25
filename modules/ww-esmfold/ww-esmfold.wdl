## WILDS WDL for predicting protein structures using ESMFold.
## ESMFold uses the ESM-2 protein language model for fast, single-sequence
## protein structure prediction without requiring multiple sequence alignments.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

#### TASK DEFINITIONS ####

task esmfold_predict {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Predict protein 3D structures from amino acid sequences using ESMFold (ESM-2 language model)"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-esmfold/ww-esmfold.wdl"
    outputs: {
        pdb_output: "Directory of predicted PDB structure files with pLDDT confidence scores in the B-factor column"
    }
  }

  parameter_meta {
    fasta_file: "Input FASTA file containing one or more protein sequences to predict"
    output_prefix: "Prefix for naming the output tarball"
    num_recycles: "Number of recycling iterations for structure refinement (higher may improve accuracy)"
    chunk_size: "Chunk size for memory-efficient inference on long sequences (e.g., 128, 64, 32). Set to 0 to disable chunking."
    cpu_only: "Run inference on CPU instead of GPU"
    cpu_offload: "Enable CPU offloading to handle longer sequences with limited GPU memory"
    max_tokens_per_batch: "Maximum tokens per batch for GPU inference (0 for automatic batching)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
    gpu_enabled: "Enable GPU for prediction (recommended for production use)"
  }

  input {
    File fasta_file
    String output_prefix
    Int num_recycles = 4
    Int chunk_size = 128
    Boolean cpu_only = false
    Boolean cpu_offload = false
    Int max_tokens_per_batch = 0
    Int cpu_cores = 4
    Int memory_gb = 16
    Boolean gpu_enabled = true
  }

  command <<<
    set -eo pipefail

    # Redirect model/temp caches to working directory to avoid filling container root filesystem
    TORCH_HOME="$(pwd)/.cache/torch"
    HF_HOME="$(pwd)/.cache/huggingface"
    TMPDIR="$(pwd)/.cache/tmp"
    export TORCH_HOME HF_HOME TMPDIR
    mkdir -p "${TORCH_HOME}" "${HF_HOME}" "${TMPDIR}" pdb_output

    # Build the esm-fold command
    esm-fold \
      -i ~{fasta_file} \
      -o pdb_output/ \
      --num-recycles ~{num_recycles} \
      ~{if chunk_size > 0 then "--chunk-size " + chunk_size else ""} \
      ~{if cpu_only then "--cpu-only" else ""} \
      ~{if cpu_offload then "--cpu-offload" else ""} \
      ~{if max_tokens_per_batch > 0 then "--max-tokens-per-batch " + max_tokens_per_batch else ""}

    # Package all PDB outputs into a compressed tarball
    tar -czf "~{output_prefix}_esmfold_results.tar.gz" pdb_output/
  >>>

  output {
    File pdb_output = "~{output_prefix}_esmfold_results.tar.gz"
  }

  runtime {
    docker: "getwilds/esmfold:2.0.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
    gpus: if gpu_enabled then "1" else "0"
  }
}
