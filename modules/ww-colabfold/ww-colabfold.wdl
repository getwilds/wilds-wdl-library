## WILDS WDL for predicting protein structures using ColabFold.
## ColabFold combines AlphaFold2 with MMseqs2 for fast, accessible protein structure prediction.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

#### TASK DEFINITIONS ####

task download_weights {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Download AlphaFold2 model weights for use with ColabFold predictions"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-colabfold/ww-colabfold.wdl"
    outputs: {
        weights_tarball: "Compressed tarball containing AlphaFold2 model weights (~15-20 GB)"
    }
  }

  parameter_meta {
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    Int cpu_cores = 2
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Force CPU-only mode for weight download (no GPU needed)
    export JAX_PLATFORMS=cpu

    # Download AlphaFold2 model weights to the default cache location
    python -m colabfold.download

    # Package the weights cache into a tarball for reuse across predictions
    tar -czf colabfold_weights.tar.gz -C /cache .
  >>>

  output {
    File weights_tarball = "colabfold_weights.tar.gz"
  }

  runtime {
    docker: "ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task colabfold_predict {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Predict protein structures from amino acid sequences using ColabFold (AlphaFold2 + MMseqs2)"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-colabfold/ww-colabfold.wdl"
    outputs: {
        results_tarball: "Compressed tarball containing all ColabFold outputs including PDB structure files, PAE plots, pLDDT plots, coverage plots, and prediction metrics"
    }
  }

  parameter_meta {
    fasta_file: "Input FASTA file containing one or more protein sequences to predict"
    weights_tarball: "Compressed tarball of AlphaFold2 model weights from the download_weights task"
    output_prefix: "Prefix for naming the output tarball"
    num_recycle: "Number of prediction recycles (higher values may improve accuracy)"
    num_models: "Number of models to generate per sequence"
    num_seeds: "Number of random seeds per model"
    msa_mode: "MSA generation mode: mmseqs2_uniref_env (default, uses MMseqs2 server), mmseqs2_uniref, or single_sequence"
    use_amber: "Enable AMBER force field for structure relaxation/energy minimization"
    stop_at_score: "Stop prediction early when this pLDDT confidence score is reached"
    use_gpu_relax: "Run AMBER relaxation on GPU instead of CPU"
    model_type: "Model type: auto (recommended), alphafold2, or alphafold2_multimer_v3"
    extra_flags: "Additional command-line flags to pass to colabfold_batch"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
    gpu_enabled: "Enable GPU for prediction (controls JAX CPU/GPU mode and PROOF GPU allocation)"
  }

  input {
    File fasta_file
    File weights_tarball
    String output_prefix
    Int num_recycle = 3
    Int num_models = 5
    Int num_seeds = 1
    String msa_mode = "mmseqs2_uniref_env"
    Boolean use_amber = true
    Int stop_at_score = 85
    Boolean use_gpu_relax = true
    String model_type = "auto"
    String extra_flags = ""
    Int cpu_cores = 8
    Int memory_gb = 48
    Boolean gpu_enabled = true
  }

  command <<<
    set -eo pipefail

    # Force CPU-only execution if GPU is not available
    if [ "~{gpu_enabled}" = "false" ]; then
      export JAX_PLATFORMS=cpu
    fi

    # Extract pre-downloaded model weights into the cache directory
    mkdir -p /cache
    tar -xzf ~{weights_tarball} -C /cache

    # Run colabfold_batch prediction
    colabfold_batch \
      ~{fasta_file} \
      results/ \
      --num-recycle ~{num_recycle} \
      --num-models ~{num_models} \
      --num-seeds ~{num_seeds} \
      --msa-mode ~{msa_mode} \
      ~{if use_amber then "--amber" else ""} \
      ~{if use_gpu_relax then "--use-gpu-relax" else ""} \
      --stop-at-score ~{stop_at_score} \
      --model-type ~{model_type} \
      ~{extra_flags}

    # Package all results into a compressed tarball
    tar -czf "~{output_prefix}_colabfold_results.tar.gz" results/
  >>>

  output {
    File results_tarball = "~{output_prefix}_colabfold_results.tar.gz"
  }

  runtime {
    docker: "ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
    gpus: if gpu_enabled then "1" else "0"
  }
}
