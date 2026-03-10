# ww-colabfold Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for predicting protein 3D structures from amino acid sequences using [ColabFold](https://github.com/sokrypton/ColabFold), which combines AlphaFold2 with the fast homology search of MMseqs2.

## Overview

This module wraps the `colabfold_batch` command-line tool for protein structure prediction. ColabFold provides a 40-60x speedup over standard AlphaFold2 by replacing the homology search step with MMseqs2, while matching AlphaFold2's prediction accuracy on CASP14 benchmarks. It supports monomeric proteins, multimers/complexes, and single-sequence predictions.

**Key benefits over using ColabFold in a Google Colab notebook:**
- Potentially higher throughput
- Files can be loaded from and saved to your own filesystem

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-colabfold.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `download_weights`

Download AlphaFold2 model weights for use with ColabFold predictions. Weights are ~15-20 GB and only need to be downloaded once, then can be reused across multiple prediction tasks.

**Inputs:**
- `cpu_cores` (Int, default=2): CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `weights_tarball` (File): Compressed tarball containing AlphaFold2 model weights (~15-20 GB)

### `colabfold_predict`

Predict protein structures from amino acid sequences using ColabFold (AlphaFold2 + MMseqs2).

**Inputs:**
- `fasta_file` (File): Input FASTA file containing one or more protein sequences
- `weights_tarball` (File): Model weights from `download_weights` task
- `output_prefix` (String): Prefix for naming the output tarball
- `num_recycle` (Int, default=3): Number of prediction recycles
- `num_models` (Int, default=5): Number of models to generate per sequence
- `num_seeds` (Int, default=1): Number of random seeds per model
- `msa_mode` (String, default="mmseqs2_uniref_env"): MSA mode (mmseqs2_uniref_env, mmseqs2_uniref, single_sequence)
- `use_amber` (Boolean, default=true): Enable AMBER relaxation
- `stop_at_score` (Int, default=85): Stop early at this pLDDT score
- `use_gpu_relax` (Boolean, default=true): Run AMBER relaxation on GPU
- `model_type` (String, default="auto"): Model type (auto, alphafold2, alphafold2_multimer_v3)
- `extra_flags` (String, default=""): Additional colabfold_batch flags
- `cpu_cores` (Int, default=8): CPU cores
- `memory_gb` (Int, default=48): Memory in GB
- `gpu_enabled` (Boolean, default=true): Enable GPU for prediction (controls JAX CPU/GPU mode and PROOF GPU allocation)

**Outputs:**
- `results_tarball` (File): Compressed tarball containing all ColabFold outputs (PDB files, PAE/pLDDT plots, coverage plots, prediction metrics)

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-colabfold/ww-colabfold.wdl" as colabfold_tasks

workflow protein_structure_prediction {
  input {
    File protein_fasta
    String sample_name
  }

  # Download weights once
  call colabfold_tasks.download_weights { }

  # Predict structures
  call colabfold_tasks.colabfold_predict {
    input:
      fasta_file = protein_fasta,
      weights_tarball = download_weights.weights_tarball,
      output_prefix = sample_name
  }

  output {
    File predicted_structures = colabfold_predict.results_tarball
  }
}
```

### Batch Prediction (Multiple Proteins)

```wdl
workflow batch_structure_prediction {
  input {
    Array[File] protein_fastas
    Array[String] sample_names
  }

  # Download weights once, reuse for all predictions
  call colabfold_tasks.download_weights { }

  scatter (idx in range(length(protein_fastas))) {
    call colabfold_tasks.colabfold_predict {
      input:
        fasta_file = protein_fastas[idx],
        weights_tarball = download_weights.weights_tarball,
        output_prefix = sample_names[idx]
    }
  }

  output {
    Array[File] predicted_structures = colabfold_predict.results_tarball
  }
}
```

### Advanced Usage Examples

**Single-sequence prediction (no MSA, for de novo designed proteins):**
```wdl
call colabfold_tasks.colabfold_predict {
  input:
    fasta_file = designed_protein_fasta,
    weights_tarball = download_weights.weights_tarball,
    output_prefix = "designed_protein",
    msa_mode = "single_sequence",
    num_recycle = 6
}
```

**Multimer/complex prediction:**
```wdl
call colabfold_tasks.colabfold_predict {
  input:
    fasta_file = complex_fasta,
    weights_tarball = download_weights.weights_tarball,
    output_prefix = "protein_complex",
    model_type = "alphafold2_multimer_v3"
}
```

**High-accuracy prediction with more models and recycles:**
```wdl
call colabfold_tasks.colabfold_predict {
  input:
    fasta_file = protein_fasta,
    weights_tarball = download_weights.weights_tarball,
    output_prefix = "high_accuracy",
    num_recycle = 6,
    num_models = 5,
    num_seeds = 3,
    stop_at_score = 95,
    memory_gb = 64
}
```

### Integration Examples

This module integrates with other WILDS components:
- **ww-diamond**: Use DIAMOND for sequence homology searches before structure prediction
- **ww-sra**: Download protein FASTAs from the Sequence Read Archive (SRA)
- **ww-ena**: Download protein FASTAs from the European Nucleotide Archive (ENA)
- **ww-aws-sso**: Interact with AWS S3 buckets

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that downloads weights and predicts the structure of a tiny peptide (Trp-cage, 20 residues) using CPU-only mode:

```bash
# Using miniWDL
miniwdl run modules/ww-colabfold/testrun.wdl

# Using Sprocket
sprocket run modules/ww-colabfold/testrun.wdl --entrypoint colabfold_example

# Using Cromwell
java -jar cromwell.jar run modules/ww-colabfold/testrun.wdl
```

### GPU Note

The test workflow runs on CPU (`gpu_enabled = false`) with minimal settings for CI compatibility. Production use should always enable GPU for practical performance. A single protein prediction that takes minutes on GPU can take hours on CPU.

## GPU Configuration

The `gpu_enabled` input controls two things:
1. **JAX execution mode**: When `false`, sets `JAX_PLATFORMS=cpu` to force CPU-only execution
2. **PROOF GPU allocation**: Translates to the `gpus` runtime attribute (`"1"` or `"0"`), which is specific to PROOF's Cromwell configuration. Other executors like miniWDL and Sprocket silently ignore this attribute.

### For PROOF Users
GPU allocation works automatically with the default (`gpu_enabled = true`). ColabFold only supports a single GPU for structure prediction.

### For Other Executors
Configure GPU passthrough at the engine level (e.g., Docker `--gpus` flag). Set `gpu_enabled = true` to ensure ColabFold uses the GPU via JAX.

### CPU-Only Mode
Set `gpu_enabled = false` to force CPU execution via `JAX_PLATFORMS=cpu`. This is functional but significantly slower and intended primarily for testing.

## Docker Container

This module uses the custom WILDS Docker image [`getwilds/colabfold:1.5.5`](https://github.com/getwilds/wilds-docker-library), built on `nvidia/cuda:11.8.0-base-ubuntu22.04` with Miniforge/mamba. It includes:
- ColabFold 1.5.5 (compatible with AlphaFold 2.3.2)
- `colabfold_batch` and `colabfold_search` CLI tools
- JAX 0.4.x with CUDA 11.8 support for GPU acceleration
- MMseqs2 for fast sequence searching
- AMBER tools for structure relaxation
- dm-haiku 0.0.12 (pinned for JAX 0.4.x compatibility)

AlphaFold2 model weights (~15-20 GB) are downloaded separately via the `download_weights` task and passed to `colabfold_predict` as input.

## Citation

If you use this module for published work, please cite ColabFold:

> Mirdita M, Schutze K, Moriwaki Y, Heo L, Ovchinnikov S, Steinegger M.
> ColabFold: making protein folding accessible to all.
> *Nature Methods* 19, 679-682 (2022).
> DOI: [10.1038/s41592-022-01488-1](https://doi.org/10.1038/s41592-022-01488-1)

## Parameters and Resource Requirements

### `download_weights` Resources
- **CPU**: 2 cores
- **Memory**: 8 GB
- **GPU**: Not required
- **Network**: Requires internet access to download ~15-20 GB of model weights
- **Runtime**: 5-15 minutes depending on network speed

### `colabfold_predict` Resources
- **CPU**: 8 cores
- **Memory**: 48 GB
- **GPU**: 1 GPU required for production use (NVIDIA, CUDA 11.8+ compatible)
- **Runtime**: Minutes per protein on GPU; hours on CPU

### Resource Scaling
- `cpu_cores`: Increase for faster MSA generation (8-16 recommended)
- `memory_gb`: Scale with protein size (48 GB for most proteins, 64+ GB for large complexes)
- `gpu_enabled`: Must be `true` for production use; `false` only for testing/debugging
- Protein length significantly impacts memory and runtime requirements


## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[ColabFold GitHub](https://github.com/sokrypton/ColabFold)**: Official ColabFold repository
- **[ColabFold Docker Wiki](https://github.com/sokrypton/ColabFold/wiki/Running-ColabFold-in-Docker)**: Docker usage guide
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
