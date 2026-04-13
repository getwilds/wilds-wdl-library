# ww-esmfold Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for predicting protein 3D structures using ESMFold. ESMFold leverages the ESM-2 protein language model for fast, single-sequence structure prediction without requiring multiple sequence alignments (MSAs), making it approximately 60x faster than AlphaFold2-based methods.

## Overview

ESMFold predicts protein structures directly from amino acid sequences using a large protein language model. Unlike AlphaFold2, it does not require MSA computation, which dramatically reduces prediction time. The predicted structures include pLDDT confidence scores in the B-factor column of the output PDB files.

**Key advantages over MSA-based methods:**
- No MSA or template database required
- Much faster inference (~60x compared to AlphaFold2)
- Single-sequence input is sufficient

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-esmfold.wdl` - Contains task definitions for ESMFold prediction
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `esmfold_predict`

Predict protein 3D structures from amino acid sequences using ESMFold.

**Inputs:**
- `fasta_file` (File): Input FASTA file containing one or more protein sequences to predict
- `output_prefix` (String): Prefix for naming the output tarball
- `num_recycles` (Int, default=4): Number of recycling iterations for structure refinement
- `chunk_size` (Int, default=128): Chunk size for memory-efficient inference on long sequences (0 to disable)
- `cpu_only` (Boolean, default=false): Run inference on CPU instead of GPU
- `cpu_offload` (Boolean, default=false): Enable CPU offloading for longer sequences with limited GPU memory
- `max_tokens_per_batch` (Int, default=0): Maximum tokens per batch for GPU inference (0 for automatic)
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=16): Memory allocated for the task in GB
- `gpu_enabled` (Boolean, default=true): Enable GPU for prediction

**Outputs:**
- `pdb_output` (File): Compressed tarball containing predicted PDB structure files with pLDDT confidence scores

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-esmfold/ww-esmfold.wdl" as esmfold_tasks

workflow my_structure_prediction {
  input {
    File protein_fasta
  }

  call esmfold_tasks.esmfold_predict { input:
      fasta_file = protein_fasta,
      output_prefix = "my_prediction"
  }

  output {
    File predicted_structures = esmfold_predict.pdb_output
  }
}
```

### Advanced Usage Examples

**CPU-only prediction for testing:**
```wdl
call esmfold_tasks.esmfold_predict { input:
    fasta_file = protein_fasta,
    output_prefix = "cpu_test",
    cpu_only = true,
    gpu_enabled = false,
    num_recycles = 1,
    cpu_cores = 4,
    memory_gb = 16
}
```

**Long sequence prediction with memory optimization:**
```wdl
call esmfold_tasks.esmfold_predict { input:
    fasta_file = long_protein_fasta,
    output_prefix = "long_sequence",
    chunk_size = 64,
    cpu_offload = true,
    memory_gb = 32
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of test protein sequences for demonstrations
- **ww-colabfold**: Alternative structure prediction using AlphaFold2 + MMseqs2

## Testing the Module

> **Note:** The ESMFold test workflow requires **24 GB of memory** to load the full 3B-parameter ESM-2 model, which exceeds the resources available on GitHub Actions runners (~16 GB). As a result, **this module's `testrun.wdl` is excluded from GitHub CI/CD** and is instead validated on the Fred Hutch HPC on a monthly basis. Linting checks still run in CI as normal.

The module includes a test workflow (`testrun.wdl`) that can be run on an HPC or local machine with sufficient memory:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint esmfold_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Creates a minimal test protein FASTA (Trp-cage miniprotein, 20 residues) using `ww-testdata`
2. Runs ESMFold prediction with minimal settings (CPU-only, 1 recycle)
3. Validates that PDB output files were generated

## Docker Container

This module uses the `getwilds/esmfold:2.0.0` container image, which includes:
- ESM-2 protein language model via the `fair-esm` package
- OpenFold dependencies for structure module
- All necessary Python dependencies for inference

## Citation

> Lin, Z., Akin, H., Rao, R., Hie, B., Zhu, Z., Lu, W., ... & Rives, A. (2023).
> Evolutionary-scale prediction of atomic-level protein structure with a language model.
> Science, 379(6637), 1123-1130.
> DOI: https://doi.org/10.1126/science.ade2574

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 4 cores
- **Memory**: 16 GB
- **GPU**: 1 GPU recommended for production use

### Resource Scaling
- `cpu_cores`: Increase for CPU-only mode or parallel processing
- `memory_gb`: Increase for longer protein sequences (>500 residues may need 32+ GB)
- `gpu_enabled`: Set to `true` for production use; `false` for testing
- `chunk_size`: Reduce (64 or 32) for very long sequences to manage memory usage
- `cpu_offload`: Enable for sequences that exceed GPU memory

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for detailed guidelines

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
- **[ESMFold GitHub](https://github.com/facebookresearch/esm)**: Official ESMFold source code and documentation
