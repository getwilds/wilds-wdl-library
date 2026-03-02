# ww-colabfold Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for predicting protein 3D structures from amino acid sequences using [ColabFold](https://github.com/sokrypton/ColabFold), which combines AlphaFold2 with the fast homology search of MMseqs2.

## Overview

This module wraps the `colabfold_batch` command-line tool for protein structure prediction. ColabFold provides a 40-60x speedup over standard AlphaFold2 by replacing the homology search step with MMseqs2, while matching AlphaFold2's prediction accuracy on CASP14 benchmarks. It supports monomeric proteins, multimers/complexes, and single-sequence predictions.

**Key features:**
- Fast MSA generation via MMseqs2 server (no local database required by default)
- AlphaFold2-quality structure predictions
- AMBER force field relaxation for refined structures
- Batch processing of multiple sequences
- Support for protein complexes via multimer mode

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-colabfold.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `colabfold_predict`

Predict protein structures from amino acid sequences using ColabFold (AlphaFold2 + MMseqs2).

**Inputs:**
- `fasta_file` (File): Input FASTA file containing one or more protein sequences
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
- `gpu_enabled` (Boolean, default=true): Whether GPU is available

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

  call colabfold_tasks.colabfold_predict {
    input:
      fasta_file = protein_fasta,
      output_prefix = sample_name
  }

  output {
    File predicted_structures = colabfold_predict.results_tarball
  }
}
```

### Advanced Usage Examples

**Single-sequence prediction (no MSA, for de novo designed proteins):**
```wdl
call colabfold_tasks.colabfold_predict {
  input:
    fasta_file = designed_protein_fasta,
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
    output_prefix = "protein_complex",
    model_type = "alphafold2_multimer_v3"
}
```

**High-accuracy prediction with more models and recycles:**
```wdl
call colabfold_tasks.colabfold_predict {
  input:
    fasta_file = protein_fasta,
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
- **Other WILDS modules**: Can be used downstream of sequence analysis workflows

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that predicts the structure of a tiny peptide (Trp-cage, 20 residues) using CPU-only mode:

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

## Docker Container

This module uses the official ColabFold Docker image [`ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2`](https://github.com/sokrypton/ColabFold/pkgs/container/colabfold), which includes:
- ColabFold 1.5.5 (compatible with AlphaFold 2.3.2)
- `colabfold_batch` and `colabfold_search` CLI tools
- CUDA 12.2.2 runtime for GPU acceleration
- MMseqs2 for fast sequence searching
- AMBER tools for structure relaxation
- AlphaFold2 model weights are downloaded at runtime

## Citation

If you use this module for published work, please cite ColabFold:

> Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S, Steinegger M.
> ColabFold: making protein folding accessible to all.
> *Nature Methods* 19, 679–682 (2022).
> DOI: [10.1038/s41592-022-01488-1](https://doi.org/10.1038/s41592-022-01488-1)

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 8 cores
- **Memory**: 48 GB
- **GPU**: 1 GPU required (NVIDIA, CUDA 12.2+ compatible)
- **Runtime**: Minutes per protein on GPU; hours on CPU

### Resource Scaling
- `cpu_cores`: Increase for faster MSA generation (8-16 recommended)
- `memory_gb`: Scale with protein size (48 GB for most proteins, 64+ GB for large complexes)
- `gpu_enabled`: Must be `true` for production use; `false` only for testing/debugging
- Protein length significantly impacts memory and runtime requirements

### Recommended GPU Hardware
- NVIDIA A100 (preferred for large proteins)
- NVIDIA V100 or T4 (suitable for most predictions)
- Minimum 16 GB GPU memory

## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[ColabFold GitHub](https://github.com/sokrypton/ColabFold)**: Official ColabFold repository
- **[ColabFold Docker Wiki](https://github.com/sokrypton/ColabFold/wiki/Running-ColabFold-in-Docker)**: Docker usage guide
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
