# ww-starling Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for [STARLING](https://github.com/idptools/starling), a latent-space probabilistic denoising diffusion model for predicting coarse-grained structural ensembles of intrinsically disordered protein (IDP) regions.

## Overview

STARLING uses a Vision Transformer architecture operating in VAE latent space to generate 3D structural conformations from amino acid sequences. This module wraps STARLING's core functionality for ensemble generation and format conversion, making it available as a reusable component in WDL workflows.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-starling.wdl` - Contains task definitions for ensemble generation and metadata querying
- **Test workflows**: `testrun.wdl` for CI/CD (CPU-only); `testrun_hpc.wdl` for monthly HPC validation with GPU enabled and the module's default 400 conformations
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `generate_ensemble`

Generate a structural ensemble for an intrinsically disordered protein sequence.

**Inputs:**
- `sequence` (String): Amino acid sequence string for the disordered protein region
- `sample_name` (String): Name identifier for the output files
- `num_conformations` (Int, default=400): Number of conformations to generate in the ensemble
- `gpu_enabled` (Boolean, default=true): Enable GPU for STARLING inference (sets device to cuda and requests GPU in runtime). See [GPU Execution Across Environments](#gpu-execution-across-environments) for backend-specific configuration.
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `starling_file` (File): STARLING ensemble file containing predicted conformations
- `pdb_file` (File): PDB topology file converted from the STARLING ensemble
- `xtc_file` (File): XTC trajectory file converted from the STARLING ensemble

### `generate_ensemble_batch`

Generate structural ensembles for multiple protein sequences from a FASTA file.

**Inputs:**
- `fasta_file` (File): FASTA file containing one or more protein sequences
- `num_conformations` (Int, default=400): Number of conformations to generate per sequence in the ensemble
- `gpu_enabled` (Boolean, default=true): Enable GPU for STARLING inference (sets device to cuda and requests GPU in runtime). See [GPU Execution Across Environments](#gpu-execution-across-environments) for backend-specific configuration.
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `starling_files` (Array[File]): Array of STARLING ensemble files, one per sequence
- `pdb_files` (Array[File]): Array of PDB topology files
- `xtc_files` (Array[File]): Array of XTC trajectory files

### `split_fasta`

Split a multi-sequence FASTA file into smaller batch files for parallel processing.

**Inputs:**
- `fasta_file` (File): Input FASTA file containing multiple protein sequences
- `sequences_per_batch` (Int, default=10): Maximum number of sequences per output batch file
- `cpu_cores` (Int, default=1): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=2): Memory allocated for the task in GB

**Outputs:**
- `batch_files` (Array[File]): Array of FASTA files, each containing up to `sequences_per_batch` sequences

### `ensemble_info`

Query metadata and summary information from a STARLING ensemble file.

**Inputs:**
- `starling_file` (File): STARLING ensemble file to query
- `sample_name` (String): Name identifier for the output file
- `cpu_cores` (Int, default=1): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=4): Memory allocated for the task in GB

**Outputs:**
- `info_file` (File): Text file containing ensemble metadata and summary statistics

## Usage as a Module

### Importing into Your Workflow

The example below imports from `refs/heads/main`, which always points at the latest version of the module. For reproducible workflows, you can pin the import to a specific release tag (`refs/tags/v0.3.0`) or commit SHA (`<full-sha>`) by swapping `refs/heads/main` in the URL — just make sure the tag or commit you pin to actually contains this module.

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-starling/ww-starling.wdl" as starling_tasks

workflow idp_analysis {
  input {
    String sequence
    String sample_name
  }

  call starling_tasks.generate_ensemble { input:
    sequence = sequence,
    sample_name = sample_name
  }

  call starling_tasks.ensemble_info { input:
    starling_file = generate_ensemble.starling_file,
    sample_name = sample_name
  }

  output {
    File starling_file = generate_ensemble.starling_file
    File pdb_file = generate_ensemble.pdb_file
    File xtc_file = generate_ensemble.xtc_file
    File info_file = ensemble_info.info_file
  }
}
```

### Advanced Usage Examples

**Custom ensemble size and resources:**
```wdl
call starling_tasks.generate_ensemble { input:
  sequence = "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKES",
  sample_name = "alpha_synuclein_nterm",
  num_conformations = 1000,
  cpu_cores = 8,
  memory_gb = 16
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of test protein sequences for demonstrations
- **Other WILDS modules**: Can be used alongside structure prediction tools like ColabFold for comparative analysis

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint starling_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow exercises all four module tasks:
1. Generates a single-sequence ensemble for the p53 N-terminal transactivation domain (residues 1-39) using `generate_ensemble`
2. Queries ensemble metadata with `ensemble_info`
3. Creates a multi-sequence IDP FASTA (4 sequences) using `ww-testdata` and processes it with `generate_ensemble_batch`
4. Splits the multi-sequence FASTA into batches of 2 using `split_fasta`

### HPC Test Workflow

`testrun_hpc.wdl` mirrors the regular testrun's coverage but enables GPU on the inference tasks (`generate_ensemble`, `generate_ensemble_batch`) and uses the module's default 400 conformations instead of the CI-trimmed 50. It is intended to be exercised on Fred Hutch HPC (via Sprocket directly or via PROOF after applying the WDL 1.0 / `gpus` swap described above) and includes a validation task that confirms expected ensemble, conversion, info, and split outputs are present.

## Docker Container

This module uses the `getwilds/starling:2.0.0a3` container image, which should include:
- STARLING (`idptools-starling`) and its dependencies
- `starling2pdb`, `starling2xtc`, and `starling2info` CLI utilities

## Citation

> Lotthammer, J.M., Ginell, G.M., Griffith, D., Emenecker, R.J., & Holehouse, A.S. (2024)
> Direct prediction of intrinsically disordered protein conformational properties from sequences.
> *Nature Methods*, 21, 465-476.
> DOI: [10.1038/s41592-023-02159-5](https://doi.org/10.1038/s41592-023-02159-5)

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 4 cores (generate_ensemble), 1 core (ensemble_info)
- **Memory**: 8 GB (generate_ensemble), 4 GB (ensemble_info)
- **Runtime**: Varies with sequence length and number of conformations

### Resource Scaling
- `cpu_cores`: Increase for longer sequences or larger ensembles
- `memory_gb`: Increase for very long disordered regions (>200 residues)
- `num_conformations`: More conformations provide better sampling but take longer

## GPU Execution Across Environments

The `gpu_enabled` input controls two things in the GPU-using tasks (`generate_ensemble` and `generate_ensemble_batch`):
1. **STARLING device selection**: When `true`, passes `-d cuda` to STARLING; when `false`, passes `-d cpu`.
2. **WDL runtime GPU attribute**: Translates to the standard `gpu: Boolean` key in the `runtime` block.

GPU scheduling is not standardized across WDL executors, so the way GPUs are requested in the `runtime` section depends on where the workflow is run. This module ships configured for the standard `gpu: Boolean` key under WDL 1.2, which works with recent versions of miniWDL, Sprocket, and Cromwell in cloud/local environments — and with Sprocket on the Fred Hutch HPC.

If you'd rather run through PROOF (the Fred Hutch point-and-click interface for Cromwell on HPC), two modifications are required because PROOF's Cromwell deployment targets WDL 1.0:

1. **Downgrade the WDL version** from `version 1.2` to `version 1.0` at the top of `ww-starling.wdl`.
2. **Switch the runtime key** from `gpu` to `gpus` (an integer count passed as a string) in both `generate_ensemble` and `generate_ensemble_batch`. Each task already includes the alternate line as a comment for convenience:

   ```wdl
   runtime {
     docker: "getwilds/starling:2.0.0a3"
     # gpu: gpu_enabled
     gpus: if gpu_enabled then "1" else "0"
     cpu: cpu_cores
     memory: "~{memory_gb} GB"
   }
   ```

If you're on Fred Hutch HPC and don't need the PROOF UI, you can keep the module as-is and submit via Sprocket. For other HPC or cloud backends, check that backend's documentation for its expected GPU runtime attribute — some engines also accept `gpuCount`, `gpuType`, or backend-specific keys via `hints`.

## Contributing

1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for detailed guidelines

## Related Resources

- **[STARLING GitHub](https://github.com/idptools/starling)**: Source code and documentation
- **[STARLING Docs](https://idptools-starling.readthedocs.io)**: Full documentation and tutorials
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
