# ww-starling-batch Pipeline

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline for batch structural ensemble generation of intrinsically disordered proteins using [STARLING](https://github.com/idptools/starling). This pipeline splits a large multi-sequence FASTA file into batches and processes them in parallel.

## Overview

This pipeline enables high-throughput ensemble generation for large sets of intrinsically disordered protein sequences. It splits the input FASTA into configurable batch sizes, scatters ensemble generation across batches for parallel processing, and collects all outputs.

**Complexity Level**: Basic (1 module)

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and demonstrates:

- **Scatter-Gather**: Splitting input data into batches for parallel processing
- **Native Batching**: Leveraging STARLING's built-in multi-sequence FASTA support within each batch
- **Resource Efficiency**: Balancing per-batch parallelism with within-batch native processing

## Pipeline Steps

1. **FASTA Splitting** (using `ww-starling` module):
   - Splits the input multi-sequence FASTA into smaller batch files
   - Each batch contains up to `sequences_per_batch` sequences

2. **Parallel Ensemble Generation** (using `ww-starling` module):
   - Scatters `generate_ensemble_batch` across all batch files
   - Each batch processes its sequences natively using STARLING's multi-sequence support
   - Generates `.starling` ensemble files plus PDB and XTC format conversions

## Module Dependencies

This pipeline imports and uses:
- **ww-starling module**: For FASTA splitting (`split_fasta`), batch ensemble generation (`generate_ensemble_batch`)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for module imports

### Input Configuration

Create an inputs JSON file with your FASTA file and batch parameters:

```json
{
  "starling_batch.fasta_file": "/path/to/protein_sequences.fasta",
  "starling_batch.sequences_per_batch": 10,
  "starling_batch.num_conformations": 400,
  "starling_batch.cpu_cores": 4,
  "starling_batch.memory_gb": 8
}
```

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-starling-batch.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-starling-batch.wdl -i inputs.json

# Using Sprocket
sprocket run ww-starling-batch.wdl inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this pipeline directly to the on-premise HPC cluster.

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `fasta_file` | Input FASTA file containing protein sequences | File | Yes | - |
| `sequences_per_batch` | Number of sequences per batch for parallel processing | Int | No | 10 |
| `num_conformations` | Number of conformations to generate per sequence | Int | No | 400 |
| `gpu_enabled` | Enable GPU for STARLING inference in each batch task | Boolean | No | true |
| `cpu_cores` | CPU cores allocated per batch task | Int | No | 4 |
| `memory_gb` | Memory in GB allocated per batch task | Int | No | 8 |

## Output Files

| Output | Description | Source Module |
|--------|-------------|---------------|
| `starling_files` | STARLING ensemble files (one per sequence, grouped by batch) | ww-starling |
| `pdb_files` | PDB topology files (one per sequence, grouped by batch) | ww-starling |
| `xtc_files` | XTC trajectory files (one per sequence, grouped by batch) | ww-starling |

## Resource Considerations

### Compute Requirements
- **Memory**: 8 GB per batch recommended; increase for very long sequences (>200 residues)
- **CPUs**: 4 cores per batch recommended
- **Storage**: Each ensemble generates `.starling`, `.pdb`, and `.xtc` files per sequence

### Optimization Tips
- Adjust `sequences_per_batch` to balance parallelism vs. task overhead (smaller batches = more parallelism, larger batches = less overhead)
- For short sequences (<50 residues), larger batches (20-50) are efficient
- For long sequences (>200 residues), smaller batches (2-5) with more memory are recommended
- Use call caching to avoid reprocessing if the pipeline is interrupted

## Testing the Pipeline

The pipeline includes a test workflow using well-known IDP sequences:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint starling_batch_example
```

The test workflow automatically:
1. Creates a test FASTA with 4 well-known IDP sequences (p53 NTAD, ASH1, NUPR1, p27Kip1) using `ww-testdata`
2. Calls the `starling_batch` pipeline with batches of 2, 50 conformations, and GPU disabled for CI compatibility
3. The pipeline splits the FASTA and scatters ensemble generation across 2 batches in parallel

## Citation

> Lotthammer, J.M., Ginell, G.M., Griffith, D., Emenecker, R.J., & Holehouse, A.S. (2024)
> Direct prediction of intrinsically disordered protein conformational properties from sequences.
> *Nature Methods*, 21, 465-476.
> DOI: [10.1038/s41592-023-02159-5](https://doi.org/10.1038/s41592-023-02159-5)

## Related WILDS Components

- **ww-starling module**: STARLING ensemble generation tasks
- **ww-testdata module**: Test data provisioning

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
