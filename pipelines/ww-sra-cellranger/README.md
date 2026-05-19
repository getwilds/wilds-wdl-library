# ww-sra-cellranger Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline for downloading single-cell RNA-seq data from SRA and processing with Cell Ranger count.

## Overview

This pipeline combines the `ww-sra` and `ww-cellranger` modules to download scRNA-seq FASTQ data from NCBI's Sequence Read Archive and run Cell Ranger `count` for gene expression quantification. It supports both public SRA data and controlled-access dbGaP data via optional NGC authentication.

The pipeline automatically handles FASTQ renaming to satisfy Cell Ranger's strict naming convention, so users only need to provide SRA accession IDs and a reference transcriptome.

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and demonstrates:

- **Module Integration**: Combining `ww-sra` and `ww-cellranger` modules
- **Data Flow**: SRA download, FASTQ renaming, and Cell Ranger processing
- **Complexity Level**: Basic (2 modules)

## Pipeline Steps

1. **SRA Download** (using `ww-sra` module):
   - Downloads FASTQ files from SRA accessions using `prefetch` + `fasterq-dump`
   - Supports controlled-access dbGaP data via optional NGC key file

2. **FASTQ Renaming** (using `ww-cellranger` module):
   - Renames downloaded FASTQs to Cell Ranger's required naming convention
   - Converts `{sra_id}_1.fastq.gz` to `{sra_id}_S1_R1_001.fastq.gz`

3. **Cell Ranger Count** (using `ww-cellranger` module):
   - Runs `cellranger count` on each sample via one of `run_count` (private Docker image), `run_count_hpc_cromwell` (HPC environment module on a Cromwell-on-HPC backend), or `run_count_hpc_sprocket` (HPC environment module inside a Sprocket-managed Lua container), selected by the `execution_mode` input
   - Generates gene expression matrices, web summaries, and metrics

## Module Dependencies

This pipeline imports and uses:
- **ww-sra module**: For SRA data download (`fastqdump` task, with optional dbGaP/NGC support)
- **ww-cellranger module**: For FASTQ renaming (`rename_fastqs` task) and gene expression quantification (one of `run_count`, `run_count_hpc_cromwell`, or `run_count_hpc_sprocket`, selected per run via the `execution_mode` input)

## Cell Ranger Software Environment

Cell Ranger is not redistributable, so the WILDS Docker Library does not publish a public Cell Ranger image. This pipeline has three `execution_mode` options that run a corresponding version of the `ww-cellranger/run_count` task:

| `execution_mode` | Task called | Use when |
| --- | --- | --- |
| `"docker"` (default) | `run_count` | You have a private Cell Ranger Docker image (built from the [WILDS Dockerfile recipe](https://github.com/getwilds/wilds-docker-library/blob/main/cellranger/Dockerfile_latest) or supplied by your organization). Override `docker_image` to point at it. |
| `"hpc_cromwell"` | `run_count_hpc_cromwell` | You are running on a Cromwell-on-HPC backend that executes tasks directly on the compute node (e.g., Fred Hutch HPC via PROOF). Cell Ranger is loaded via `module load` on the host and made available on `PATH`. Override `cellranger_module` if your site uses a different module name/version. |
| `"hpc_sprocket"` | `run_count_hpc_sprocket` | You are running on a Sprocket-on-HPC backend, which always runs tasks inside a container under Apptainer. Cell Ranger is loaded inside a minimal Lua container with the host's Lmod and Cell Ranger software trees bind-mounted in (see [`.github/configs/sprocket-hpc.toml`](../../.github/configs/sprocket-hpc.toml)). Override `cellranger_module` if your site uses a different module name/version. |

Unused inputs are ignored — e.g., `docker_image` has no effect when `execution_mode = "hpc_sprocket"`. An invalid `execution_mode` causes all dispatch branches to skip and the pipeline's `select_first` calls to fail loudly. See the [ww-cellranger module README](../../modules/ww-cellranger/README.md) for more detail on the three task variants.

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- One of: Docker/Apptainer support with a private Cell Ranger image (default `execution_mode = "docker"`), a Cromwell-on-HPC backend that runs tasks directly on the compute node and provides Cell Ranger via Lmod (`execution_mode = "hpc_cromwell"`), or a Sprocket-on-HPC config that bind-mounts the host's Lmod and Cell Ranger software trees into the container (`execution_mode = "hpc_sprocket"`)
- Internet access for SRA downloads
- Sufficient compute resources for Cell Ranger (64GB+ RAM recommended)
- **Platform**: Cell Ranger requires Linux x86_64 with AVX support (not compatible with Apple Silicon)

### Input Configuration

Create an inputs JSON file with your SRA accessions and Cell Ranger reference:

```json
{
  "sra_cellranger.sra_id_list": ["SRR12345678"],
  "sra_cellranger.ref_gex": "/path/to/cellranger/reference.tar.gz",
  "sra_cellranger.ncpu": 8,
  "sra_cellranger.memory_gb": 64,
  "sra_cellranger.execution_mode": "hpc_cromwell",
  "sra_cellranger.docker_image": "ghcr.io/getwilds/cellranger:10.0.0",
  "sra_cellranger.cellranger_module": "CellRanger/10.0.0"
}
```

> **Note:** This template sets `execution_mode: "hpc_cromwell"` so it works out of the box for Fred Hutch users on PROOF (where Cell Ranger is provided as the `CellRanger/10.0.0` environment module under Cromwell). The WDL default is `"docker"` — if you are running on a container-based backend, change it to `"docker"` and override `docker_image` to point at your private Cell Ranger image. If you are submitting to Sprocket-on-HPC instead of Cromwell-on-HPC, change it to `"hpc_sprocket"`. See [Cell Ranger Software Environment](#cell-ranger-software-environment) above for details.
>
> The `ngc_file` parameter is optional. Include it when downloading controlled-access dbGaP data.

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-sra-cellranger.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-sra-cellranger.wdl -i inputs.json

# Using Sprocket
sprocket run ww-sra-cellranger.wdl @inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/datademos/proof-how-to/) to submit this pipeline directly to the on-premise HPC cluster. For controlled-access dbGaP data, use [PROOF Regulated](https://sciwiki.fredhutch.org/datademos/proof-regulated/) to ensure compliance with data use agreements.

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `sra_id_list` | List of SRA accession IDs to download | Array[String] | Yes | - |
| `ref_gex` | Cell Ranger GEX reference transcriptome tarball | File | Yes | - |
| `ncpu` | Number of CPU cores | Int | No | 8 |
| `memory_gb` | Memory allocation in GB | Int | No | 64 |
| `max_reads` | Maximum reads to download per sample (for testing) | Int | No | all reads |
| `ngc_file` | NGC repository key file for controlled-access dbGaP data | File | No | - |
| `create_bam` | Whether Cell Ranger should generate a BAM file | Boolean | No | true |
| `expect_cells` | Expected number of recovered cells per sample | Int | No | - |
| `chemistry` | Assay configuration (e.g., SC3Pv2, SC3Pv3) | String | No | auto-detect |
| `execution_mode` | Which Cell Ranger task to dispatch to: `"docker"`, `"hpc_cromwell"`, or `"hpc_sprocket"`. See [Cell Ranger Software Environment](#cell-ranger-software-environment). | String | No | `"docker"` |
| `docker_image` | Private Cell Ranger Docker image used by `run_count`. Ignored unless `execution_mode = "docker"`. | String | No | `ghcr.io/getwilds/cellranger:10.0.0` |
| `cellranger_module` | HPC environment module used by the `run_count_hpc_*` tasks. Ignored unless `execution_mode` starts with `hpc_`. | String | No | `CellRanger/10.0.0` |

### Cell Ranger Reference

Cell Ranger requires a pre-built reference transcriptome tarball. You can:
- Download pre-built references from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads)
- Build custom references using `cellranger mkref`

## Output Files

| Output | Description | Source Module |
|--------|-------------|---------------|
| `cellranger_results` | Compressed tarballs of Cell Ranger count output directories | ww-cellranger |
| `cellranger_web_summaries` | Web summary HTML files for each sample | ww-cellranger |
| `cellranger_metrics` | Metrics summary CSV files for each sample | ww-cellranger |

## Resource Considerations

### Compute Requirements
- **Memory**: 64GB+ recommended for human reference (6GB minimum for small test datasets)
- **CPUs**: 8+ cores recommended for efficient processing
- **Storage**: Sufficient space for SRA downloads, FASTQ files, and Cell Ranger outputs
- **Network**: Stable internet connection for SRA downloads

### Platform Requirements
- Cell Ranger requires Linux x86_64 with AVX instruction support
- Works natively on Linux x86_64 systems
- Works on Intel-based Macs via Docker
- **Not supported** on Apple Silicon (M1/M2/M3/M4) Macs due to lack of AVX support in Docker's x86_64 emulation

## Testing the Pipeline

The pipeline includes a test workflow that automatically runs with minimal test data:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

The test workflow automatically:
1. Downloads a small Cell Ranger reference
2. Downloads a 10x Chromium 3' v2 scRNA-seq dataset from SRA (SRR9169219, limited to 100k reads)
3. Renames FASTQs for Cell Ranger compatibility
4. Runs Cell Ranger count
5. Outputs results, web summary, and metrics

### HPC Test Workflow

`testrun_hpc.wdl` mirrors the regular testrun's coverage but dispatches with `execution_mode = "hpc_sprocket"` so it exercises the Sprocket + module-load path end-to-end. It is intended to be exercised on Fred Hutch HPC via Sprocket on a monthly basis to validate the HPC dispatch path; the Cromwell-on-HPC path (`execution_mode = "hpc_cromwell"`) is exercised by users running Cromwell on their HPC by hand. Uses the same tiny SRA download and Cell Ranger reference as `testrun.wdl`.

## Related WILDS Components

- **ww-sra module**: SRA download functionality
- **ww-cellranger module**: Cell Ranger processing functionality
- **ww-sra-star pipeline**: Alternative pipeline for bulk RNA-seq alignment
- **ww-sra-salmon pipeline**: Alternative pipeline for bulk RNA-seq quantification

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
