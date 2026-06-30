# ww-sra-cellranger Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline for downloading single-cell RNA-seq data from SRA, processing with Cell Ranger count, and removing ambient RNA with CellBender.

## Overview

This pipeline combines the `ww-sra`, `ww-cellranger`, and `ww-cellbender` modules to download scRNA-seq FASTQ data from NCBI's Sequence Read Archive, run Cell Ranger `count` for gene expression quantification, and remove ambient RNA contamination from the resulting count matrices. It supports both public SRA data and controlled-access dbGaP data via optional NGC authentication.

The pipeline automatically handles FASTQ renaming to satisfy Cell Ranger's strict naming convention, so users only need to provide SRA accession IDs and a reference transcriptome.

## Scope and Data Governance

This pipeline is scoped as an on-ramp for users who are new to computational workflows, but want Cell Ranger count outputs from SRA samples. To match that scope, BAM output is disabled: the pipeline hardcodes `create_bam = false` on all Cell Ranger calls and does not expose it as a workflow input. This avoids producing the largest Cell Ranger output for users whose downstream tooling (Seurat, Scanpy) reads directly from the `.h5` or matrix files, and avoids accidental persistence of BAMs derived from controlled-access dbGaP data, where redistribution restrictions typically prohibit sharing them.

When pulling controlled-access data from dbGaP via the `ngc_file` input, you are bound by the data use agreement attached to that study, which typically restricts where derived data (FASTQs, count matrices, BAMs, etc.) may be stored. **Make sure to run this pipeline in a location approved for regulated data storage** and avoid persisting outputs to general-purpose or shared filesystems.

**Fred Hutch users:** Use [PROOF Regulated](https://sciwiki.fredhutch.org/datademos/proof-regulated/) to submit this pipeline. PROOF Regulated stages analysis data under `/fh/regulated`, which is set up for compliance with dbGaP and similar data use agreements.

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and demonstrates:

- **Module Integration**: Combining `ww-sra`, `ww-cellranger`, and `ww-cellbender` modules
- **Data Flow**: SRA download, FASTQ renaming, Cell Ranger processing, and ambient RNA removal
- **Complexity Level**: Basic (3 modules)

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

4. **Ambient RNA Removal** (using `ww-cellbender` module):
   - Runs CellBender `remove-background` on each sample that produced a raw feature-barcode matrix
   - Samples skipped by Cell Ranger (chemistry detection failure) are automatically excluded
   - GPU acceleration is enabled by default; set `cellbender_gpu_enabled = false` for CPU-only execution

## Module Dependencies

This pipeline imports and uses:
- **ww-sra module**: For SRA data download (`fastqdump` task, with optional dbGaP/NGC support)
- **ww-cellranger module**: For FASTQ renaming (`rename_fastqs` task) and gene expression quantification (one of `run_count`, `run_count_hpc_cromwell`, or `run_count_hpc_sprocket`, selected per run via the `execution_mode` input)
- **ww-cellbender module**: For ambient RNA removal (`remove_background` task), run on each sample that produced a raw feature-barcode matrix

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
- Sufficient compute resources for Cell Ranger (64GB+ RAM recommended) and CellBender (32GB+ RAM; GPU strongly recommended, set `cellbender_gpu_enabled = false` for CPU-only)
- **Platform**: Cell Ranger requires Linux x86_64 with AVX support (not compatible with Apple Silicon)

### Input Configuration

Create an inputs JSON file with your SRA accessions and Cell Ranger reference. Provide the samples either as a text file of accessions (`sra_id_file`, shown below) or as an inline array (`sra_id_list`):

```json
{
  "sra_cellranger.sra_id_file": "/path/to/SRR_Acc_List.txt",
  "sra_cellranger.ref_gex": "/path/to/cellranger/reference.tar.gz",
  "sra_cellranger.ncpu": 8,
  "sra_cellranger.memory_gb": 64,
  "sra_cellranger.execution_mode": "hpc_cromwell",
  "sra_cellranger.docker_image": "ghcr.io/getwilds/cellranger:10.0.0",
  "sra_cellranger.cellranger_module": "CellRanger/10.0.0",
  "sra_cellranger.cellbender_gpu_enabled": true,
  "sra_cellranger.cellbender_expected_cells": 5000
}
```

> **Note:** This template sets `execution_mode: "hpc_cromwell"` so it works out of the box for Fred Hutch users on PROOF (where Cell Ranger is provided as the `CellRanger/10.0.0` environment module under Cromwell). The WDL default is `"docker"` — if you are running on a container-based backend, change it to `"docker"` and override `docker_image` to point at your private Cell Ranger image. If you are submitting to Sprocket-on-HPC instead of Cromwell-on-HPC, change it to `"hpc_sprocket"`. See [Cell Ranger Software Environment](#cell-ranger-software-environment) above for details.
>
> The `ngc_file` parameter is optional. Include it when downloading controlled-access dbGaP data.

### Selecting Samples with SRA Run Selector

`sra_id_file` is the one-accession-per-line text file produced by the **Accession List** button in NCBI's [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/). For any study with more than a handful of samples, this is the recommended way to provide inputs. It should look like this:

```
SRR7722937
SRR1039508
SRR13777504
```

One accession per line, with no header row, commas, or quotes: exactly what the **Accession List** button exports (scroll to the **Select** table > find the **Download** column > click the **Accession List** button).

**Curate your list in Run Selector before running.** Many SRA studies group multiple assay types under a single project: it is common to find single-cell RNA-seq, multiome RNA and ATAC, bulk RNA-seq, and whole-exome runs all under one accession. Use Run Selector's filters to select only the runs you intend to feed to Cell Ranger, then export that subset as your accession list. If a non-single-cell run does slip through, `skip_on_chemistry_failure` is the safety net rather than a substitute for curation. See [Mixed single-cell / non-single-cell input](#mixed-single-cell--non-single-cell-input) below.

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
| `sra_id_file` | Text file of SRA accession IDs, one per line | File | One of `sra_id_file` / `sra_id_list` | - |
| `sra_id_list` | List of SRA accession IDs to download. Alternative to `sra_id_file`. | Array[String] | One of `sra_id_file` / `sra_id_list` | - |
| `ref_gex` | Cell Ranger GEX reference transcriptome tarball | File | Yes | - |
| `ncpu` | Number of CPU cores | Int | No | 8 |
| `memory_gb` | Memory allocation in GB | Int | No | 64 |
| `max_reads` | Maximum reads to download per sample (for testing) | Int | No | all reads |
| `ngc_file` | NGC repository key file for controlled-access dbGaP data | File | No | - |
| `expect_cells` | Expected number of recovered cells per sample | Int | No | - |
| `chemistry` | Assay configuration (e.g., SC3Pv2, SC3Pv3) | String | No | auto-detect |
| `skip_on_chemistry_failure` | If true, samples Cell Ranger can't auto-detect a chemistry for are skipped instead of failing the workflow. See the module [README](../../modules/ww-cellranger/README.md#graceful-chemistry-detection-skip). | Boolean | No | `false` |
| `execution_mode` | Which Cell Ranger task to dispatch to: `"docker"`, `"hpc_cromwell"`, or `"hpc_sprocket"`. See [Cell Ranger Software Environment](#cell-ranger-software-environment). | String | No | `"docker"` |
| `docker_image` | Private Cell Ranger Docker image used by `run_count`. Ignored unless `execution_mode = "docker"`. | String | No | `ghcr.io/getwilds/cellranger:10.0.0` |
| `cellranger_module` | HPC environment module used by the `run_count_hpc_*` tasks. Ignored unless `execution_mode` starts with `hpc_`. | String | No | `CellRanger/10.0.0` |
| `organize_results` | When true, package all Cell Ranger outputs into a tarball organized by sample subdirectory. | Boolean | No | `false` |
| `output_prefix` | Prefix for the organized results tarball filename. | String | No | `cellranger_results` |
| `cellbender_gpu_enabled` | Enable GPU acceleration for CellBender. Set to `false` for CPU-only execution. | Boolean | No | `true` |
| `cellbender_expected_cells` | Expected number of real cells per sample passed to CellBender. | Int | No | auto |
| `cellbender_total_droplets_included` | Total number of droplets for CellBender to analyze per sample. | Int | No | auto |
| `cellbender_epochs` | Number of CellBender training epochs. | Int | No | `150` |
| `cellbender_low_count_threshold` | Droplets with total UMI count below this value are excluded from CellBender analysis. | Int | No | `5` |
| `cellbender_cpu_cores` | Number of CPU cores for CellBender. | Int | No | `4` |
| `cellbender_memory_gb` | Memory in GB for CellBender. | Int | No | `32` |

### Cell Ranger Reference

Cell Ranger requires a pre-built reference transcriptome tarball. You can:
- Download pre-built references from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads)
- Build custom references using `cellranger mkref`

## Output Files

| Output | Description | Source Module |
|--------|-------------|---------------|
| `single_cell_sample_list` | List of sample IDs that Cell Ranger ran successfully | pipeline |
| `skipped_sample_list` | List of sample IDs that Cell Ranger failed to auto-detect chemistry for (empty unless `skip_on_chemistry_failure = true`) | pipeline |
| `cellranger_results` | Compressed tarballs of Cell Ranger count output directories | ww-cellranger |
| `cellranger_web_summaries` | Web summary HTML files | ww-cellranger |
| `cellranger_metrics` | Metrics summary CSV files | ww-cellranger |
| `cellranger_filtered_h5s` | Filtered feature-barcode matrix HDF5 files | ww-cellranger |
| `cellranger_raw_h5s` | Raw feature-barcode matrix HDF5 files | ww-cellranger |
| `cellbender_output_h5s` | CellBender cleaned count matrices (all barcodes retained), one per successful sample | ww-cellbender |
| `cellbender_filtered_h5s` | CellBender filtered count matrices (barcodes with >50% cell probability), one per successful sample | ww-cellbender |
| `organized_results` | Tarball of all Cell Ranger and CellBender outputs organized into per-sample subdirectories (absent unless `organize_results = true`) | pipeline |

### Mixed single-cell / non-single-cell input

It's not uncommon for an SRA study to contain a mix of single-cell and bulk RNA-seq samples. When Cell Ranger can't assign a chemistry to a sample — typically because it isn't single-cell at all, but also when the reads are too few or too short for any chemistry's minimum requirements — it exits with an error. By default this fails the workflow for that scatter shard.

Setting `skip_on_chemistry_failure = true` (default `false`) makes such samples succeed with no count outputs instead. Their IDs are written to `skipped_sample_list`, the IDs of samples that completed normally go to `single_cell_sample_list`, and the `cellranger_*` output arrays contain results only for the successful samples. See [the module README](../../modules/ww-cellranger/README.md#graceful-chemistry-detection-skip) for detection details and the `chemistry` interaction. Note that very small `max_reads` values (a few hundred thousand or less) can cause legitimate single-cell samples to be skipped if the downsampled FASTQs fall below Cell Ranger's per-chemistry read-length minimums.

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
2. Downloads two SRA datasets, limited to 1M reads each: SRR7722937 (10x Chromium 3' v3 scRNA-seq demo dataset, expected to run successfully) and SRR1039508 (bulk RNA-seq, expected to trigger the `skip_on_chemistry_failure` path)
3. Renames FASTQs for Cell Ranger compatibility
4. Runs Cell Ranger count
5. Runs CellBender `remove-background` on the successful sample (CPU-only with reduced epochs and droplets for CI speed)
6. Validates that CellBender produced output for the single-cell sample and was correctly skipped for the bulk sample

### HPC Test Workflow

`testrun_hpc.wdl` mirrors the regular testrun's coverage but dispatches with `execution_mode = "hpc_sprocket"` so it exercises the Sprocket + module-load path end-to-end. It is intended to be exercised on Fred Hutch HPC via Sprocket on a monthly basis to validate the HPC dispatch path; the Cromwell-on-HPC path (`execution_mode = "hpc_cromwell"`) is exercised by users running Cromwell on their HPC by hand. Uses the same tiny SRA download and Cell Ranger reference as `testrun.wdl`.

## Related WILDS Components

- **ww-sra module**: SRA download functionality
- **ww-cellranger module**: Cell Ranger processing functionality
- **ww-cellbender module**: Ambient RNA removal functionality
- **ww-sra-star pipeline**: Alternative pipeline for bulk RNA-seq alignment
- **ww-sra-salmon pipeline**: Alternative pipeline for bulk RNA-seq quantification

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

Special thanks to [Hrishi Venkatesh](https://github.com/hvenkat94) for extensive testing of this pipeline on HPC infrastructure and for identifying important edge cases that have shaped its current scope and design. Thank you for your contributions!

If you would like to contribute to this WILDS WDL pipeline, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
