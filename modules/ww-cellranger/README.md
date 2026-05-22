# ww-cellranger

[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for single-cell RNA-seq analysis using Cell Ranger.

## Overview

This module provides reusable WDL tasks for preparing and processing single-cell RNA-seq data using Cell Ranger from 10x Genomics.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `run_count`, `run_count_hpc_cromwell`, `run_count_hpc_sprocket`, `rename_fastqs`
- **Test workflows**: `testrun.wdl` for CI/CD; `testrun_hpc.wdl` for monthly HPC validation of the Sprocket module-load path
- **Container**: Cell Ranger is not redistributable, so the WILDS Docker Library does not publish a public Cell Ranger image. A [Dockerfile recipe](https://github.com/getwilds/wilds-docker-library/blob/main/cellranger/Dockerfile_latest) is provided so users can build their own private image. Fred Hutch users running on institutional HPC can instead use `run_count_hpc_cromwell` or `run_count_hpc_sprocket` with the `CellRanger/10.0.0` environment module under Fred Hutch's institutional Cell Ranger license.

## Important Requirements

### Platform Requirements

**Cell Ranger requires Linux x86_64 architecture with AVX instruction support.** When using `run_count`, the binary runs inside a Docker container built from the [WILDS Cell Ranger Dockerfile recipe](https://github.com/getwilds/wilds-docker-library/blob/main/cellranger/Dockerfile_latest), which works on:
- Linux x86_64 systems (native)
- Intel-based Macs (via Docker)
- Cloud platforms and HPC clusters (GitHub Actions, AWS, Google Cloud, etc.)

When using either `run_count_hpc_cromwell` or `run_count_hpc_sprocket`, the same architecture requirement applies to the HPC compute nodes that load the Cell Ranger environment module.

**Note for Apple Silicon (M1/M2/M3) Mac users:** Local testing is not supported. Docker's x86_64 emulation on Apple Silicon does not support the AVX instructions that Cell Ranger requires. The CI/CD pipeline will test this module on compatible infrastructure.

### FASTQ Naming Convention

**Your FASTQ files must follow Cell Ranger's naming convention:**
- Format: `SampleName_S1_L001_R1_001.fastq.gz` (for Read 1)
- Format: `SampleName_S1_L001_R2_001.fastq.gz` (for Read 2)

If your files don't follow this convention, you can use the `rename_fastqs` task from this module to rename them automatically.

### Current Limitations

**This module currently only accepts samples that are gene expression (GEX) only.** Feature barcoding data (e.g., antibody-derived tags, CRISPR guides) is not supported at this time.

If you need support for feature barcoding or other Cell Ranger features, please file a feature request on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Tasks

### Picking a `run_count` variant

Cell Ranger's license prevents WILDS from publishing a public Docker image, and the two HPC engines this library supports (Cromwell and Sprocket) load host environment modules in different ways. This module therefore exposes three near-identical tasks; pick the one that matches your environment:

| Task | Use when | How Cell Ranger is obtained |
| --- | --- | --- |
| `run_count` | You have a private Cell Ranger Docker image (built locally from the [WILDS Dockerfile recipe](https://github.com/getwilds/wilds-docker-library/blob/main/cellranger/Dockerfile_latest) or supplied by your organization). Right choice for cloud, CI, and most local runs. | `docker` runtime key pulls a private Cell Ranger image |
| `run_count_hpc_cromwell` | You are running on an institutional HPC via Cromwell (e.g., Fred Hutch HPC via PROOF), and your Cromwell backend executes tasks directly on the compute node. | No `docker` runtime key; the task runs on the host and `module load` makes the licensed Cell Ranger binary available on `PATH` |
| `run_count_hpc_sprocket` | You are running on an institutional HPC via Sprocket, where tasks always run inside a container under Apptainer. | Runs inside a minimal Lua container (`getwilds/lua:5.3.6`); the host's Lmod tree, modulefiles, and Cell Ranger software tree are bind-mounted in via the Sprocket HPC config so `module load` works inside the container |

All three tasks produce identical outputs and accept the same scientific parameters; they differ only in how they obtain the Cell Ranger binary. The Sprocket variant requires a Sprocket configuration that bind-mounts the host's Lmod and Cell Ranger software trees into the container — see [`.github/configs/sprocket-hpc.toml`](../../.github/configs/sprocket-hpc.toml) for the reference setup used in the monthly HPC test run.

### `run_count`

Run `cellranger count` on gene expression reads from one GEM well using a private Cell Ranger Docker image.

**Inputs:**
- `r1_fastqs` (Array[File]): Array of R1 FASTQ files (contain cell barcodes and UMIs)
- `r2_fastqs` (Array[File]): Array of R2 FASTQ files (contain cDNA sequences)
- `ref_gex` (File): GEX reference transcriptome tarball (e.g., from 10x Genomics)
- `sample_id` (String): Sample ID for output naming
- `create_bam` (Boolean, default=true): Generate BAM file
- `cpu_cores` (Int, default=8): Number of CPU cores to use
- `memory_gb` (Int, default=64): Memory allocation in GB
- `expect_cells` (Int, optional): Expected number of recovered cells
- `chemistry` (String, optional): Assay configuration (e.g., SC3Pv3)
- `skip_on_chemistry_failure` (Boolean, default=false): If true, samples Cell Ranger can't assign a chemistry to (typically non-single-cell data) succeed with absent count outputs and `chemistry_status="skipped_non_single_cell"` instead of failing. Other Cell Ranger failures always re-raise. See [Graceful chemistry-detection skip](#graceful-chemistry-detection-skip) for details.
- `docker_image` (String, default=`ghcr.io/getwilds/cellranger:10.0.0`): Private Cell Ranger Docker image. Cell Ranger is not redistributable, so you must supply your own image (see the WILDS Dockerfile recipe linked above) and override the default if it is not available in your registry.

**Important:** All input FASTQs must be from one GEM well. If you have multiple GEM wells, use `run_count` separately for each well. The task validates that FASTQ filenames follow Cell Ranger's naming convention and will fail with a helpful error message if they don't.

**Outputs:**
- `chemistry_status` (File): One-line marker file with contents `ok` (Cell Ranger ran to completion) or `skipped_non_single_cell` (chemistry detection failed and `skip_on_chemistry_failure=true`). Always present.
- `results_tar` (File?): Compressed tarball of Cell Ranger count output directory. Absent when the sample was skipped.
- `web_summary` (File?): Web summary HTML file with QC metrics. Absent when the sample was skipped.
- `metrics_summary` (File?): Metrics summary CSV file with key statistics. Absent when the sample was skipped.
- `filtered_h5` (File?): Filtered feature-barcode matrix HDF5 file. Absent when the sample was skipped.

### Graceful chemistry-detection skip

When `cellranger count` cannot determine which single-cell chemistry an input came from (a common symptom of non-single-cell data being fed in by mistake — e.g. bulk RNA-seq mixed into an SRA study), it exits with an error before producing any of its usual outputs. By default, the task re-raises that failure. Setting `skip_on_chemistry_failure = true` makes the task exit cleanly with `chemistry_status = "skipped_non_single_cell"` and no count outputs, so a scatter over many samples can complete even when a subset turn out not to be single-cell.

The detection is heuristic — it greps Cell Ranger's stderr for substrings like "could not detect", "ambiguous chemistry", and `NO_INPUT_ANTIBODY_READS`. Any failure that doesn't match these markers (OOM, corrupt reference, malformed FASTQ, etc.) still re-raises so real bugs aren't silently swallowed. The marker list may need tuning across Cell Ranger versions; this module has been validated against Cell Ranger 10.0.0.

> **Important:** `skip_on_chemistry_failure` only triggers when Cell Ranger is allowed to auto-detect the chemistry. If you also pass `chemistry` (e.g. `chemistry = "SC3Pv3"`), Cell Ranger skips detection entirely and a non-single-cell sample will fail with a different error (typically barcode-validation related) that this heuristic does not catch. Leave `chemistry` unset when you want the graceful-skip path to be available.

### `run_count_hpc_cromwell`

Run `cellranger count` on gene expression reads from one GEM well using the host's environment-module system instead of a Docker image. The task omits the `docker` runtime key entirely so that Cromwell's HPC backend executes it directly on the compute node, where `module load` makes the licensed Cell Ranger binary available on `PATH`. Intended for Cromwell-on-HPC deployments such as the Fred Hutch HPC via PROOF; on container-based backends, no Cell Ranger binary will be on `PATH` and the task will fail.

**Inputs:** Same as `run_count`, except `docker_image` is replaced by:
- `cellranger_module` (String, default=`"CellRanger/10.0.0"`): HPC environment module to load before invoking `cellranger count`.

**Outputs:** Same as `run_count`.

Example, Cromwell-on-HPC with environment modules:
```wdl
call cellranger_tasks.run_count_hpc_cromwell {
  input:
    r1_fastqs = r1_fastqs,
    r2_fastqs = r2_fastqs,
    ref_gex = reference,
    sample_id = "my_sample",
    cellranger_module = "CellRanger/10.0.0"
}
```

### `run_count_hpc_sprocket`

Run `cellranger count` on gene expression reads from one GEM well from inside a minimal Lua container, with the host's Lmod tree and Cell Ranger software tree bind-mounted in. Sprocket always runs tasks inside a container under Apptainer, so this variant cannot omit the `docker` runtime key the way `run_count_hpc_cromwell` does; instead the container ships only the bits needed to execute `module load` (Lua) and the Cell Ranger binary itself comes in via host bind-mounts configured in the Sprocket HPC config.

**Inputs:** Same as `run_count_hpc_cromwell`.

**Outputs:** Same as `run_count`.

Example, Sprocket-on-HPC with environment modules:
```wdl
call cellranger_tasks.run_count_hpc_sprocket {
  input:
    r1_fastqs = r1_fastqs,
    r2_fastqs = r2_fastqs,
    ref_gex = reference,
    sample_id = "my_sample",
    cellranger_module = "CellRanger/10.0.0"
}
```

> **Sprocket configuration note:** This task assumes a Sprocket config that bind-mounts `/app/lmod`, the host's modulefiles, and the host's Cell Ranger software tree into the container. The reference setup used by the monthly HPC test run lives in [`.github/configs/sprocket-hpc.toml`](../../.github/configs/sprocket-hpc.toml).

### `rename_fastqs`

Renames FASTQ files to match the Cell Ranger naming convention. Useful for preparing SRA downloads or other non-standard FASTQ files for Cell Ranger input.

**Inputs:**
- `r1_fastq` (File): R1 FASTQ file to rename
- `r2_fastq` (File): R2 FASTQ file to rename
- `sample_id` (String): Sample ID to use as the prefix in the renamed file

**Outputs:**
- `r1_renamed` (File): R1 FASTQ file renamed to `{sample_id}_S1_R1_001.fastq.gz`
- `r2_renamed` (File): R2 FASTQ file renamed to `{sample_id}_S1_R2_001.fastq.gz`

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellranger/ww-cellranger.wdl" as cellranger_tasks

workflow my_single_cell_pipeline {
  input {
    Array[File] r1_fastqs  # Must follow Cell Ranger naming convention
    Array[File] r2_fastqs  # Must follow Cell Ranger naming convention
    String sample_name
    File gex_reference
  }

  # Run Cell Ranger count
  call cellranger_tasks.run_count {
    input:
      r1_fastqs = r1_fastqs,
      r2_fastqs = r2_fastqs,
      ref_gex = gex_reference,
      sample_id = sample_name
  }

  output {
    File results = run_count.results_tar
    File web_summary = run_count.web_summary
    File metrics = run_count.metrics_summary
    File filtered_h5 = run_count.filtered_h5
  }
}
```

### Advanced Usage Examples

**Custom resource allocation:**
```wdl
call cellranger_tasks.run_count {
  input:
    r1_fastqs = r1_fastqs,
    r2_fastqs = r2_fastqs,
    ref_gex = reference,
    sample_id = "my_sample",
    cpu_cores = 16,
    memory_gb = 128
}
```

**Specifying chemistry and expected cells:**
```wdl
call cellranger_tasks.run_count {
  input:
    r1_fastqs = r1_fastqs,
    r2_fastqs = r2_fastqs,
    ref_gex = reference,
    sample_id = "my_sample",
    chemistry = "SC3Pv3",
    expect_cells = 5000
}
```

**Without BAM generation (faster, less storage):**
```wdl
call cellranger_tasks.run_count {
  input:
    r1_fastqs = r1_fastqs,
    r2_fastqs = r2_fastqs,
    ref_gex = reference,
    sample_id = "my_sample",
    create_bam = false
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-sra**: Download public single-cell sequencing data before analysis

## Testing the Module

The module includes a demonstration workflow that can be tested independently. The workflow in `testrun.wdl` automatically downloads test data and runs without requiring input files.

**Fred Hutch users**: If you would like run `testrun.wdl` on PROOF, first edit it to use the `run_count_hpc` task instead of `run_count`. That way the Cell Ranger module will be used instead of the private Docker image used by CI/CD.

The test workflow (`cellranger_example`) automatically:
1. Downloads a small GEX reference using `ww-testdata`
2. Downloads test FASTQ data with proper naming convention using `ww-testdata`
3. Runs Cell Ranger count analysis
4. Validates all outputs

### HPC Test Workflow

`testrun_hpc.wdl` mirrors the regular testrun's coverage but dispatches to `run_count_hpc_sprocket` instead of `run_count`. It is intended to be exercised on Fred Hutch HPC via Sprocket on a monthly basis to validate the module-load path; the Cromwell variant (`run_count_hpc_cromwell`) is exercised by users running Cromwell on their HPC by hand. Uses the same tiny test data as `testrun.wdl` so the only thing this run validates is the HPC dispatch path.

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 16-128 GB recommended (scales with dataset size)
  - Small datasets (<5000 cells): 16-32 GB
  - Medium datasets (5000-10000 cells): 32-64 GB
  - Large datasets (>10000 cells): 64-128 GB
- **CPUs**: 4-16 cores recommended; Cell Ranger benefits from multi-threading
- **Disk space**: Ensure sufficient space for output files
  - With BAM: ~2-5x the size of input FASTQs
  - Without BAM: ~0.5-1x the size of input FASTQs

### Advanced Parameters

- `chemistry`: Specify the assay configuration if auto-detection fails (e.g., "SC3Pv2", "SC3Pv3")
- `expect_cells`: Provide expected cell count if desired
- `create_bam`: Set to `false` to save storage and processing time if BAM files are truly not needed.

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- One of: Docker/Apptainer support with a private Cell Ranger image (for `run_count`), a Cromwell-on-HPC backend that runs tasks directly on the compute node and provides Cell Ranger via Lmod (for `run_count_hpc_cromwell`), or a Sprocket-on-HPC config that bind-mounts the host's Lmod and Cell Ranger software trees into the container (for `run_count_hpc_sprocket`)
- Sufficient computational resources (Cell Ranger can be memory-intensive)
- 10x Genomics reference transcriptome (available from [10x Genomics Download Center](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads))

## Performance Considerations

- **Memory usage**: Cell Ranger requires substantial RAM; allocate 64GB+ for typical datasets
- **CPU scaling**: Cell Ranger benefits significantly from multiple cores (8-16 recommended)
- **Storage requirements**:
  - Ensure sufficient space for outputs (2-5x input size with BAM)
  - Consider setting `create_bam = false` if BAM files are truly not needed
- **Runtime**: Typical runtime ranges greatly from 1-12 hours depending on dataset size and resources

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Test data from the ww-testdata module
- Comprehensive validation of all outputs

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to Cell Ranger usage or configuration, please refer to the [Cell Ranger documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
