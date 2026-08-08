# ww-scater
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for single-cell RNA-seq dimensionality reduction and QC visualization using scater.

## Overview

This module provides a reusable WDL task for PCA and UMAP dimensionality reduction, plus QC and reduced-dimension diagnostic plotting, of single-cell RNA-seq data using [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) only (no `scran`). It accepts a normalized `SingleCellExperiment` RDS object with a `logcounts` assay (e.g. produced by [`ww-scran`](../ww-scran/)) and runs the following steps:

1. **Load data** — reads the `SingleCellExperiment` RDS object and checks for a `logcounts` assay
2. **Feature selection** — ranks genes by variance of log-expression and selects the top highly variable genes for dimensionality reduction
3. **Dimensionality reduction** — runs PCA (`runPCA`) on the selected genes, then UMAP (`runUMAP`) on the resulting PCA embedding
4. **QC metrics** — computes per-cell QC metrics with `perCellQCMetrics`
5. **Visualization** — generates PCA, UMAP, and per-cell QC scatter plots with `plotReducedDim`/`plotColData`, and saves the updated `SingleCellExperiment` object

This module is intended to run downstream of [`ww-scran`](../ww-scran/) in the standard single-cell post-processing chain (load matrix -> QC filter -> normalize -> HVG -> PCA -> UMAP -> ...).

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `run_scater`
- **Scripts**: `scater_analysis.R` (fetched via curl at runtime)
- **Test workflow**: `testrun.wdl`
- **Container**: `getwilds/scater:1.40.2`

## Scripts

All scripts are fetched from this repository at runtime via `curl`.

| Script | Used by | Language | Description |
|--------|---------|----------|-------------|
| [`scater_analysis.R`](scater_analysis.R) | `run_scater` | R | Loads a normalized `SingleCellExperiment` RDS object, selects highly variable genes, runs PCA and UMAP, computes per-cell QC metrics, and generates diagnostic plots |

## Tasks

### `run_scater`

Performs PCA and UMAP dimensionality reduction and QC/reduced-dimension visualization on a normalized `SingleCellExperiment` RDS object.

**Inputs:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `sce_rds` | File | — | Normalized `SingleCellExperiment` RDS object (e.g. from `ww-scran`) containing a `logcounts` assay |
| `sample_name` | String | — | Sample name used for output file prefixes |
| `n_pcs` | Int | 50 | Number of principal components to compute |
| `n_hvgs` | Int | 2000 | Number of most variable genes (by log-expression variance) to use for PCA |
| `random_seed` | Int | 100 | Random seed for UMAP, for reproducibility |
| `memory_gb` | Int | 8 | Memory allocated for the task in GB |
| `cpu_cores` | Int | 2 | Number of CPU cores allocated for the task |
| `docker_image` | String | `getwilds/scater:1.40.2` | Docker image to use for this task |

**Outputs:**

| Name | Type | Description |
|------|------|-------------|
| `sce_object` | File | `SingleCellExperiment` RDS object with PCA and UMAP embeddings and per-cell QC metrics |
| `pca_plot` | File | Scatter plot of the first two PCA components |
| `umap_plot` | File | Scatter plot of the UMAP embedding |
| `qc_plot` | File | Per-cell QC scatter plot (library size vs. detected genes, colored by total counts) |

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-scater/ww-scater.wdl" as scater_tasks

workflow my_scrna_analysis {
  input {
    File sce_rds
    String sample_name
  }

  call scater_tasks.run_scater {
    input:
      sce_rds     = sce_rds,
      sample_name = sample_name
  }

  output {
    File sce_object = run_scater.sce_object
    File pca_plot    = run_scater.pca_plot
    File umap_plot   = run_scater.umap_plot
    File qc_plot     = run_scater.qc_plot
  }
}
```

### Advanced Usage Examples

**More principal components and a larger HVG set:**
```wdl
call scater_tasks.run_scater {
  input:
    sce_rds     = sce_rds,
    sample_name = "sample",
    n_pcs       = 30,
    n_hvgs      = 5000
}
```

## Testing the Module

The test workflow (`testrun.wdl`) automatically downloads a public 10x Genomics PBMC dataset via `ww-testdata`, loads it into a `SingleCellExperiment` via `ww-dropletutils`, normalizes it with `ww-scran`, and runs `run_scater` on it, then validates all output files.

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint scater_example
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Internet access for fetching scripts from GitHub at runtime
- R environment with scater (provided by `getwilds/scater:1.40.2` container)

## Citation

> McCarthy DJ, Campbell KR, Lun ATL, Wills QF (2017). "Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R." Bioinformatics, 33(8), 1179-1186.

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
