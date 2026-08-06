# ww-scran
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for single-cell RNA-seq QC, normalization, and feature selection using scran.

## Overview

This module provides a reusable WDL task for per-cell QC filtering, deconvolution-based normalization, and highly variable gene (HVG) selection of single-cell RNA-seq data using [scran](https://bioconductor.org/packages/release/bioc/html/scran.html) (and its hard dependencies, e.g. `scuttle`, only — no `scater`). It accepts a `SingleCellExperiment` RDS object (e.g. produced by [`ww-dropletutils`](../ww-dropletutils/)) and runs the following steps:

1. **Load data** — reads the `SingleCellExperiment` RDS object
2. **QC filtering** — computes per-cell QC metrics (library size, detected genes, mitochondrial percentage) directly from the counts matrix and flags outliers using MAD-based thresholds
3. **Normalization** — estimates size factors via scran's pooling/deconvolution method (`quickCluster` + `computeSumFactors`) and applies `logNormCounts`
4. **Feature selection** — models per-gene mean-variance trends with `modelGeneVar` and selects the top highly variable genes
5. **Dimensionality reduction** — runs PCA (base R `prcomp`) on the selected HVGs and saves the resulting `SingleCellExperiment` object

This module is intended to run downstream of [`ww-dropletutils`](../ww-dropletutils/) in the standard single-cell post-processing chain (load matrix -> QC filter -> normalize -> HVG -> PCA -> ...).

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `run_scran`
- **Scripts**: `scran_analysis.R` (fetched via curl at runtime)
- **Test workflow**: `testrun.wdl`
- **Container**: `getwilds/scran:1.40.0`

## Scripts

All scripts are fetched from this repository at runtime via `curl`.

| Script | Used by | Language | Description |
|--------|---------|----------|-------------|
| [`scran_analysis.R`](scran_analysis.R) | `run_scran` | R | Loads a `SingleCellExperiment` RDS object, filters cells by QC thresholds, normalizes with scran's deconvolution method, models mean-variance trends, and runs PCA on highly variable genes |

## Tasks

### `run_scran`

Performs QC filtering, normalization, and highly variable gene selection on a `SingleCellExperiment` RDS object.

**Inputs:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `sce_rds` | File | — | `SingleCellExperiment` RDS object (e.g. from `ww-dropletutils`) containing the raw counts to QC and normalize |
| `sample_name` | String | — | Sample name used for output file prefixes |
| `mito_pattern` | String | `^MT-` | Regex pattern identifying mitochondrial gene symbols |
| `nmads` | Float | 3.0 | Number of median absolute deviations from the median used to flag low-quality cells |
| `min_mean` | Float | 0.1 | Minimum mean expression for a gene to be used in size factor estimation and HVG modeling |
| `n_hvgs` | Int | 2000 | Number of top highly variable genes to select for PCA |
| `memory_gb` | Int | 8 | Memory allocated for the task in GB |
| `cpu_cores` | Int | 2 | Number of CPU cores allocated for the task |
| `docker_image` | String | `getwilds/scran:1.40.0` | Docker image to use for this task |

**Outputs:**

| Name | Type | Description |
|------|------|-------------|
| `sce_object` | File | Normalized `SingleCellExperiment` RDS object with size factors and PCA |
| `qc_plot` | File | Per-cell QC scatter plot (library size vs. detected genes, colored by discard status) |
| `size_factor_plot` | File | Histogram of scran deconvolution size factors |
| `mean_variance_plot` | File | Mean-variance trend plot from highly variable gene modeling |
| `hvg_table` | File | CSV of all genes ranked by biological variance component |

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-scran/ww-scran.wdl" as scran_tasks

workflow my_scrna_analysis {
  input {
    File sce_rds
    String sample_name
  }

  call scran_tasks.run_scran {
    input:
      sce_rds     = sce_rds,
      sample_name = sample_name
  }

  output {
    File sce_object          = run_scran.sce_object
    File qc_plot              = run_scran.qc_plot
    File size_factor_plot     = run_scran.size_factor_plot
    File mean_variance_plot   = run_scran.mean_variance_plot
    File hvg_table            = run_scran.hvg_table
  }
}
```

### Advanced Usage Examples

**Non-human data (mouse mitochondrial gene prefix):**
```wdl
call scran_tasks.run_scran {
  input:
    sce_rds      = sce_rds,
    sample_name  = "mouse_sample",
    mito_pattern = "^mt-"
}
```

**Stricter QC and larger HVG set:**
```wdl
call scran_tasks.run_scran {
  input:
    sce_rds     = sce_rds,
    sample_name = "sample",
    nmads       = 2.5,
    n_hvgs      = 5000
}
```

## Testing the Module

The test workflow (`testrun.wdl`) automatically downloads a public 10x Genomics PBMC dataset via `ww-testdata`, loads it into a `SingleCellExperiment` via `ww-dropletutils`, and runs `run_scran` on it, then validates all output files.

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint scran_example
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Internet access for fetching scripts from GitHub at runtime
- R environment with scran (provided by `getwilds/scran:1.40.0` container)

## Citation

> Lun ATL, McCarthy DJ, Marioni JC (2016). "A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor." F1000Research, 5, 2122.
> Lun ATL, Bach K, Marioni JC (2016). "Pooling across cells to normalize single-cell RNA sequencing data with many zero counts." Genome Biology, 17, 75.

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
