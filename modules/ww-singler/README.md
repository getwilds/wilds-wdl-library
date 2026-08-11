# ww-singler
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for single-cell RNA-seq clustering, marker gene identification, and cell-type annotation using scran and SingleR.

## Overview

This module provides a reusable WDL task for graph-based clustering, per-cluster marker gene identification, and reference-based cell-type annotation of single-cell RNA-seq data using [scran](https://bioconductor.org/packages/release/bioc/html/scran.html) (`clusterCells`, `findMarkers`) and [SingleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html), with reference datasets from [celldex](https://bioconductor.org/packages/release/bioc/html/celldex.html). It accepts a dimensionality-reduced `SingleCellExperiment` RDS object with a PCA embedding (e.g. produced by [`ww-scater`](../ww-scater/)) and runs the following steps:

1. **Load data** — reads the `SingleCellExperiment` RDS object and checks for a `logcounts` assay and a `PCA` reducedDim
2. **Clustering** — assigns cells to clusters with `scran::clusterCells` using the PCA embedding
3. **Marker gene identification** — finds per-cluster marker genes with `scran::findMarkers`
4. **Cell-type annotation** — annotates each cluster with `SingleR`, using either a named [celldex](https://bioconductor.org/packages/release/bioc/html/celldex.html) reference dataset (default) or a user-supplied reference `SingleCellExperiment`, and saves the updated object

This module is intended to run downstream of [`ww-scater`](../ww-scater/) in the standard single-cell post-processing chain (load matrix -> QC filter -> normalize -> HVG -> PCA -> UMAP -> clustering -> marker genes -> cell-type annotation).

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `run_singler`
- **Scripts**: `singler_analysis.R` (fetched via curl at runtime)
- **Test workflow**: `testrun.wdl`
- **Container**: `getwilds/singler:2.14.1`

## Scripts

All scripts are fetched from this repository at runtime via `curl`.

| Script | Used by | Language | Description |
|--------|---------|----------|-------------|
| [`singler_analysis.R`](singler_analysis.R) | `run_singler` | R | Loads a dimensionality-reduced `SingleCellExperiment` RDS object, clusters cells, identifies per-cluster marker genes, and annotates cell types with SingleR |

## Tasks

### `run_singler`

Performs clustering, marker gene identification, and cell-type annotation on a dimensionality-reduced `SingleCellExperiment` RDS object.

**Inputs:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `sce_rds` | File | — | Dimensionality-reduced `SingleCellExperiment` RDS object (e.g. from `ww-scater`) containing a PCA reducedDim |
| `sample_name` | String | — | Sample name used for output file prefixes |
| `reference_rds` | File? | — | Optional custom reference `SingleCellExperiment` RDS object with labeled cell types. If omitted, a celldex reference dataset is fetched by name |
| `reference_dataset` | String | `HumanPrimaryCellAtlasData` | Name of the celldex reference dataset function to fetch (e.g. `HumanPrimaryCellAtlasData`, `MouseRNAseqData`), used when `reference_rds` is not provided |
| `label_column` | String | `label.main` | Column name in the reference object's colData containing cell-type labels |
| `n_top_markers` | Int | 10 | Number of top marker genes to report per cluster |
| `random_seed` | Int | 100 | Random seed for clustering, for reproducibility |
| `memory_gb` | Int | 16 | Memory allocated for the task in GB |
| `cpu_cores` | Int | 2 | Number of CPU cores allocated for the task |
| `docker_image` | String | `getwilds/singler:2.14.1` | Docker image to use for this task |

**Outputs:**

| Name | Type | Description |
|------|------|-------------|
| `sce_object` | File | `SingleCellExperiment` RDS object with cluster assignments and predicted cell types |
| `cluster_table` | File | CSV mapping each cell barcode to its assigned cluster |
| `marker_table` | File | CSV of top marker genes per cluster |
| `prediction_table` | File | CSV of SingleR cell-type predictions and scores per cluster |

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-singler/ww-singler.wdl" as singler_tasks

workflow my_scrna_analysis {
  input {
    File sce_rds
    String sample_name
  }

  call singler_tasks.run_singler {
    input:
      sce_rds     = sce_rds,
      sample_name = sample_name
  }

  output {
    File sce_object      = run_singler.sce_object
    File cluster_table    = run_singler.cluster_table
    File marker_table     = run_singler.marker_table
    File prediction_table = run_singler.prediction_table
  }
}
```

### Advanced Usage Examples

**Mouse reference dataset:**
```wdl
call singler_tasks.run_singler {
  input:
    sce_rds            = sce_rds,
    sample_name        = "mouse_sample",
    reference_dataset  = "MouseRNAseqData"
}
```

**Custom reference SingleCellExperiment:**
```wdl
call singler_tasks.run_singler {
  input:
    sce_rds        = sce_rds,
    sample_name    = "sample",
    reference_rds  = custom_reference_sce,
    label_column   = "cell_type"
}
```

## Testing the Module

The test workflow (`testrun.wdl`) automatically downloads a public 10x Genomics PBMC dataset via `ww-testdata`, loads it into a `SingleCellExperiment` via `ww-dropletutils`, normalizes it with `ww-scran`, runs dimensionality reduction with `ww-scater`, and runs `run_singler` on it, then validates all output files. Note that the demo dataset is rat PBMCs annotated against the default human reference (`HumanPrimaryCellAtlasData`), so predicted labels are not biologically meaningful for this test data — the test only validates that the pipeline runs end-to-end and produces well-formed outputs.

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint singler_example
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Internet access for fetching scripts from GitHub and reference datasets from celldex at runtime
- R environment with scran, SingleR, and celldex (provided by `getwilds/singler:2.14.1` container)

## Citation

> Aran D, Looney AP, Liu L, Wu E, Fong V, Hsu A, Chak S, et al. (2019). "Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage." Nature Immunology, 20(2), 163-172.
> Lun ATL, McCarthy DJ, Marioni JC (2016). "A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor." F1000Research, 5, 2122.

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
