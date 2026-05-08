# ww-seurat
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for single-cell RNA-seq analysis using Seurat.

## Overview

This module provides a reusable WDL task for QC, normalization, clustering, and marker gene identification of single-cell RNA-seq data using [Seurat](https://satijalab.org/seurat/). It accepts a Cell Ranger HDF5 output file and runs the following steps:

1. **Load data** — reads the filtered feature-barcode matrix from a `.h5` file
2. **QC filtering** — calculates mitochondrial gene percentage and filters cells by minimum feature count and maximum mitochondrial percentage
3. **Normalization** — runs SCTransform (regressing out mitochondrial percentage) followed by PCA and UMAP (top 30 PCs)
4. **Clustering** — Louvain clustering at a user-specified resolution
5. **Marker genes** — runs `FindAllMarkers` (positive markers only) and saves the top 30 markers per cluster; generates a heatmap of the top 8 markers per cluster

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `run_seurat`
- **Scripts**: `seurat_analysis.R` (fetched via curl at runtime)
- **Test workflow**: `testrun.wdl`
- **Container**: `getwilds/seurat:5.2.1`

## Scripts

All scripts are fetched from this repository at runtime via `curl`.

| Script | Used by | Language | Description |
|--------|---------|----------|-------------|
| [`seurat_analysis.R`](seurat_analysis.R) | `run_seurat` | R | Loads a Cell Ranger `.h5` matrix, filters cells by QC thresholds, normalizes with SCTransform, runs PCA/UMAP, performs Louvain clustering, and identifies marker genes per cluster |

## Tasks

### `run_seurat`

Performs QC, normalization, clustering, and marker gene identification on a Cell Ranger HDF5 matrix file.

**Inputs:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `input_h5` | File | — | Cell Ranger filtered feature-barcode matrix HDF5 file (`filtered_feature_bc_matrix.h5`) |
| `sample_name` | String | — | Sample name used for output file prefixes |
| `min_cells` | Int | 3 | Minimum number of cells a gene must be detected in to be retained |
| `min_features` | Int | 200 | Minimum number of features (genes) a cell must have to pass QC |
| `max_percent_mt` | Float | 10.0 | Maximum mitochondrial gene percentage allowed per cell |
| `resolution` | Float | 0.5 | Louvain clustering resolution (higher values produce more clusters) |
| `memory_gb` | Int | 16 | Memory allocated for the task in GB |
| `cpu_cores` | Int | 4 | Number of CPU cores allocated for the task |

**Outputs:**

| Name | Type | Description |
|------|------|-------------|
| `seurat_object` | File | Processed Seurat RDS object |
| `qc_plot` | File | QC violin plot showing nFeature_RNA, nCount_RNA, and percent.mt (pre-filtering) |
| `umap_plot` | File | UMAP plot colored by Louvain cluster |
| `heatmap_plot` | File | Heatmap of top 8 marker genes per cluster |
| `cluster_markers` | File | CSV of top 30 marker genes per cluster from `FindAllMarkers` |

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-seurat/ww-seurat.wdl" as seurat_tasks

workflow my_scrna_analysis {
  input {
    File cellranger_h5
    String sample_name
  }

  call seurat_tasks.run_seurat {
    input:
      input_h5    = cellranger_h5,
      sample_name = sample_name
  }

  output {
    File seurat_object   = run_seurat.seurat_object
    File qc_plot         = run_seurat.qc_plot
    File umap_plot       = run_seurat.umap_plot
    File heatmap_plot    = run_seurat.heatmap_plot
    File cluster_markers = run_seurat.cluster_markers
  }
}
```

## Testing the Module

The test workflow (`testrun.wdl`) automatically downloads a public 10x Genomics PBMC dataset via `ww-testdata` and runs `run_seurat` on it, then validates all output files.

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint seurat_example
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Internet access for fetching scripts from GitHub at runtime
- R environment with Seurat and supporting packages (provided by `getwilds/seurat:5.2.1` container)

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
