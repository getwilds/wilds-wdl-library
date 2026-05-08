# ww-seurat
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for single-cell RNA-seq analysis using Seurat.

## Overview

TODO: Add overview description.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `run_seurat`
- **Scripts**: `seurat_analysis.R` (fetched via curl at runtime)
- **Test workflow**: `testrun.wdl`
- **Containers**: `getwilds/seurat:TODO`

## Scripts

All scripts are fetched from this repository at runtime via `curl`.

| Script | Used by | Language | Description |
|--------|---------|----------|-------------|
| [`seurat_analysis.R`](seurat_analysis.R) | `run_seurat` | R | TODO |

## Tasks

### `run_seurat`

Performs single-cell RNA-seq analysis using Seurat.

**Inputs:**
- `input_matrix_dir` (Directory): Cell Ranger output directory (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
- `sample_name` (String): Sample name used for output file prefixes
- `min_cells` (Int): Minimum number of cells a gene must be detected in (default: 3)
- `min_features` (Int): Minimum number of features per cell (default: 200)
- `max_features` (Int): Maximum number of features per cell (default: 2500)
- `max_percent_mt` (Float): Maximum mitochondrial gene percentage per cell (default: 5.0)
- `n_variable_features` (Int): Number of highly variable features for PCA (default: 2000)
- `n_dims` (Int): Number of PCA dimensions for clustering and UMAP (default: 10)
- `resolution` (Float): Clustering resolution (default: 0.5)
- `memory_gb` (Int): Memory allocation in GB (default: 16)
- `cpu_cores` (Int): Number of CPU cores (default: 4)

**Outputs:**
- `seurat_object` (File): Processed Seurat RDS object
- `umap_plot` (File): UMAP plot showing cell clusters
- `cluster_markers` (File): Marker genes identified for each cluster

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-seurat/ww-seurat.wdl" as seurat_tasks

workflow my_scrna_analysis {
  input {
    Directory cellranger_output
    String sample_name
  }

  call seurat_tasks.run_seurat {
    input:
      input_matrix_dir = cellranger_output,
      sample_name = sample_name
  }

  output {
    File seurat_object   = run_seurat.seurat_object
    File umap_plot       = run_seurat.umap_plot
    File cluster_markers = run_seurat.cluster_markers
  }
}
```

## Testing the Module

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl --inputs testrun_inputs.json

# Using miniWDL
miniwdl run testrun.wdl input_matrix_dir=<path>

# Using Sprocket
sprocket run testrun.wdl --entrypoint seurat_example
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Internet access for fetching scripts from GitHub at runtime
- R environment with Seurat package (provided by `getwilds/seurat:TODO` container)

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
