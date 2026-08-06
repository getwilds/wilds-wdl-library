# ww-dropletutils Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for loading 10x Genomics feature-barcode matrices and calling real cells from empty droplets using the Bioconductor [DropletUtils](https://bioconductor.org/packages/DropletUtils/) package. This is the first module in the standard single-cell post-processing chain (load matrix -> QC filter -> normalize -> HVG -> PCA -> neighbors -> UMAP -> clustering -> marker genes), producing a `SingleCellExperiment` object ready for downstream Bioconductor modules such as `ww-scran`, `ww-scater`, and `ww-singler`.

## Overview

Cell Ranger's `count` pipeline outputs raw and filtered feature-barcode matrices, but distinguishing real cells from empty droplets containing only ambient RNA is a standard first QC step before any downstream analysis. This module wraps two DropletUtils operations:

- Loading a 10x HDF5 matrix directly into a `SingleCellExperiment` object (for already-filtered matrices)
- Calling real cells from a raw, unfiltered matrix using the `emptyDrops` statistical test, along with a barcode-rank QC plot

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-dropletutils.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `read10x_counts`

Loads a 10x Genomics feature-barcode matrix (H5 format) into a `SingleCellExperiment` object using `DropletUtils::read10xCounts`. Use this for matrices that are already filtered (e.g. Cell Ranger's `filtered_feature_bc_matrix.h5`).

**Inputs:**
- `h5_matrix` (File): 10x Genomics feature-barcode matrix in HDF5 format (filtered or raw)
- `sample_name` (String): Sample name used as output file prefix and SCE column metadata
- `cpu_cores` (Int, default=2): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB
- `docker_image` (String, default=`getwilds/dropletutils:1.32.0`): Docker image to use for this task

**Outputs:**
- `sce_rds` (File): `SingleCellExperiment` object (RDS) containing the loaded counts, gene, and barcode metadata

### `empty_drops_filter`

Calls real cells from a raw (unfiltered) 10x feature-barcode matrix using `DropletUtils::emptyDrops`, filters the `SingleCellExperiment` to called cells, and produces a barcode-rank QC plot. Use this for raw matrices (e.g. Cell Ranger's `raw_feature_bc_matrix.h5`) that include all detected barcodes, including empty droplets.

**Inputs:**
- `raw_h5_matrix` (File): Raw (unfiltered) 10x Genomics feature-barcode matrix in HDF5 format, containing all detected barcodes
- `sample_name` (String): Sample name used as output file prefix and SCE column metadata
- `fdr_threshold` (Float, default=0.01): False discovery rate threshold below which a barcode is called a real cell
- `lower_umi_threshold` (Int, default=500): UMI count threshold below which barcodes are assumed to be empty droplets and used to estimate the ambient RNA profile. Also controls how many barcodes above this threshold get tested by `emptyDrops`'s Monte Carlo procedure; raising it reduces runtime on matrices with a large raw barcode population at the cost of testing fewer ambiguous barcodes
- `random_seed` (Int, default=100): Random seed for emptyDrops Monte Carlo p-value computation, for reproducibility
- `cpu_cores` (Int, default=2): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB
- `docker_image` (String, default=`getwilds/dropletutils:1.32.0`): Docker image to use for this task

**Outputs:**
- `filtered_sce_rds` (File): `SingleCellExperiment` object (RDS) filtered to barcodes called as cells by `emptyDrops`
- `empty_drops_csv` (File): Per-barcode `emptyDrops` statistics (Total, LogProb, PValue, FDR) for all tested barcodes
- `barcode_rank_pdf` (File): Barcode-rank plot (log-log UMI count vs. rank) with knee and inflection points marked

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-dropletutils/ww-dropletutils.wdl" as dropletutils_tasks

workflow my_analysis_pipeline {
  input {
    File raw_h5_matrix
    String sample_name
  }

  call dropletutils_tasks.empty_drops_filter {
    input:
      raw_h5_matrix = raw_h5_matrix,
      sample_name = sample_name
  }

  output {
    File filtered_sce_rds = empty_drops_filter.filtered_sce_rds
  }
}
```

### Advanced Usage Examples

**Custom resource allocation and thresholds:**
```wdl
call dropletutils_tasks.empty_drops_filter {
  input:
    raw_h5_matrix = large_raw_matrix,
    sample_name = "large_sample",
    fdr_threshold = 0.001,
    lower_umi_threshold = 200,
    cpu_cores = 4,
    memory_gb = 16
}
```

**Loading an already-filtered matrix (skip emptyDrops):**
```wdl
call dropletutils_tasks.read10x_counts {
  input:
    h5_matrix = filtered_matrix,
    sample_name = "prefiltered_sample"
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of 10x filtered and raw H5 test matrices for demonstrations
- **ww-cellbender**: Can be used upstream of or alongside `ww-cellbender` for ambient RNA removal; `emptyDrops` filtering and CellBender's `remove-background` both accept the same raw feature-barcode matrix
- **Downstream Bioconductor modules**: The `SingleCellExperiment` RDS output is intended as input for `ww-scran` (normalization, HVG selection, clustering), `ww-scater` (QC metrics, PCA/UMAP), and `ww-singler` (cell-type annotation)

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads a filtered 10x H5 matrix and a raw 10x H5 matrix using `ww-testdata`
2. Loads the filtered matrix directly into a `SingleCellExperiment` with `read10x_counts`
3. Calls real cells from the raw matrix with `empty_drops_filter`, demonstrating both module tasks in a realistic workflow context

## Docker Container

This module uses the `getwilds/dropletutils:1.32.0` container image, which includes:
- Bioconductor DropletUtils (version 1.32.0)
- R and required Bioconductor/CRAN dependencies (SingleCellExperiment, Matrix, HDF5Array, etc.)
- All necessary system dependencies for reading 10x HDF5 matrices

## Citation

> Lun, A.T.L., Riesenfeld, S., Andrews, T. et al.
> EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data.
> Genome Biology 20, 63 (2019).
> DOI: [10.1186/s13059-019-1662-y](https://doi.org/10.1186/s13059-019-1662-y)

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 2 cores
- **Memory**: 8 GB
- **Runtime**: A few minutes per sample for demo data; scales with matrix size and number of barcodes tested by `emptyDrops`

### Resource Scaling
- `cpu_cores`: DropletUtils itself is largely single-threaded per sample; increase mainly to scatter across more samples concurrently
- `memory_gb`: Both tasks load the full count matrix into memory right after `read10xCounts` (needed for `saveRDS` and for `emptyDrops` performance), so increase for matrices with very large numbers of raw barcodes (e.g. >1M droplets) or high gene counts

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for detailed guidelines

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
- **[OSCA Book](https://bioconductor.org/books/release/OSCA/)**: Orchestrating Single-Cell Analysis with Bioconductor, the reference workflow this module chain follows

