# ww-bioc-sc Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline for standard single-cell RNA-seq post-processing using the Bioconductor/OSCA (Orchestrating Single-Cell Analysis) ecosystem.

## Overview

This pipeline combines four Bioconductor-based WILDS modules to take 10x Genomics feature-barcode matrices from raw counts to annotated cell types: `ww-dropletutils`, `ww-scran`, `ww-scater`, and `ww-singler`. It follows the standard OSCA single-cell post-processing chain: load matrix -> QC filter -> normalize -> HVG -> PCA -> UMAP -> clustering -> marker genes -> cell-type annotation.

This pipeline serves as both a functional workflow and a demonstration of modular WDL design patterns within the WILDS ecosystem, including reuse of an upstream module's intermediate output (`ww-scran`'s highly variable gene list) by a downstream module (`ww-scater`) to avoid redundant computation.

**Complexity level:** Intermediate (4 modules)

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and demonstrates:

- **Module Integration**: Chaining four single-cell Bioconductor modules into one workflow
- **Data Flow**: Passing a `SingleCellExperiment` RDS object through QC, normalization, dimensionality reduction, and annotation steps
- **Cross-Module Reuse**: `ww-scater` reuses `ww-scran`'s highly variable gene selection via the `hvg_list` output/input, instead of recomputing it
- **Scatter-Gather**: Parallel per-sample processing of multiple single-cell samples

## Pipeline Steps

1. **Load Matrix** (using `ww-dropletutils` module):
   - For pre-filtered samples (`h5_matrix`): loads the 10x Genomics filtered feature-barcode matrix (H5 format) directly into a `SingleCellExperiment` object
   - For raw samples (`raw_h5_matrix`): calls real cells from the unfiltered matrix with `emptyDrops` first, then loads the result into a `SingleCellExperiment` object

2. **QC, Normalization, and HVG Selection** (using `ww-scran` module):
   - Computes per-cell QC metrics and flags low-quality cells
   - Normalizes with scran's pooling/deconvolution method
   - Models mean-variance trends and selects highly variable genes

3. **Dimensionality Reduction** (using `ww-scater` module):
   - Runs PCA on the highly variable genes selected by `ww-scran` (reused via `hvg_list`)
   - Runs UMAP on the PCA embedding
   - Generates QC and reduced-dimension diagnostic plots

4. **Clustering, Marker Genes, and Cell-Type Annotation** (using `ww-singler` module):
   - Clusters cells using the PCA embedding
   - Identifies per-cluster marker genes
   - Annotates each cluster's cell type against a celldex reference dataset (or a custom reference)

## Module Dependencies

This pipeline imports and uses:
- **ww-dropletutils module**: For loading 10x feature-barcode matrices (`read10x_counts` task) and calling real cells from raw matrices (`empty_drops_filter` task)
- **ww-scran module**: For QC filtering, normalization, and HVG selection (`run_scran` task)
- **ww-scater module**: For PCA/UMAP dimensionality reduction and QC visualization (`run_scater` task)
- **ww-singler module**: For clustering, marker gene identification, and cell-type annotation (`run_singler` task)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for module imports and celldex reference dataset downloads

### Input Configuration

Create an inputs JSON file with your samples:

```json
{
  "bioc_sc.samples": [
    {
      "name": "sample1",
      "h5_matrix": "/path/to/sample1_filtered_feature_bc_matrix.h5"
    },
    {
      "name": "sample2",
      "raw_h5_matrix": "/path/to/sample2_raw_feature_bc_matrix.h5"
    }
  ],
  "bioc_sc.mito_pattern": "^MT-",
  "bioc_sc.reference_dataset": "HumanPrimaryCellAtlasData"
}
```

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-bioc-sc.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-bioc-sc.wdl -i inputs.json

# Using Sprocket
sprocket run ww-bioc-sc.wdl @inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/datademos/proof-how-to/) to submit this pipeline directly to the on-premise HPC cluster:

1. Ensure your input files are accessible by the cluster
2. Update the inputs JSON with your sample information
3. Submit through the PROOF interface

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `samples` | List of SingleCellSample objects, each with a name and either a pre-filtered (`h5_matrix`) or raw (`raw_h5_matrix`) 10x feature-barcode matrix (H5) | Array[SingleCellSample] | Yes | - |
| `mito_pattern` | Regex pattern identifying mitochondrial gene symbols (passed to `ww-scran`) | String | No | `^MT-` |
| `nmads` | Number of median absolute deviations used to flag low-quality cells (passed to `ww-scran`) | Float | No | 3.0 |
| `n_hvgs` | Number of top highly variable genes to select (passed to `ww-scran` and used as a fallback in `ww-scater`) | Int | No | 2000 |
| `n_pcs` | Number of principal components to compute (passed to `ww-scater`) | Int | No | 50 |
| `reference_dataset` | Name of the celldex reference dataset function to fetch for cell-type annotation (passed to `ww-singler`) | String | No | `HumanPrimaryCellAtlasData` |
| `reference_ensembl` | Fetch the celldex reference with Ensembl gene IDs instead of gene symbols (passed to `ww-singler`) | Boolean | No | `true` |
| `label_column` | Column name in the reference object's colData containing cell-type labels (passed to `ww-singler`) | String | No | `label.main` |
| `cpu_cores` | Number of CPU cores allocated for `ww-dropletutils`, `ww-scran`, and `ww-scater` tasks | Int | No | 2 |
| `memory_gb` | Memory allocated in GB for `ww-dropletutils`, `ww-scran`, and `ww-scater` tasks | Int | No | 8 |
| `singler_memory_gb` | Memory allocated in GB for `ww-singler` tasks (higher default: celldex reference loading and clustering need more headroom) | Int | No | 16 |

### SingleCellSample Structure

Provide exactly one of `h5_matrix` or `raw_h5_matrix` per sample.

**Pre-filtered matrix:**
```json
{
  "name": "sample_name",
  "h5_matrix": "/path/to/filtered_feature_bc_matrix.h5"
}
```

**Raw (unfiltered) matrix** (real cells are called with `emptyDrops` first):
```json
{
  "name": "sample_name",
  "raw_h5_matrix": "/path/to/raw_feature_bc_matrix.h5"
}
```

## Output Files

The pipeline produces arrays of outputs (one per sample) from all four modules:

| Output | Description | Source Module |
|--------|-------------|---------------|
| `sce_rds` | Loaded `SingleCellExperiment` RDS objects, prior to QC/normalization | ww-dropletutils |
| `empty_drops_csv` | Per-barcode `emptyDrops` statistics, for samples provided as a raw matrix (null for pre-filtered samples) | ww-dropletutils |
| `barcode_rank_plot` | Barcode-rank QC plots, for samples provided as a raw matrix (null for pre-filtered samples) | ww-dropletutils |
| `scran_sce_object` | Normalized `SingleCellExperiment` RDS objects with size factors and PCA | ww-scran |
| `scran_qc_plot` | Per-cell QC scatter plots from scran | ww-scran |
| `scran_size_factor_plot` | Size factor histograms | ww-scran |
| `scran_mean_variance_plot` | Mean-variance trend plots | ww-scran |
| `scran_hvg_table` | CSVs of genes ranked by biological variance | ww-scran |
| `scater_sce_object` | `SingleCellExperiment` RDS objects with PCA and UMAP embeddings | ww-scater |
| `scater_pca_plot` | PCA scatter plots | ww-scater |
| `scater_umap_plot` | UMAP scatter plots | ww-scater |
| `scater_qc_plot` | Per-cell QC scatter plots from scater | ww-scater |
| `singler_sce_object` | Final `SingleCellExperiment` RDS objects with cluster assignments and predicted cell types | ww-singler |
| `singler_cluster_table` | CSVs mapping cell barcodes to clusters | ww-singler |
| `singler_marker_table` | CSVs of top marker genes per cluster | ww-singler |
| `singler_prediction_table` | CSVs of SingleR cell-type predictions per cluster | ww-singler |

## Resource Considerations

### Compute Requirements
- **Memory**: 8GB recommended for `ww-dropletutils`/`ww-scran`/`ww-scater` tasks; 16GB recommended for `ww-singler` tasks (celldex reference loading)
- **CPUs**: 2 cores per task is sufficient for most single-cell datasets; these tools are largely single-threaded per sample
- **Network**: Stable internet connection for module imports and celldex reference dataset downloads (cached after first fetch by the underlying ExperimentHub, but not across separate task invocations)

### Optimization Tips
- Scatter parallelism means resource requirements apply per sample, not per pipeline run; total resource usage scales with the number of samples processed concurrently
- Provide a custom reference `SingleCellExperiment` via `ww-singler`'s underlying `reference_rds` input (not exposed at the pipeline level in v1) if celldex's built-in references don't fit your organism or cell types

## Testing the Pipeline

The pipeline includes a test workflow that downloads a public 10x Genomics PBMC dataset:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint bioc_sc_example
```

The test workflow automatically:
1. Downloads a raw (unfiltered) 10x Genomics H5 matrix (human PBMC sample) via `ww-testdata`
2. Runs the full `ww-bioc-sc` pipeline on it via the `raw_h5_matrix` path (calling real cells with `emptyDrops` before QC)
3. Validates all output files

**Test Dataset Details:**
- Uses the human PBMC raw-matrix fixture from `ww-testdata`, paired with the default `HumanPrimaryCellAtlasData` celldex reference, so `ww-singler`'s cell-type predictions are biologically meaningful for this demo run

## Integration Patterns

This pipeline demonstrates several key WDL patterns:
- **Module Composition**: Chaining four analysis modules into a single coherent workflow
- **Cross-Module Data Reuse**: Passing an intermediate result (`ww-scran`'s `hvg_list`) between modules to avoid redundant computation, while keeping each module independently usable on its own
- **Struct Usage**: Grouping per-sample inputs (name, matrix file) into a single struct
- **Scatter-Gather**: Parallel processing of multiple samples with array-based outputs

## Extending the Pipeline

This pipeline can be extended by:
- Exposing `ww-singler`'s custom `reference_rds` input at the pipeline level for organism/cell-type combinations not covered by celldex
- Adding multi-sample integration (e.g. batch correction) before or after clustering
- Adding doublet detection upstream of QC filtering
- Adding ambient RNA correction (e.g. `ww-cellbender`) upstream of `ww-dropletutils`

## Related WILDS Components

- **ww-dropletutils module**: 10x matrix loading and empty droplet calling
- **ww-scran module**: QC filtering, normalization, and HVG selection
- **ww-scater module**: PCA/UMAP dimensionality reduction and visualization
- **ww-singler module**: Clustering, marker genes, and cell-type annotation
- **Other pipelines**: Additional integration examples

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
