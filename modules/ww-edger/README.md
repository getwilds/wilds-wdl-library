# ww-edger
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for differential expression analysis using edgeR.

## Overview

This module provides a reusable WDL task for differential gene expression analysis using edgeR's quasi-likelihood (QL) pipeline. It takes a combined count matrix and sample metadata as input, applies TMM normalization, filters low-count genes with `filterByExpr`, fits a QL negative binomial model, and outputs comprehensive results and visualizations.

The module is designed to accept the same count matrix format produced by `ww-deseq2`'s `combine_count_matrices` task, making it easy to run both tools side-by-side for comparison.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `run_edger`
- **Scripts**: `edger_analysis.R` (fetched via curl at runtime)
- **Test workflow**: `testrun.wdl` (demonstration workflow using Pasilla test data)
- **Containers**: `getwilds/edger:4.10.0` (edgeR analysis), `getwilds/python-utils:0.1.0` (count matrix prep and validation)
- **Dependencies**: Integrates with `ww-testdata` module for test data generation

## Scripts

All scripts are fetched from this repository at runtime via `curl`. This keeps the Docker images generic and makes it easy to iterate on analysis logic without rebuilding containers.

| Script | Used by | Language | Description |
|--------|---------|----------|-------------|
| [`edger_analysis.R`](edger_analysis.R) | `run_edger` | R | Runs the full edgeR QL pipeline: reads counts and metadata, applies `filterByExpr` filtering, calculates TMM normalization factors, estimates dispersion, fits a QL model, performs QL F-tests, and generates MD/BCV/volcano/heatmap plots with results CSVs |

## Tasks

### `run_edger`

Performs differential expression analysis using edgeR's quasi-likelihood pipeline with TMM normalization.

The analysis follows the recommended edgeR QL workflow:
1. `DGEList` construction with TMM normalization via `calcNormFactors`
2. Low-count gene filtering via `filterByExpr`
3. Dispersion estimation via `estimateDisp`
4. QL model fitting via `glmQLFit`
5. QL F-test via `glmQLFTest`
6. Results extraction via `topTags`

**Inputs:**
- `counts_matrix` (File): Combined matrix of gene-level counts (tab-delimited, genes as rows, samples as columns)
- `sample_metadata` (File): Sample metadata with experimental conditions (tab-delimited)
- `condition_column` (String): Column name containing experimental conditions (default: "condition")
- `reference_level` (String): Reference/control level for the contrast (default: first alphabetically)
- `contrast` (String): Contrast in format `treatment,control`; empty = use second vs. first factor level (default: "")
- `fdr_threshold` (Float): FDR threshold for calling significant genes (default: 0.05)
- `lfc_threshold` (Float): Minimum absolute log2 fold change for calling significant genes (default: 1.0)
- `min_count` (Int): Minimum count per sample passed to `filterByExpr` (default: 10)
- `min_total_count` (Int): Minimum total count across all samples passed to `filterByExpr` (default: 15)
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores (default: 2)
- `docker_image` (String): Docker image (default: `getwilds/edger:4.10.0`)

**Outputs:**
- `edger_results` (File): Complete edgeR results with statistics for all genes (logFC, logCPM, F, PValue, FDR)
- `edger_significant` (File): Filtered results for genes meeting FDR and LFC thresholds
- `edger_normalized_counts` (File): TMM-normalized CPM values for all samples
- `edger_md_plot` (File): MD plot showing log fold change vs. mean log CPM
- `edger_volcano_plot` (File): Volcano plot showing log fold change vs. statistical significance
- `edger_heatmap` (File): Heatmap of top differentially expressed genes (up to 50)
- `edger_bcv_plot` (File): Biological coefficient of variation plot showing dispersion estimates

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-edger/ww-edger.wdl" as edger_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deseq2/ww-deseq2.wdl" as deseq2_tasks

workflow my_rna_seq_analysis {
  input {
    Array[File] star_count_files
    Array[String] sample_names
    Array[String] conditions
    String treatment_condition
    String control_condition
  }

  call deseq2_tasks.combine_count_matrices {
    input:
      gene_count_files = star_count_files,
      sample_names = sample_names,
      sample_conditions = conditions
  }

  call edger_tasks.run_edger {
    input:
      counts_matrix = combine_count_matrices.counts_matrix,
      sample_metadata = combine_count_matrices.sample_metadata,
      condition_column = "condition",
      reference_level = control_condition,
      contrast = "${treatment_condition},${control_condition}"
  }

  output {
    File de_results = run_edger.edger_results
    File significant_genes = run_edger.edger_significant
    File md_plot = run_edger.edger_md_plot
    File volcano_plot = run_edger.edger_volcano_plot
    File heatmap = run_edger.edger_heatmap
  }
}
```

### Advanced Usage Examples

**Explicit contrast specification:**
```wdl
call edger_tasks.run_edger {
  input:
    counts_matrix = my_counts,
    sample_metadata = my_metadata,
    condition_column = "treatment_group",
    reference_level = "vehicle_control",
    contrast = "drug_treated,vehicle_control"
}
```

**Relaxed significance thresholds:**
```wdl
call edger_tasks.run_edger {
  input:
    counts_matrix = my_counts,
    sample_metadata = my_metadata,
    fdr_threshold = 0.10,
    lfc_threshold = 0.5
}
```

**Resource optimization for large datasets:**
```wdl
call edger_tasks.run_edger {
  input:
    counts_matrix = my_counts,
    sample_metadata = my_metadata,
    memory_gb = 16,
    cpu_cores = 4
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-deseq2**: Use `combine_count_matrices` to prepare input; run edgeR and DESeq2 side-by-side for robust DE calling
- **ww-star**: Use gene count outputs from STAR alignment as input after combining with `combine_count_matrices`
- **ww-testdata**: Automatic provisioning of Pasilla test count data and metadata
- **Downstream analysis**: Results CSVs can feed into pathway enrichment or other analysis modules

## Testing the Module

The module includes a test workflow that runs with test data and requires no inputs:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

### Test Data Workflow

The test workflow automatically:
1. Generates the Pasilla test dataset (7 samples, 10,000 genes) using `ww-testdata`
2. Combines individual STAR count files into a count matrix
3. Performs differential expression analysis with edgeR
4. Generates all visualizations and results files
5. Validates all outputs

No input files or parameters are required.

## Configuration Guidelines

### Resource Allocation

- **Memory**: 4-16 GB recommended (scales with dataset size)
- **CPUs**: 1-4 cores recommended; edgeR is primarily single-threaded

### Statistical Parameters

- **contrast**: Specify as `treatment,control`; leave empty to compare second vs. first factor level
- **reference_level**: Set control/baseline condition for proper fold change direction
- **fdr_threshold**: 0.05 is standard; relax to 0.10 for exploratory analyses
- **lfc_threshold**: 1.0 (2-fold change) is a common practical cutoff; set to 0 for FDR-only filtering

### filterByExpr Parameters

- **min_count**: Minimum count in the smallest group (default 10 matches edgeR's recommendation)
- **min_total_count**: Minimum summed count across all samples (default 15); raise for larger datasets

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Internet access for fetching scripts from GitHub at runtime
- R environment with edgeR, ggplot2, pheatmap, RColorBrewer (provided by `getwilds/edger:4.10.0`)
- Python 3 with pandas (provided by `getwilds/python-utils:0.1.0` for count matrix prep and validation)

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to edgeR usage or interpretation, refer to the [edgeR user guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRusersguide.pdf). Please cite their work if you use edgeR in your analyses:

Chen Y, Lun ATL, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. F1000Research. 2016;5:1438.

Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26(1):139-140.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
