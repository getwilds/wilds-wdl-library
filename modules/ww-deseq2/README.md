# ww-deseq2
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for differential expression analysis using DESeq2.

## Overview

This module provides reusable WDL tasks for comprehensive differential gene expression analysis using DESeq2. It handles the complete workflow from count matrix preparation through statistical analysis and visualization generation. The module supports both single-factor and complex experimental designs with comprehensive output validation.

The module can run completely standalone with automatic test data generation. Individual tasks can be imported and used with existing count matrices and metadata files for production analyses.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `combine_count_matrices`, `run_deseq2`, `validate_outputs`
- **Workflow**: `deseq2_example` (demonstration workflow using test data)
- **Containers**: `getwilds/combine-counts:0.1.0`, `getwilds/deseq2:1.40.2`
- **Dependencies**: Integrates with `ww-testdata` module for test data generation
- **Test Data**: Uses Pasilla test dataset with 7 samples and 10,000 genes

## Tasks

### `combine_count_matrices`
Combines STAR gene count files from multiple samples into a single count matrix with accompanying metadata.

**Inputs:**
- `gene_count_files` (Array[File]): Array of STAR gene count files (ReadsPerGene.out.tab)
- `sample_names` (Array[String]): Sample names corresponding to the gene_count_files
- `sample_conditions` (Array[String]): Experimental conditions for each sample
- `memory_gb` (Int): Memory allocation in GB (default: 4)
- `cpu_cores` (Int): Number of CPU cores (default: 1)
- `count_column` (Int): Column number to extract from STAR files (default: 2)
  - 2 = unstranded counts
  - 3 = stranded counts, first read forward
  - 4 = stranded counts, first read reverse

**Outputs:**
- `counts_matrix` (File): Combined matrix of gene-level counts from all samples
- `sample_metadata` (File): Metadata file containing sample names and conditions

### `run_deseq2`
Performs differential expression analysis using DESeq2 with comprehensive statistical analysis and visualization via a prewritten R script: [deseq2_analysis.R](https://github.com/getwilds/wilds-docker-library/blob/main/deseq2/deseq2_analysis.R)

**Inputs:**
- `counts_matrix` (File): Combined matrix of gene-level counts
- `sample_metadata` (File): Sample metadata with experimental conditions
- `condition_column` (String): Column name containing experimental conditions (default: "condition")
- `reference_level` (String): Reference level for DESeq2 contrast (typically control condition)
- `contrast` (String): DESeq2 contrast string in format 'condition,treatment,control'
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores (default: 2)

**Outputs:**
- `deseq2_results` (File): Complete DESeq2 differential expression results with statistics (used to make volcano plot)
- `deseq2_significant` (File): Filtered results containing only statistically significant genes (used to make heatmap)
- `deseq2_normalized_counts` (File): DESeq2 normalized count values for all samples
- `deseq2_pca_plot` (File): Principal Component Analysis plot showing sample clustering
- `deseq2_volcano_plot` (File): Volcano plot showing log fold change vs. statistical significance
- `deseq2_heatmap` (File): Heatmap visualization of differentially expressed genes

### `validate_outputs`
Validates DESeq2 analysis outputs for correctness and completeness.

**Inputs:**
- `deseq2_results` (File): DESeq2 results file to validate
- `deseq2_significant` (File): Significant results file to validate
- `normalized_counts` (File): Normalized counts file to validate
- `expected_samples` (Int): Expected number of samples in the analysis
- `expected_genes_min` (Int): Minimum expected number of genes in results (default: 1000)

**Outputs:**
- `validation_report` (File): Report summarizing validation results and any issues found

## Usage as a Module

### Importing into Your Workflow

```wdl
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
  
  call deseq2_tasks.run_deseq2 {
    input:
      counts_matrix = combine_count_matrices.counts_matrix,
      sample_metadata = combine_count_matrices.sample_metadata,
      condition_column = "condition",
      reference_level = control_condition,
      contrast = "condition,${treatment_condition},${control_condition}"
  }
  
  call deseq2_tasks.validate_outputs {
    input:
      deseq2_results = run_deseq2.deseq2_results,
      deseq2_significant = run_deseq2.deseq2_significant,
      normalized_counts = run_deseq2.deseq2_normalized_counts,
      expected_samples = length(sample_names)
  }
  
  output {
    File differential_expression_results = run_deseq2.deseq2_results
    File significant_genes = run_deseq2.deseq2_significant
    File pca_plot = run_deseq2.deseq2_pca_plot
    File volcano_plot = run_deseq2.deseq2_volcano_plot
    File heatmap = run_deseq2.deseq2_heatmap
  }
}
```

### Advanced Usage Examples

**Using pre-combined count matrix:**
```wdl
call deseq2_tasks.run_deseq2 {
  input:
    counts_matrix = my_existing_matrix,
    sample_metadata = my_existing_metadata,
    condition_column = "treatment_group",
    reference_level = "control"
}
```

**Complex contrast specification:**
```wdl
call deseq2_tasks.run_deseq2 {
  input:
    contrast = "treatment,drug_treated,vehicle_control",
    reference_level = "vehicle_control"
}
```

**Resource optimization for large datasets:**
```wdl
call deseq2_tasks.run_deseq2 {
  input:
    memory_gb = 16,
    cpu_cores = 4
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-star**: Use gene count outputs from STAR alignment for differential expression
- **ww-testdata**: Automatic provisioning of test count data and metadata
- **Custom workflows**: Foundation for complete RNA-seq analysis pipelines
- **Downstream analysis**: Outputs can feed into pathway analysis or other modules

## Testing the Module

The module includes a demonstration workflow that runs with test data and requires no inputs:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-deseq2.wdl

# Using miniWDL
miniwdl run ww-deseq2.wdl

# Using Sprocket
sprocket run ww-deseq2.wdl
```

### Test Data Workflow

The `deseq2_example` workflow automatically:
1. Generates Pasilla test dataset using `ww-testdata`
2. Combines individual count files into a count matrix
3. Performs differential expression analysis with DESeq2
4. Generates comprehensive visualizations
5. Validates all outputs

No input files or parameters are required - the workflow uses pre-configured test data and parameters.

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 4-16 GB recommended (scales with dataset size and number of samples)
- **CPUs**: 1-4 cores recommended; DESeq2 benefits from modest parallelization
- **Count Column**: Choose appropriate column from STAR output:
  - Column 2 (unstranded): Most common for standard protocols
  - Column 3/4 (stranded): For stranded RNA-seq protocols

### Statistical Parameters

- **Contrast**: Specify comparisons in DESeq2 format (condition,numerator,denominator)
- **Reference Level**: Set control or baseline condition for proper fold change calculation
- **Condition Column**: Ensure metadata column name matches your experimental design

### Test Workflow Configuration

The `deseq2_example` workflow:
- Uses Pasilla dataset with 7 samples and 10,000 genes
- Demonstrates basic treatment vs. control comparison
- Runs with fixed parameters optimized for test data
- Requires no user inputs or configuration

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- R environment with DESeq2 package (provided in container)
- Python environment for count matrix combination (provided in container)

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to DESeq2 usage or interpretation, please refer to the [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). Please make sure to cite their work if you use DESeq2 in your analyses:

Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15(12):550.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
