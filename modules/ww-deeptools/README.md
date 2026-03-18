# ww-deeptools Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for [deepTools](https://deeptools.readthedocs.io/) — a suite of tools for exploring deep sequencing data. deepTools provides utilities for generating normalized coverage tracks, comparing signal between samples, computing matrices over genomic regions, and creating publication-quality heatmaps, profiles, and QC plots.

## Overview

This module wraps the most commonly used deepTools subcommands for analyzing high-throughput sequencing data, particularly ChIP-seq, ATAC-seq, and RNA-seq experiments. It supports the full workflow from BAM coverage generation through matrix computation to visualization.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-deeptools.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `bam_coverage`

Generates a normalized coverage track (bigWig or bedGraph) from a BAM file.

**Inputs:**
- `bam` (File): Input BAM file (sorted and indexed)
- `bai` (File): BAM index file
- `sample_name` (String): Name identifier for the sample
- `output_format` (String, default="bigwig"): Output format (bigwig or bedgraph)
- `bin_size` (Int, default=50): Output bin size in bases
- `normalize_using` (String, default="RPKM"): Normalization method (RPKM, CPM, BPM, RPGC, None)
- `effective_genome_size` (Int, optional): Effective genome size for RPGC normalization
- `extend_reads` (Int, optional): Fragment length to extend reads to
- `ignore_duplicates` (Boolean, default=false): Ignore duplicate reads
- `min_mapping_quality` (Int, default=0): Minimum mapping quality filter
- `extra_args` (String, default=""): Additional bamCoverage arguments
- `cpu_cores` (Int, default=4): CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `coverage_file` (File): Normalized coverage track file

### `bam_compare`

Compares two BAM files (e.g., treatment vs. control) and generates a comparison track.

**Inputs:**
- `treatment_bam` (File): Treatment/ChIP BAM file
- `treatment_bai` (File): Treatment BAM index
- `control_bam` (File): Control/input BAM file
- `control_bai` (File): Control BAM index
- `sample_name` (String): Output name identifier
- `operation` (String, default="log2"): Comparison operation (log2, ratio, subtract, add, mean, etc.)
- `normalize_using` (String, default="RPKM"): Normalization method
- `extra_args` (String, default=""): Additional bamCompare arguments
- `cpu_cores` (Int, default=4): CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `comparison_file` (File): Coverage comparison track file

### `compute_matrix`

Computes a matrix of signal values from bigWig files over genomic regions for visualization.

**Inputs:**
- `bigwig_files` (Array[File]): BigWig score files
- `regions_file` (File): BED or GTF file defining regions of interest
- `sample_name` (String): Output name identifier
- `mode` (String, default="reference-point"): Mode (reference-point or scale-regions)
- `reference_point` (String, default="TSS"): Reference point (TSS, TES, center)
- `before_region` (Int, default=3000): Distance upstream (bp)
- `after_region` (Int, default=3000): Distance downstream (bp)
- `bin_size` (Int, default=10): Signal averaging bin size
- `extra_args` (String, default=""): Additional computeMatrix arguments
- `cpu_cores` (Int, default=4): CPU cores
- `memory_gb` (Int, default=16): Memory in GB

**Outputs:**
- `matrix_gz` (File): Compressed matrix for plotting
- `matrix_tab` (File): Tab-separated matrix values

### `plot_heatmap`

Creates a heatmap visualization of signal enrichment over genomic regions.

**Inputs:**
- `matrix_gz` (File): Compressed matrix from computeMatrix
- `sample_name` (String): Output name identifier
- `plot_format` (String, default="png"): Image format (png, pdf, svg, eps)
- `color_map` (String, default="RdYlBu"): Matplotlib colormap name
- `extra_args` (String, default=""): Additional plotHeatmap arguments
- `cpu_cores` (Int, default=1): CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `heatmap` (File): Heatmap image file

### `plot_profile`

Creates a profile plot of average signal enrichment over genomic regions.

**Inputs:**
- `matrix_gz` (File): Compressed matrix from computeMatrix
- `sample_name` (String): Output name identifier
- `plot_format` (String, default="png"): Image format (png, pdf, svg, eps)
- `extra_args` (String, default=""): Additional plotProfile arguments
- `cpu_cores` (Int, default=1): CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `profile` (File): Profile plot image file

### `multi_bam_summary`

Computes read coverage summary across multiple BAM files for sample comparison.

**Inputs:**
- `bam_files` (Array[File]): Sorted and indexed BAM files
- `bai_files` (Array[File]): Corresponding BAM index files
- `sample_name` (String): Output name identifier
- `mode` (String, default="bins"): Mode (bins or BED-file)
- `bin_size` (Int, default=10000): Bin size for bins mode
- `regions_file` (File, optional): BED file for BED-file mode
- `extra_args` (String, default=""): Additional multiBamSummary arguments
- `cpu_cores` (Int, default=4): CPU cores
- `memory_gb` (Int, default=16): Memory in GB

**Outputs:**
- `summary_npz` (File): Compressed numpy array for analysis
- `raw_counts` (File): Tab-separated raw read counts

### `plot_correlation`

Generates a correlation heatmap or scatterplot from multiBamSummary output.

**Inputs:**
- `summary_npz` (File): Numpy array from multiBamSummary
- `sample_name` (String): Output name identifier
- `correlation_method` (String, default="spearman"): Method (spearman or pearson)
- `plot_type` (String, default="heatmap"): Plot type (heatmap or scatterplot)
- `plot_format` (String, default="png"): Image format
- `extra_args` (String, default=""): Additional plotCorrelation arguments
- `cpu_cores` (Int, default=1): CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `correlation_plot` (File): Correlation plot image
- `correlation_matrix` (File): Tab-separated correlation values

### `plot_pca`

Generates a PCA plot from multiBamSummary output.

**Inputs:**
- `summary_npz` (File): Numpy array from multiBamSummary
- `sample_name` (String): Output name identifier
- `plot_format` (String, default="png"): Image format
- `extra_args` (String, default=""): Additional plotPCA arguments
- `cpu_cores` (Int, default=1): CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `pca_plot` (File): PCA plot image
- `pca_data` (File): Tab-separated PCA eigenvalues and coordinates

### `plot_fingerprint`

Generates a fingerprint plot to assess ChIP enrichment quality.

**Inputs:**
- `bam_files` (Array[File]): Sorted and indexed BAM files
- `bai_files` (Array[File]): Corresponding BAM index files
- `sample_name` (String): Output name identifier
- `labels` (String, optional): Space-separated labels for each BAM file
- `plot_format` (String, default="png"): Image format
- `extra_args` (String, default=""): Additional plotFingerprint arguments
- `cpu_cores` (Int, default=4): CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `fingerprint_plot` (File): Fingerprint plot image
- `quality_metrics` (File): Tab-separated quality metrics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deeptools/ww-deeptools.wdl" as deeptools_tasks

workflow chip_seq_analysis {
  input {
    File treatment_bam
    File treatment_bai
    File control_bam
    File control_bai
    File regions_bed
  }

  # Generate coverage track
  call deeptools_tasks.bam_coverage { input:
      bam = treatment_bam,
      bai = treatment_bai,
      sample_name = "treatment",
      normalize_using = "CPM"
  }

  # Compare treatment vs control
  call deeptools_tasks.bam_compare { input:
      treatment_bam = treatment_bam,
      treatment_bai = treatment_bai,
      control_bam = control_bam,
      control_bai = control_bai,
      sample_name = "chip_vs_input"
  }

  # Compute matrix and generate heatmap
  call deeptools_tasks.compute_matrix { input:
      bigwig_files = [bam_coverage.coverage_file],
      regions_file = regions_bed,
      sample_name = "analysis"
  }

  call deeptools_tasks.plot_heatmap { input:
      matrix_gz = compute_matrix.matrix_gz,
      sample_name = "analysis"
  }
}
```

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint deeptools_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads reference data and test BAM files using `ww-testdata`
2. Generates normalized coverage tracks with `bamCoverage`
3. Computes a signal matrix over genomic regions with `computeMatrix`
4. Creates heatmap and profile visualizations
5. Runs multi-sample QC with `multiBamSummary`, `plotCorrelation`, `plotPCA`, and `plotFingerprint`

## Docker Container

This module uses the `getwilds/deeptools:3.5.6` container image, which includes:
- deepTools 3.5.6
- Python with numpy, scipy, matplotlib, and other dependencies
- All necessary system dependencies for coverage computation and plotting

## Citation

> Ramírez F, Ryan DP, Grüning B, Bhatt V, Luo Y, Lee BT, Silla-Martínez JM, Nöbel S, Muiño JM, Wurmus R, Dündar F, Manke T
> deepTools2: a next generation web server for deep-sequencing data analysis
> Nucleic Acids Research (2016) 44(W1):W160-W165
> DOI: 10.1093/nar/gkw257

## Parameters and Resource Requirements

### Default Resources
- **bamCoverage/bamCompare**: 4 CPU, 8 GB RAM
- **computeMatrix**: 4 CPU, 16 GB RAM
- **multiBamSummary**: 4 CPU, 16 GB RAM
- **Plotting tasks**: 1 CPU, 8 GB RAM
- **plotFingerprint**: 4 CPU, 8 GB RAM

### Resource Scaling
- `cpu_cores`: bamCoverage, bamCompare, computeMatrix, and multiBamSummary benefit from multiple cores
- `memory_gb`: Increase for large genomes or many samples; computeMatrix and multiBamSummary are memory-intensive

## Contributing

1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for detailed guidelines

## Related Resources

- **[deepTools Documentation](https://deeptools.readthedocs.io/)**: Official deepTools documentation
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
