# ww-cellbender Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for removing ambient RNA contamination and barcode swapping artifacts from single-cell RNA sequencing data using [CellBender](https://cellbender.readthedocs.io/).

## Overview

CellBender uses a deep generative model to distinguish real cell signal from technical noise sources in UMI-based scRNA-seq, snRNA-seq, and CITE-seq experiments. Its `remove-background` command takes the raw (unfiltered) feature-barcode matrix produced by Cell Ranger and outputs a cleaned count matrix suitable for downstream analysis.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-cellbender.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `remove_background`

Runs CellBender `remove-background` on a raw 10x Genomics feature-barcode matrix to remove ambient RNA and barcode swapping noise.

**Inputs:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_h5` | File | required | Raw feature-barcode matrix HDF5 file from Cell Ranger (`raw_feature_bc_matrix.h5`) |
| `sample_name` | String | required | Sample name used as output file prefix |
| `expected_cells` | Int? | auto | Estimated number of real cells; helps CellBender set priors |
| `total_droplets_included` | Int? | auto | Total droplets to analyze (ranked by UMI count) |
| `fpr` | String | `"0.01"` | False positive rate target(s), space-separated for multiple (e.g. `"0.01 0.05"`) |
| `epochs` | Int | `150` | Number of training epochs |
| `learning_rate` | Float | `0.0001` | Base learning rate for the variational inference optimizer |
| `model` | String | `"full"` | CellBender model variant: `naive`, `simple`, `ambient`, `swapping`, or `full` |
| `low_count_threshold` | Int | `5` | Droplets below this UMI count are excluded from analysis |
| `exclude_feature_types` | String? | none | Space-separated feature types to exclude (e.g. `"Antibody Capture"`) |
| `checkpoint_mins` | Float | `7.0` | How frequently (in minutes) to save a training checkpoint |
| `gpu_enabled` | Boolean | `true` | Enable GPU acceleration (`--cuda`); set to `false` for CPU-only execution |
| `cpu_cores` | Int | `4` | Number of CPU cores allocated for the task |
| `memory_gb` | Int | `32` | Memory allocated for the task in GB |
| `docker_image` | String | `getwilds/cellbender:0.3.2` | Docker image to use for this task |

**Outputs:**

| Output | Type | Description |
|--------|------|-------------|
| `output_h5` | File | Full cleaned count matrix in H5 format (all input barcodes retained) |
| `filtered_h5` | File | Filtered cleaned count matrix (barcodes with >50% cell probability only) |
| `report_html` | File | Interactive HTML report with training diagnostics and quality plots |
| `summary_pdf` | File | PDF summary of the inference procedure |
| `log_file` | File | Run log with diagnostics and parameter summary |
| `cell_barcodes_csv` | File | Plain-text list of cell barcodes exceeding 50% cell probability threshold |
| `metrics_csv` | File | Run metrics CSV (one row per FPR value) |
| `checkpoint_tar` | File | Trained model checkpoint tarball (can be used to resume or inspect training) |

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellbender/ww-cellbender.wdl" as cellbender_tasks

workflow my_scrna_pipeline {
  input {
    File raw_h5
    String sample_name
    Int expected_cells
  }

  call cellbender_tasks.remove_background {
    input:
      input_h5 = raw_h5,
      sample_name = sample_name,
      expected_cells = expected_cells
  }

  output {
    File cleaned_h5      = remove_background.filtered_h5
    File qc_report       = remove_background.report_html
    File cell_barcodes   = remove_background.cell_barcodes_csv
  }
}
```

### Advanced Usage: Multiple FPR Values

```wdl
call cellbender_tasks.remove_background {
  input:
    input_h5 = raw_h5,
    sample_name = "my_sample",
    expected_cells = 5000,
    total_droplets_included = 15000,
    fpr = "0.01 0.05 0.1",
    epochs = 200,
    cpu_cores = 8,
    memory_gb = 64
}
```

### Integration with ww-cellranger

CellBender is typically run after Cell Ranger. The `raw_h5` output from `ww-cellranger` can be passed directly as input:

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellranger/ww-cellranger.wdl" as cellranger_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellbender/ww-cellbender.wdl" as cellbender_tasks

workflow scrna_pipeline {
  # ... Cell Ranger step ...
  call cellranger_tasks.cellranger_count { ... }

  call cellbender_tasks.remove_background {
    input:
      input_h5 = cellranger_count.raw_h5,
      sample_name = sample_name
  }
}
```

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that downloads a public 10x Genomics raw feature-barcode matrix and runs CellBender on it.

```bash
# Using Sprocket
sprocket run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

## Docker Container

This module uses the `getwilds/cellbender:0.3.2` container image from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library), which includes:

- CellBender 0.3.2
- Python 3 with PyTorch and all CellBender dependencies
- CPU-only runtime (GPU is not required but will accelerate training if available)

## Citation

> CellBender remove-background: a deep generative model for unsupervised removal of background noise from scRNA-seq datasets
> Stephen J Fleming, Mark D Chaffin, Alessandro Arduini, Amer-Denis Akkad, Eric Banks, John C Marioni, Anthony A Philippakis, Patrick T Ellinor, Mehrtash Babadi
> bioRxiv 2019.12.9.869750; doi: https://doi.org/10.1101/791699

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 4 cores
- **Memory**: 32 GB
- **Runtime**: Variable depending on dataset size and epoch count; typically 30-120 minutes per sample on CPU

### Resource Scaling
- For large datasets (>20,000 cells or >100,000 droplets), increase `memory_gb` to 64 GB or more
- Increasing `cpu_cores` up to 8 can speed up posterior computation
- The `epochs` parameter controls training time: fewer epochs are faster but may reduce accuracy; the default of 150 is suitable for most datasets

## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org

## Related Resources

- **[CellBender Documentation](https://cellbender.readthedocs.io/)**: Official CellBender documentation
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
