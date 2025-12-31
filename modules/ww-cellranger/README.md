# ww-cellranger

[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for single-cell RNA-seq analysis using Cell Ranger.

## Overview

This module provides reusable WDL tasks for preparing and processing single-cell RNA-seq data using Cell Ranger from 10x Genomics.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `run_count`, `prepare_fastqs`
- **Test workflow**: `testrun.wdl` (demonstration workflow with automatic test data support)
- **Container**: `getwilds/cellranger:10.0.0` (WILDS Docker image with Cell Ranger installed)

## Important Requirements

### Platform Requirements

**Cell Ranger requires Linux x86_64 architecture with AVX instruction support.** This module runs Cell Ranger inside a Docker container (`getwilds/cellranger:10.0.0`), which works on:
- Linux x86_64 systems (native)
- Intel-based Macs (via Docker)
- Cloud platforms and HPC clusters (GitHub Actions, AWS, Google Cloud, etc.)

**Note for Apple Silicon (M1/M2/M3) Mac users:** Local testing is not supported. Docker's x86_64 emulation on Apple Silicon does not support the AVX instructions that Cell Ranger requires. The CI/CD pipeline will test this module on compatible infrastructure.

### FASTQ Naming Convention

**Your FASTQ files must follow Cell Ranger's naming convention:**
- Format: `SampleName_S1_L001_R1_001.fastq.gz` (for Read 1)
- Format: `SampleName_S1_L001_R2_001.fastq.gz` (for Read 2)

If your files don't follow this convention and you have just one pair of files, you use the `prepare_fastqs` task to rename them automatically.

### Current Limitations

**This module currently only accepts samples that are gene expression (GEX) only.** Feature barcoding data (e.g., antibody-derived tags, CRISPR guides) is not supported at this time.

If you need support for feature barcoding or other Cell Ranger features, please file a feature request on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Tasks

### `run_count`

Run `cellranger count` on gene expression reads from one GEM well.

**Inputs:**
- `r1_fastqs` (Array[File]): Array of R1 FASTQ files (contain cell barcodes and UMIs)
- `r2_fastqs` (Array[File]): Array of R2 FASTQ files (contain cDNA sequences)
- `ref_gex` (File): GEX reference transcriptome tarball (e.g., from 10x Genomics)
- `sample_id` (String): Sample ID for output naming
- `create_bam` (Boolean, default=true): Generate BAM file
- `cpu_cores` (Int, default=8): Number of CPU cores to use
- `memory_gb` (Int, default=64): Memory allocation in GB
- `expect_cells` (Int, optional): Expected number of recovered cells
- `chemistry` (String, optional): Assay configuration (e.g., SC3Pv3)

**Important:** All input FASTQs must be from one GEM well. If you have multiple GEM wells, use `run_count` separately for each well. The task validates that FASTQ filenames follow Cell Ranger's naming convention and will fail with a helpful error message if they don't.

**Outputs:**
- `results_tar` (File): Compressed tarball of Cell Ranger count output directory
- `web_summary` (File): Web summary HTML file with QC metrics
- `metrics_summary` (File): Metrics summary CSV file with key statistics

### `prepare_fastqs`

Rename FASTQs to Cell Ranger convention: `<sample_name>_S1_L001_R1_001.fastq.gz`. Useful when working with FASTQs from SRA or other sources that don't follow Cell Ranger's naming requirements.

**Inputs:**
- `r1_fastqs` (Array[File]): Array of R1 FASTQ files
- `r2_fastqs` (Array[File]): Array of R2 FASTQ files
- `sample_name` (String): Sample name for FASTQ naming

**Outputs:**
- `renamed_r1_fastqs` (Array[File]): R1 FASTQ files renamed to Cell Ranger convention
- `renamed_r2_fastqs` (Array[File]): R2 FASTQ files renamed to Cell Ranger convention

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellranger/ww-cellranger.wdl" as cellranger_tasks

workflow my_single_cell_pipeline {
  input {
    Array[File] r1_fastqs  # Must follow Cell Ranger naming convention
    Array[File] r2_fastqs  # Must follow Cell Ranger naming convention
    String sample_name
    File gex_reference
  }

  # Run Cell Ranger count
  call cellranger_tasks.run_count {
    input:
      r1_fastqs = r1_fastqs,
      r2_fastqs = r2_fastqs,
      ref_gex = gex_reference,
      sample_id = sample_name
  }

  output {
    File results = run_count.results_tar
    File web_summary = run_count.web_summary
    File metrics = run_count.metrics_summary
  }
}
```

### Advanced Usage Examples

**Custom resource allocation:**
```wdl
call cellranger_tasks.run_count {
  input:
    r1_fastqs = r1_fastqs,
    r2_fastqs = r2_fastqs,
    ref_gex = reference,
    sample_id = "my_sample",
    cpu_cores = 16,
    memory_gb = 128
}
```

**Specifying chemistry and expected cells:**
```wdl
call cellranger_tasks.run_count {
  input:
    r1_fastqs = r1_fastqs,
    r2_fastqs = r2_fastqs,
    ref_gex = reference,
    sample_id = "my_sample",
    chemistry = "SC3Pv3",
    expect_cells = 5000
}
```

**Without BAM generation (faster, less storage):**
```wdl
call cellranger_tasks.run_count {
  input:
    r1_fastqs = r1_fastqs,
    r2_fastqs = r2_fastqs,
    ref_gex = reference,
    sample_id = "my_sample",
    create_bam = false
}
```

**Using prepare_fastqs for SRA data:**
```wdl
# If your FASTQs don't follow Cell Ranger naming convention
call cellranger_tasks.prepare_fastqs {
  input:
    r1_fastqs = sra_r1_files,
    r2_fastqs = sra_r2_files,
    sample_name = "my_sample"
}

call cellranger_tasks.run_count {
  input:
    r1_fastqs = prepare_fastqs.renamed_r1_fastqs,
    r2_fastqs = prepare_fastqs.renamed_r2_fastqs,
    ref_gex = reference,
    sample_id = "my_sample"
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-sra**: Download public single-cell sequencing data before analysis

## Testing the Module

The module includes a demonstration workflow that can be tested independently. The workflow in `testrun.wdl` automatically downloads test data and runs without requiring input files:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

The test workflow (`cellranger_example`) automatically:
1. Downloads a small GEX reference using `ww-testdata`
2. Downloads demonstration FASTQ data from SRA using `ww-sra`
3. Renames FASTQs to Cell Ranger naming convention using `prepare_fastqs`
4. Runs Cell Ranger count analysis
5. Validates all outputs

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 16-128 GB recommended (scales with dataset size)
  - Small datasets (<5000 cells): 16-32 GB
  - Medium datasets (5000-10000 cells): 32-64 GB
  - Large datasets (>10000 cells): 64-128 GB
- **CPUs**: 4-16 cores recommended; Cell Ranger benefits from multi-threading
- **Disk space**: Ensure sufficient space for output files
  - With BAM: ~2-5x the size of input FASTQs
  - Without BAM: ~0.5-1x the size of input FASTQs

### Advanced Parameters

- `chemistry`: Specify the assay configuration if auto-detection fails (e.g., "SC3Pv2", "SC3Pv3")
- `expect_cells`: Provide expected cell count if desired
- `create_bam`: Set to `false` to save storage and processing time if BAM files are truly not needed.

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker or Apptainer support for containerized execution
- Sufficient computational resources (Cell Ranger can be memory-intensive)
- 10x Genomics reference transcriptome (available from [10x Genomics Download Center](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads))

## Performance Considerations

- **Memory usage**: Cell Ranger requires substantial RAM; allocate 64GB+ for typical datasets
- **CPU scaling**: Cell Ranger benefits significantly from multiple cores (8-16 recommended)
- **Storage requirements**:
  - Ensure sufficient space for outputs (2-5x input size with BAM)
  - Consider setting `create_bam = false` if BAM files are truly not needed
- **Runtime**: Typical runtime ranges greatly from 1-12 hours depending on dataset size and resources

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real single-cell sequencing data from SRA
- Comprehensive validation of all outputs
- Integration testing with ww-testdata and ww-sra modules

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to Cell Ranger usage or configuration, please refer to the [Cell Ranger documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
