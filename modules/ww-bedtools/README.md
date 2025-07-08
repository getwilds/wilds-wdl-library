# ww-bedtools
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for using BEDTools utilities.

## Overview

This module provides reusable WDL tasks for using **BEDTools** to work with genomic intervals. It includes comprehensive tools for coverage analysis, intersection operations, and genomic window creation with read counting.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger sequence analysis pipelines. It can run completely standalone with automatic test data download and alignment, or integrate with existing BAM files.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `coverage`, `intersect`, `makewindows`, `validate_outputs`
- **Workflow**: `bedtools_example` (demonstration workflow with automatic test data support)
- **Container**: `getwilds/bedtools:2.31.1`
- **Dependencies**: Integrates with `ww-sra`, `ww-bwa`, and `ww-testdata` modules for complete workflows
- **Test Data**: Automatically downloads reference genome, BED file, and SRA data when not provided

## Tasks

### `coverage`

Calculates average read coverage over given BED intervals for a sample.

**Inputs:**

- `bed_file` (File): BED file containing genomic intervals
- `aligned_bam` (File): Input aligned and indexed BAM file
- `sample_name` (String): Name of the sample provided for output files
- `cpu_cores` (Int): Number of CPU cores (default: 2)
- `memory_gb` (Int): Memory allocation in GB (default: 16)

**Outputs:**

- `name` (String): Sample name
- `mean_coverage` (File): Output file with mean coverage values

### `intersect`

Uses BEDTools to find overlaps between the sample's BAM and a BED file.

**Inputs:**

- `bed_file` (File): BED file to intersect with
- `aligned_bam` (File): Input aligned and indexed BAM file
- `sample_name` (String): Name of the sample provided for output files
- `flags` (String): BEDTools intersect flags (default: `-header -wo`)
- `cpu_cores` (Int): Number of CPU cores (default: 2)
- `memory_gb` (Int): Memory allocation in GB (default: 16)

**Outputs:**

- `name` (String): Sample name
- `intersect_output` (File): `bedtools intersect` result file

### `makewindows`

Creates genomic windows by chromosome and counts reads overlapping each window.

**Inputs:**

- `bed_file` (File): BED file of input intervals
- `aligned_bam` (File): Input aligned BAM file
- `bam_index` (File): Index of aligned BAM file
- `reference_fasta` (File): Reference genome FASTA
- `reference_index` (File): Reference genome index
- `list_chr` (Array[String]): List of chromosomes to analyze
- `sample_name` (String): Name of the sample provided for output files
- `tmp_dir` (String): Temporary directory path
- `cpu_cores` (Int): Number of CPU cores (default: 10)
- `memory_gb` (Int): Memory allocation in GB (default: 24)

**Outputs:**

- `name` (String): Sample name
- `counts_bed` (File): `.tar.gz` archive of per-chromosome BED files with read counts

### `validate_outputs`

Checks existence and content of all output files by sample.

**Inputs:**

- `intersect_files` (Array[File]): Intersect output files
- `coverage_files` (Array[File]): Mean coverage files
- `window_count_files` (Array[File]): BED file tarballs with read counts
- `sample_names` (Array[String]): Sample names corresponding to each file

**Outputs:**

- `report` (File): Human-readable validation report

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/main/modules/ww-bedtools/ww-bedtools.wdl" as bedtools

struct BedtoolsSample {
    String name
    File bam
    File bam_index
}

workflow my_pipeline {
  input {
    Array[BedtoolsSample] samples
    File bed_file
    File reference_fasta
    File reference_index
  }

  scatter (sample in samples) {
    call bedtools.coverage {
      input:
        bed_file = bed_file,
        aligned_bam = sample.bam,
        sample_name = sample.name
    }

    call bedtools.intersect {
      input:
        bed_file = bed_file,
        aligned_bam = sample.bam,
        sample_name = sample.name
    }
  }

  output {
    Array[File] coverage_results = coverage.mean_coverage
    Array[File] intersect_results = intersect.intersect_output
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-sra**: Download sequencing data prior to alignment (built into demo workflow)
- **ww-bwa**: Sequence alignment (built into demo workflow)
- **ww-star**: Use aligned BAM files from STAR for BEDTools analysis
- **Custom workflows**: Combine with variant callers or other genomic analysis tools

## Testing the Module

The module includes a demonstration workflow with comprehensive testing and **automatic test data download**:

```bash
# Using Cromwell (no input files required - test data downloaded automatically)
java -jar cromwell.jar run ww-bedtools.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-bedtools.wdl -i inputs.json

# Using Sprocket
sprocket run ww-bedtools.wdl inputs.json
```

### Automatic Demo Mode

When no samples, reference files, or BED file are provided, the workflow automatically:
1. Downloads reference genome data and BED file using `ww-testdata`
2. Downloads SRA data (default: ERR1258306) using `ww-sra`
3. Builds BWA index using `ww-bwa`
4. Aligns reads using `ww-bwa`
5. Runs all BEDTools analyses (coverage, intersect, makewindows)

### Test Input Format

**Minimal input (uses automatic demo data):**
```json
{
  "bedtools_example.demo_sra_id": "ERR1258306",
  "bedtools_example.intersect_flags": "-header -wo",
  "bedtools_example.chromosomes": ["chr22"],
  "bedtools_example.tmp_dir": "/tmp",
  "bedtools_example.cpus": 2,
  "bedtools_example.memory_gb": 8
}
```

**Full input (provide your own data):**
```json
{
  "bedtools_example.samples": [
    {
      "name": "sample1",
      "bam": "/path/to/sample1.bam",
      "bam_index": "/path/to/sample1.bam.bai"
    }
  ],
  "bedtools_example.bed_file": "/path/to/regions.bed",
  "bedtools_example.reference_fasta": "/path/to/genome.fasta",
  "bedtools_example.reference_index": "/path/to/genome.fasta.fai",
  "bedtools_example.intersect_flags": "-header -wo",
  "bedtools_example.chromosomes": [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
  ],
  "bedtools_example.tmp_dir": "/tmp",
  "bedtools_example.cpus": 2,
  "bedtools_example.memory_gb": 8
}
```

**Note**: You can mix and match - provide some inputs and let others use test data.

## Workflow: `bedtools_example`

The included example workflow processes BAM files using all BEDTools tasks and validates the results.

**Inputs:**

- `bed_file` (File?): Optional BED file of intervals (test data used if not provided)
- `samples` (Array[BedtoolsSample]?): Optional structs containing BAM, index, and sample name (demo data generated if not provided)
- `reference_fasta` (File?): Optional genome reference FASTA (test data used if not provided)
- `reference_index` (File?): Optional FAI index for reference (test data used if not provided)
- `demo_sra_id` (String): SRA ID for demo data (default: "ERR1258306")
- `intersect_flags` (String): Flags for `bedtools intersect` (default: `-header -wo`)
- `chromosomes` (Array[String]): Chromosomes for `makewindows` (default: all human chromosomes)
- `tmp_dir` (String): Temporary path for intermediate files (default: "/tmp")
- `cpus` (Int): Threads for each task (default: 2)
- `memory_gb` (Int): Memory for each task (default: 8)

**Outputs:**

- `intersect_results` (Array[File]): Output from `intersect`
- `coverage_results` (Array[File]): Output from `coverage`
- `window_counts` (Array[File]): Output tarballs from `makewindows`
- `validation_report` (File): Output report from `validate_outputs`

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:
- **Memory**: 8-24 GB recommended depending on genome size and data volume
- **CPUs**: 2 cores for coverage/intersect tasks; 10 cores recommended for `makewindows` due to parallelization
- **Storage**: Sufficient space for input BAMs, reference files, and output files

### Advanced Parameters

- `intersect_flags`: Customize BEDTools intersect behavior (e.g., "-wo" for overlap amounts)
- `chromosomes`: Limit analysis to specific chromosomes for faster processing  
- `tmp_dir`: Set appropriate temporary directory with sufficient space
- Ensure the reference genome is pre-indexed and has an associated `.fai`

### Demo Configuration

- `demo_sra_id`: Change to use different SRA sample for testing
- Consider using single chromosome (e.g., ["chr22"]) for faster demo runs
- Resource parameters apply to both demo and user-provided data modes

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker or Apptainer support
- Input BAM files must be sorted and indexed (when providing your own data)
- Reference genome with FASTA index (when providing your own data)
- Sufficient memory for genome processing (varies by genome size)

## Features

- **Standalone execution**: Complete workflow with automatic test data download
- **Flexible input**: Use your own data or automatic demo data
- **BEDTools intersect analysis**: Screen for overlaps between two sets of genomic features
- **BEDTools coverage analysis**: Read coverage results across BED intervals
- **BEDTools makewindows output**: Per-chromosome BED files of read counts with parallelization
- **Validation**: Built-in output validation and reporting
- **Module integration**: Seamlessly combines with ww-sra, ww-bwa, and ww-testdata
- **Modular design**: Integrates with other WILDS workflows and tools
- **Scalable**: Supports batch analysis across many samples
- **Flexible**: Customizable resource settings per task
- **Robust**: Includes error handling and reproducible output filenames

## Advanced Usage

### Chromosome-specific Analysis

For faster processing or targeted analysis, limit to specific chromosomes:

```json
{
  "bedtools_example.chromosomes": ["chr1", "chr2", "chr3"]
}
```

### Custom Intersect Flags

Modify BEDTools intersect behavior:

```json
{
  "bedtools_example.intersect_flags": "-u -f 0.50"
}
```

### Resource Optimization

Adjust resources based on data size:

```json
{
  "bedtools_example.cpus": 4,
  "bedtools_example.memory_gb": 16
}
```

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (SRA sample ERR1258306 for integration testing)
- Comprehensive validation of all outputs
- Integration testing with ww-sra, ww-bwa, and ww-testdata modules
- Chromosome 22 subset for efficiency during CI/CD

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, contact the Fred Hutch Data Science Lab (DaSL) at [wilds@fredhutch.org](mailto:wilds@fredhutch.org) or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

To contribute to this module, please review the [WILDS Contributor Guide](https://getwilds.org/guide/) and the [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md).

## License

Distributed under the MIT License. See `LICENSE` for full details.
