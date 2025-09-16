# ww-bedtools
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for working with genomic intervals using BEDTools.

## Overview

This module provides reusable WDL tasks for genomic interval analysis using BEDTools, a powerful suite of utilities for comparing, manipulating, and analyzing genomic features. The module includes tasks for coverage analysis, interval intersection, and window-based read counting across chromosomes.

The module is designed to be a foundational component within the WILDS ecosystem and includes automatic test data downloading when no inputs are provided.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `coverage`, `intersect`, `makewindows`, `validate_outputs`
- **Workflow**: `bedtools_example` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/bedtools:2.31.1`
- **Test Data**: Automatically downloads test data when no inputs provided using `ww-testdata` module

## Tasks

### `coverage`
Calculates mean read coverage across BED intervals using BEDTools coverage.

**Inputs:**
- `bed_file` (File): BED file containing genomic intervals
- `aligned_bam` (File): Input aligned and indexed BAM file
- `sample_name` (String): Sample name for output files
- `cpu_cores` (Int): CPU cores (default: 2)
- `memory_gb` (Int): Memory allocation (default: 16)

**Outputs:**
- `name` (String): Sample name that was processed
- `mean_coverage` (File): File containing mean read coverage across BED intervals

### `intersect`
Performs BEDTools intersect between BED intervals and BAM alignments.

**Inputs:**
- `bed_file` (File): BED file containing genomic intervals
- `aligned_bam` (File): Input aligned and indexed BAM file
- `sample_name` (String): Sample name for output files
- `flags` (String): BEDTools intersect command flags (default: "-header -wo")
- `cpu_cores` (Int): CPU cores (default: 2)
- `memory_gb` (Int): Memory allocation (default: 16)

**Outputs:**
- `name` (String): Sample name that was processed
- `intersect_output` (File): BEDTools intersect results file

### `makewindows`
Creates genomic windows and counts reads per window across specified chromosomes.

**Inputs:**
- `bed_file` (File): BED file containing genomic intervals
- `reference_fasta` (File): Reference genome FASTA file
- `reference_index` (File): Reference genome index file
- `aligned_bam` (File): Input aligned BAM file
- `bam_index` (File): Index of aligned BAM file
- `list_chr` (Array[String]): Array of chromosome names to process
- `sample_name` (String): Sample name for output files
- `tmp_dir` (String): Path to temporary directory
- `cpu_cores` (Int): CPU cores (default: 10)
- `memory_gb` (Int): Memory allocation (default: 24)

**Outputs:**
- `name` (String): Sample name that was processed
- `counts_bed` (File): Tarball of per-chromosome BED files with read counts

### `validate_outputs`
Validates all BEDTools outputs and generates comprehensive statistics.

**Inputs:**
- `intersect_files` (Array[File]): BEDTools intersect output files
- `coverage_files` (Array[File]): Coverage analysis results
- `window_count_files` (Array[File]): Window count tarballs
- `sample_names` (Array[String]): Sample names that were processed

**Outputs:**
- `report` (File): Validation report summarizing file check results

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bedtools/ww-bedtools.wdl" as bedtools_tasks

struct BedtoolsSample {
    String name
    File bam
    File bam_index
}

workflow my_interval_analysis_pipeline {
  input {
    Array[BedtoolsSample] samples
    File bed_file
    File reference_fasta
    File reference_index
    Array[String] chromosomes
  }
  
  scatter (sample in samples) {
    call bedtools_tasks.coverage {
      input:
        bed_file = bed_file,
        aligned_bam = sample.bam,
        sample_name = sample.name
    }
    
    call bedtools_tasks.intersect {
      input:
        bed_file = bed_file,
        aligned_bam = sample.bam,
        sample_name = sample.name,
        flags = "-header -wo"
    }
    
    call bedtools_tasks.makewindows {
      input:
        bed_file = bed_file,
        aligned_bam = sample.bam,
        bam_index = sample.bam_index,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        sample_name = sample.name,
        list_chr = chromosomes
    }
  }
  
  call bedtools_tasks.validate_outputs {
    input:
      intersect_files = intersect.intersect_output,
      coverage_files = coverage.mean_coverage,
      window_count_files = makewindows.counts_bed,
      sample_names = coverage.name
  }
  
  output {
    Array[File] intersect_results = intersect.intersect_output
    Array[File] coverage_results = coverage.mean_coverage
    Array[File] window_counts = makewindows.counts_bed
    File validation_report = validate_outputs.report
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Provides reference genomes, BED files, and test BAM files
- **ww-star**: Use aligned BAM files from STAR for interval analysis
- **ww-bwa**: Use aligned BAM files from BWA for interval analysis
- **Custom workflows**: Any workflow producing BAM files can use these interval analysis tasks

## Testing the Module

The module includes a demonstration workflow with comprehensive testing:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-bedtools.wdl

# Using miniWDL
miniwdl run ww-bedtools.wdl

# Using Sprocket
sprocket run ww-bedtools.wdl
```

### Automatic Test Data

The demonstration workflow automatically downloads all required test data using the `ww-testdata` module, including BED files, BAM files, and reference genome data.

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration per task:
- **Coverage task**: 2 cores, 16GB RAM (typical)
- **Intersect task**: 2 cores, 16GB RAM (typical)  
- **Makewindows task**: 10 cores, 24GB RAM (parallelized across chromosomes)
- **Validation task**: 1 core, 2GB RAM

### Chromosome Selection

The `makewindows` task allows you to specify which chromosomes to analyze:
- **Default**: All autosomes plus X and Y (chr1-chr22, chrX, chrY)
- **Custom**: Specify any subset of chromosomes based on your analysis needs
- **Performance**: Fewer chromosomes = faster runtime

### Intersect Flags

Customize BEDTools intersect behavior with the `intersect_flags` parameter:
- **"-header -wo"**: Include header and write original entries plus overlap (default)
- **"-wa -wb"**: Write both original A and B entries
- **"-c"**: Count overlaps for each A entry
- **"-v"**: Report A entries with no overlap in B

### Window Analysis

The `makewindows` task creates 500kb windows and counts reads with quality filters:
- **Window size**: Fixed at 500kb for consistency
- **Quality filters**: MAPQ ≥ 20, properly paired reads only
- **Output**: Per-chromosome BED files packed in a tarball

## Advanced Usage

### Custom Intersect Analysis

```wdl
call bedtools_tasks.intersect {
  input:
    bed_file = target_regions,
    aligned_bam = sample.bam,
    sample_name = sample.name,
    flags = "-wa -wb -f 0.5",  # Require 50% overlap
    cpu_cores = 4,
    memory_gb = 8
}
```

### Targeted Coverage Analysis

```wdl
call bedtools_tasks.coverage {
  input:
    bed_file = exon_regions,
    aligned_bam = sample.bam,
    sample_name = sample.name,
    cpu_cores = 2,
    memory_gb = 16
}
```

### Chromosome-Specific Window Analysis

```wdl
call bedtools_tasks.makewindows {
  input:
    bed_file = whole_genome_bed,
    aligned_bam = sample.bam,
    bam_index = sample.bam_index,
    reference_fasta = ref_fasta,
    reference_index = ref_fasta_index,
    sample_name = sample.name,
    list_chr = ["chr1", "chr2", "chr3"],  # Focus on specific chromosomes
    cpu_cores = 6,
    memory_gb = 20
}
```

## Output Format

### Coverage Output
Tab-delimited file with BED intervals plus mean coverage column.

### Intersect Output  
Format depends on flags used, typically includes original BED entries with alignment overlaps.

### Window Counts
Tarball containing per-chromosome BED files with 500kb windows and read counts.

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Input BAM files must be sorted and indexed
- Reference genome with FASTA index for window analysis
- BED file with genomic intervals of interest
- Docker/Apptainer support
- Sufficient computational resources (varies by task)

## Features

- **Multiple analysis types**: Coverage, intersection, and window-based counting
- **Automatic test data**: Downloads all required test data automatically
- **Parallel processing**: Chromosomes processed in parallel for efficiency
- **Quality filtering**: Built-in MAPQ and pairing filters for reliable results
- **Comprehensive validation**: Detailed output validation with file size and line counts
- **Flexible configuration**: Customizable parameters for all major settings
- **Integration ready**: Designed for use with other WILDS modules

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (chromosome 1 subset for efficiency)
- Comprehensive validation of all outputs
- Cross-platform compatibility testing

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
