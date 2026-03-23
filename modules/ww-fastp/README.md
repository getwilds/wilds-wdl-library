# ww-fastp Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for [fastp](https://github.com/OpenGene/fastp), an ultra-fast all-in-one FASTQ preprocessor. fastp performs quality filtering, adapter trimming, and comprehensive QC reporting in a single pass over the data.

## Overview

This module provides WDL tasks for running fastp on both paired-end and single-end sequencing data. fastp is significantly faster than similar tools (Trimmomatic, Cutadapt) while performing multiple preprocessing operations simultaneously, including adapter auto-detection, quality filtering, length filtering, and HTML/JSON report generation.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-fastp.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `fastp_paired`

Run fastp on paired-end FASTQ files for quality filtering, adapter trimming, and QC reporting.

**Inputs:**
- `sample_name` (String): Name identifier for the sample
- `r1_fastq` (File): Read 1 input FASTQ file (gzipped or uncompressed)
- `r2_fastq` (File): Read 2 input FASTQ file (gzipped or uncompressed)
- `qualified_quality_phred` (Int, default=15): Minimum base quality score for a base to be qualified
- `length_required` (Int, default=15): Minimum read length after trimming
- `detect_adapter_for_pe` (Boolean, default=true): Enable auto-detection of adapters for paired-end data
- `adapter_fasta` (File?, optional): FASTA file with custom adapter sequences
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `r1_trimmed` (File): Trimmed and filtered R1 FASTQ file
- `r2_trimmed` (File): Trimmed and filtered R2 FASTQ file
- `html_report` (File): HTML quality control report
- `json_report` (File): JSON quality control report

### `fastp_single`

Run fastp on single-end FASTQ files for quality filtering, adapter trimming, and QC reporting.

**Inputs:**
- `sample_name` (String): Name identifier for the sample
- `fastq` (File): Input FASTQ file (gzipped or uncompressed)
- `qualified_quality_phred` (Int, default=15): Minimum base quality score for a base to be qualified
- `length_required` (Int, default=15): Minimum read length after trimming
- `adapter_fasta` (File?, optional): FASTA file with custom adapter sequences
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `trimmed_fastq` (File): Trimmed and filtered FASTQ file
- `html_report` (File): HTML quality control report
- `json_report` (File): JSON quality control report

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-fastp/ww-fastp.wdl" as fastp_tasks

struct FastpSample {
    String name
    File r1_fastq
    File r2_fastq
}

workflow my_preprocessing_pipeline {
  input {
    Array[FastpSample] samples
  }

  scatter (sample in samples) {
    call fastp_tasks.fastp_paired {
      input:
        sample_name = sample.name,
        r1_fastq = sample.r1_fastq,
        r2_fastq = sample.r2_fastq
    }
  }

  output {
    Array[File] trimmed_r1 = fastp_paired.r1_trimmed
    Array[File] trimmed_r2 = fastp_paired.r2_trimmed
    Array[File] reports = fastp_paired.html_report
  }
}
```

### Advanced Usage Examples

**Custom quality filtering:**
```wdl
call fastp_tasks.fastp_paired {
  input:
    sample_name = "stringent_sample",
    r1_fastq = r1_file,
    r2_fastq = r2_file,
    qualified_quality_phred = 20,
    length_required = 50,
    cpu_cores = 8,
    memory_gb = 16
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-fastqc**: Can be used alongside fastp for additional QC checks
- **ww-bwa / ww-star**: Trimmed reads can be passed directly to alignment modules

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint fastp_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads paired-end test FASTQ data using `ww-testdata`
2. Runs fastp paired-end trimming and quality filtering
3. Runs fastp single-end trimming and quality filtering (using R1 as input)
4. Produces trimmed FASTQ files and HTML/JSON QC reports for both modes

## Docker Container

This module uses the `getwilds/fastp:1.1.0` container image, which includes:
- fastp v1.1.0
- All necessary system dependencies for FASTQ preprocessing

## Citation

> Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor.
> Bioinformatics. 2018;34(17):i884-i890.
> DOI: https://doi.org/10.1093/bioinformatics/bty560

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 4 cores
- **Memory**: 8 GB
- **Runtime**: Typically under 5 minutes for standard WGS samples

### Resource Scaling
- `cpu_cores`: fastp supports multi-threading; 4 cores is a good default
- `memory_gb`: 8 GB is sufficient for most use cases; increase for very large files

## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact us via the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for detailed guidelines

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
