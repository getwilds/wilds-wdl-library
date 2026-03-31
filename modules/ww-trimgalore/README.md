# ww-trimgalore Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module wrapping [Trim Galore](https://github.com/FelixKrueger/TrimGalore), a wrapper around Cutadapt and FastQC for consistent adapter and quality trimming of FASTQ files. Trim Galore auto-detects adapter sequences (Illumina, Nextera, Small RNA, and more) and performs quality trimming in a single step.

## Overview

This module provides tasks for running Trim Galore on both paired-end and single-end FASTQ data. Trim Galore automates adapter detection and quality trimming using Cutadapt under the hood, with optional FastQC reporting after trimming.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-trimgalore.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `trimgalore_paired`

Run Trim Galore on paired-end FASTQ files for adapter and quality trimming.

**Inputs:**
- `sample_name` (String): Name identifier for the sample
- `r1_fastq` (File): Read 1 input FASTQ file (gzipped or uncompressed)
- `r2_fastq` (File): Read 2 input FASTQ file (gzipped or uncompressed)
- `quality` (Int, default=20): Phred quality score threshold for trimming low-quality ends
- `length` (Int, default=20): Minimum read length after trimming; shorter reads are discarded
- `stringency` (Int, default=1): Minimum overlap with adapter sequence required to trim
- `run_fastqc` (Boolean, default=false): Run FastQC on trimmed files
- `adapter` (String?, optional): Adapter sequence for R1 (auto-detected if not specified)
- `adapter2` (String?, optional): Adapter sequence for R2 (auto-detected if not specified)
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `r1_trimmed` (File): Trimmed and validated R1 FASTQ file
- `r2_trimmed` (File): Trimmed and validated R2 FASTQ file
- `r1_report` (File): Trimming report for R1
- `r2_report` (File): Trimming report for R2

### `trimgalore_single`

Run Trim Galore on a single-end FASTQ file for adapter and quality trimming.

**Inputs:**
- `sample_name` (String): Name identifier for the sample
- `fastq` (File): Input FASTQ file (gzipped or uncompressed)
- `quality` (Int, default=20): Phred quality score threshold for trimming low-quality ends
- `length` (Int, default=20): Minimum read length after trimming; shorter reads are discarded
- `stringency` (Int, default=1): Minimum overlap with adapter sequence required to trim
- `run_fastqc` (Boolean, default=false): Run FastQC on trimmed files
- `adapter` (String?, optional): Adapter sequence (auto-detected if not specified)
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `trimmed_fastq` (File): Trimmed FASTQ file
- `trimming_report` (File): Trimming report

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-trimgalore/ww-trimgalore.wdl" as trimgalore_tasks

struct TrimGaloreSample {
    String name
    File r1_fastq
    File r2_fastq
}

workflow my_trimming_pipeline {
  input {
    Array[TrimGaloreSample] samples
  }

  scatter (sample in samples) {
    call trimgalore_tasks.trimgalore_paired {
      input:
        sample_name = sample.name,
        r1_fastq = sample.r1_fastq,
        r2_fastq = sample.r2_fastq
    }
  }

  output {
    Array[File] r1_trimmed = trimgalore_paired.r1_trimmed
    Array[File] r2_trimmed = trimgalore_paired.r2_trimmed
  }
}
```

### Advanced Usage Examples

**Custom quality threshold and adapter:**
```wdl
call trimgalore_tasks.trimgalore_paired {
  input:
    sample_name = "stringent_sample",
    r1_fastq = r1_file,
    r2_fastq = r2_file,
    quality = 30,
    length = 50,
    run_fastqc = true,
    cpu_cores = 8,
    memory_gb = 16
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-salmon**: Quantify RNA transcripts from trimmed FASTQs
- **ww-star / ww-bwa / ww-bowtie / ww-bowtie2**: Use trimmed FASTQ output as input to alignment modules
- **ww-fastqc**: Additional standalone QC if needed beyond Trim Galore's built-in FastQC option

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint trimgalore_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads test FASTQ data using `ww-testdata`
2. Runs both paired-end and single-end trimming on the test data

## Docker Container

This module uses the `getwilds/trim-galore:0.6.11` container image, which includes:
- Trim Galore v0.6.11
- Cutadapt (adapter trimming engine)
- FastQC (optional QC reporting)
- All necessary Perl and system dependencies

## Citation

> Trim Galore
> Felix Krueger
> https://github.com/FelixKrueger/TrimGalore

> Cutadapt: Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12. DOI: 10.14806/ej.17.1.200

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 4 cores
- **Memory**: 8 GB
- **Runtime**: A few minutes per sample for typical FASTQ files

### Resource Scaling
- `cpu_cores`: Trim Galore supports multi-core via `--cores`; 4 cores is a good default
- `memory_gb`: 8 GB is sufficient for most datasets; increase for very large files
- Paired-end mode processes both reads simultaneously and benefits from additional cores

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
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for detailed guidelines

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
