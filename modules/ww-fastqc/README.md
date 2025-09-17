# ww-fastqc Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for FastQC quality control analysis of high-throughput sequencing data. FastQC provides comprehensive quality control checks on raw sequence data coming from high-throughput sequencing pipelines.

## Overview

This module wraps the FastQC bioinformatics tool to generate quality control reports for FASTQ files. It supports both single-end and paired-end sequencing data and produces both HTML reports for visualization and ZIP archives containing all analysis data.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-fastqc.wdl` - Contains all tasks and demonstration workflow
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `run_fastqc`

Run FastQC quality control analysis on FASTQ files.

**Inputs:**
- `sample_name` (String): Name identifier for the sample
- `r1_fastq` (File, optional): Read 1 FASTQ file (required for analysis)
- `r2_fastq` (File, optional): Read 2 FASTQ file (optional for paired-end data)
- `cpu_cores` (Int, default=2): Number of CPU cores allocated for FastQC
- `memory_gb` (Int, default=4): Memory allocated for FastQC in GB
- `adapters` (File, optional): Adapter sequences file for contamination screening
- `limits` (File, optional): Limits file to override default warning/error thresholds
- `contaminants` (File, optional): Contaminants file for contamination screening

**Outputs:**
- `html_reports` (Array[File]): FastQC HTML quality control reports
- `zip_reports` (Array[File]): FastQC ZIP archives containing all report data

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-fastqc/ww-fastqc.wdl" as fastqc_tasks

struct FastQCSample {
    String name
    File? r1_fastq
    File? r2_fastq
}

workflow my_analysis_pipeline {
  input {
    Array[FastQCSample] samples
  }

  scatter (sample in samples) {
    call fastqc_tasks.run_fastqc {
      input:
        sample_name = sample.name,
        r1_fastq = sample.r1_fastq,
        r2_fastq = sample.r2_fastq
    }
  }

  output {
    Array[Array[File]] fastqc_html_reports = run_fastqc.html_reports
    Array[Array[File]] fastqc_zip_reports = run_fastqc.zip_reports
  }
}
```

### Advanced Usage Examples

**Custom resource allocation:**
```wdl
call fastqc_tasks.run_fastqc {
  input:
    sample_name = "large_sample",
    r1_fastq = large_r1_file,
    r2_fastq = large_r2_file,
    cpu_cores = 4,
    memory_gb = 8
}
```

**With custom adapter file:**
```wdl
call fastqc_tasks.run_fastqc {
  input:
    sample_name = "custom_adapters_sample",
    r1_fastq = input_r1,
    r2_fastq = input_r2,
    adapters = custom_adapters_file
}
```

**Single-end data:**
```wdl
call fastqc_tasks.run_fastqc {
  input:
    sample_name = "single_end_sample",
    r1_fastq = single_end_fastq
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of test data for demonstrations
- **Quality control pipelines**: Can be used as an initial QC step before alignment or assembly
- **Other WILDS modules**: Commonly used before alignment (ww-bwa, ww-bowtie2) or trimming tools

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-fastqc.wdl

# Using miniWDL
miniwdl run ww-fastqc.wdl

# Using Sprocket
sprocket run ww-fastqc.wdl
```

### Automatic Demo Mode

The workflow automatically:
1. Downloads test FASTQ data using `ww-testdata`
2. Runs FastQC analysis on both paired-end and single-end examples
3. Generates HTML and ZIP reports for quality assessment

## Docker Container

This module uses the `getwilds/fastqc:0.12.1` container image, which includes:
- FastQC v0.12.1
- Java runtime environment
- All necessary system dependencies for quality control analysis

## FastQC Information

- **Purpose**: Quality control analysis for high-throughput sequencing data
- **Tool**: FastQC v0.12.1
- **Output formats**: HTML reports (human-readable) and ZIP archives (machine-readable)
- **Demo Data**: Uses test FASTQ files from ww-testdata module

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 2 cores
- **Memory**: 4 GB
- **Runtime**: Typically 1-5 minutes per FASTQ file depending on file size

### Resource Scaling
FastQC resource requirements scale with input file size:
- `cpu_cores`: FastQC can utilize multiple threads for faster processing
- `memory_gb`: Increase for very large FASTQ files (>10GB)
- Small files (< 1GB): Default resources are sufficient
- Large files (> 5GB): Consider increasing to 4 cores and 8GB memory

## FastQC Quality Metrics

FastQC analyzes multiple quality aspects:
- **Per base sequence quality**: Quality scores across read positions
- **Per sequence quality scores**: Overall quality distribution
- **Per base sequence content**: Nucleotide composition bias
- **Per sequence GC content**: GC content distribution
- **Per base N content**: Ambiguous base calls
- **Sequence length distribution**: Read length consistency
- **Sequence duplication levels**: PCR duplicate detection
- **Overrepresented sequences**: Adapter contamination and artifacts
- **Adapter content**: Adapter sequence contamination
- **Kmer content**: Overrepresented k-mer sequences

## Common Use Cases

### Pre-alignment Quality Control
```wdl
# Run FastQC before alignment to assess raw data quality
call fastqc_tasks.run_fastqc { input: sample_name = "raw_qc", r1_fastq = raw_r1, r2_fastq = raw_r2 }
```

### Post-trimming Quality Assessment
```wdl
# Compare quality before and after read trimming
call fastqc_tasks.run_fastqc { input: sample_name = "post_trim_qc", r1_fastq = trimmed_r1, r2_fastq = trimmed_r2 }
```

### Batch Processing
```wdl
# Process multiple samples in parallel
scatter (sample in samples) {
  call fastqc_tasks.run_fastqc { input: sample_name = sample.name, r1_fastq = sample.r1, r2_fastq = sample.r2 }
}
```

## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**: Official FastQC tool documentation
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
