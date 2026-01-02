# ww-sra
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for downloading genomic data from the NCBI Sequence Read Archive (SRA).

## Overview

This module provides reusable WDL tasks for downloading sequencing data from the SRA using the SRA toolkit. It handles both single-end and paired-end reads, automatically detecting the read type and processing accordingly.

The module uses `parallel-fastq-dump` for efficient, multi-threaded downloading of FASTQ files from SRA accessions.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `fastqdump`
- **Test workflow**: `testrun.wdl` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/sra-tools:3.1.1`

## Tasks

### `fastqdump`
Downloads FASTQ files from SRA accessions with automatic paired-end detection.

**Inputs:**
- `sra_id` (String): SRA accession ID
- `ncpu` (Int): Number of CPUs for parallel download (default: 8)
- `max_reads` (Int, optional): Maximum number of reads to download for testing/downsampling

**Outputs:**
- `r1_end` (File): R1 FASTQ file
- `r2_end` (File): R2 FASTQ file (empty for single-end)
- `is_paired_end` (Boolean): Paired-end detection flag

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl" as sra_tasks

workflow my_workflow {
  input {
    Array[String] sra_accessions = ["SRR12345678"]
  }

  scatter (id in sra_accessions) {
    call sra_tasks.fastqdump {
      input:
        sra_id = id,
        ncpu = 4
    }
  }

  output {
    Array[File] r1_fastqs = fastqdump.r1_end
    Array[File] r2_fastqs = fastqdump.r2_end
    Array[Boolean] is_paired_end = fastqdump.is_paired_end
  }
}
```

### Example with Downsampling

For testing or memory-constrained environments, limit the number of reads downloaded:

```wdl
call sra_tasks.fastqdump {
  input:
    sra_id = "SRR1039508",
    ncpu = 2,
    max_reads = 1000000  # Download only first 1M reads
}
```

### Integration Examples

This module pairs well with other WILDS modules:
- **ww-star**: For RNA-seq alignment after SRA download
- **ww-bwa**: For DNA-seq alignment after SRA download
- **ww-sra-star pipeline**: Complete pipeline from SRA to alignment

## Testing the Module

The `sra_example` test workflow in `testrun.wdl` requires no input parameters and automatically downloads data with hardcoded settings:

- **SRA ID**: ERR1258306 (for efficient testing)
- **CPU cores**: 2 (balanced performance)

### Running the Test Workflow

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint sra_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:
- **Memory**: Automatically scales with CPU count (2 * ncpu + " GB")
- **CPUs**: Adjustable based on available resources and download speed requirements
- **Network**: Requires stable internet connection for SRA downloads

### SRA Accession Types

- **SRR**: Single sample run accessions
- **ERR**: ENA (European Nucleotide Archive) accessions
- **DRR**: DDBJ (DNA Data Bank of Japan) accessions
- Supports both single-end and paired-end automatically

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for SRA downloads
- Sufficient storage space for downloaded FASTQ files

## Features

- **Automatic paired-end detection**: Determines read structure automatically
- **Parallel downloading**: Multi-threaded downloads for improved performance
- **Standardized output**: Consistent naming for downstream processing
- **Cross-platform**: Works with SRA accessions from NCBI, ENA, and DDBJ
- **Optional downsampling**: Limit read count for testing and development via `max_reads` parameter

## Performance Considerations

- **Download speed**: Performance scales with CPU count and network bandwidth
- **Storage requirements**: FASTQ files can be large; ensure adequate disk space
- **Memory usage**: Scales automatically with CPU allocation
- **Downsampling**: Use `max_reads` parameter to limit data download for testing or when working with constrained resources (e.g., CI runners, local development)

## Output Description

- **R1 FASTQ files**: Forward reads in gzip-compressed FASTQ format
- **R2 FASTQ files**: Reverse reads for paired-end data (empty file for single-end)
- **Paired-end flags**: Boolean indicators for read structure per sample

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real SRA accessions for integration testing
- Comprehensive validation of all outputs including FASTQ format validation
- Tests with both single-end and paired-end data

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to SRA toolkit usage, please refer to the [NCBI SRA toolkit documentation](https://github.com/ncbi/sra-tools).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.

