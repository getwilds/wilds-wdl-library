# ww-sra
[![Project Status: Stable – Useable, full support, open to feedback, stable API.](https://getwilds.org/badges/badges/stable.svg)](https://getwilds.org/badges/#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for downloading genomic data from the NCBI Sequence Read Archive (SRA).

## Overview

This module provides reusable WDL tasks for downloading sequencing data from the SRA using the SRA toolkit. It handles both single-end and paired-end reads, automatically detecting the read type and processing accordingly.

The module uses `prefetch` + `fasterq-dump` for efficient, multi-threaded downloading of FASTQ files from SRA accessions, with optional NGC authentication for downloading controlled-access dbGaP data.

## Data Governance

When pulling controlled-access data from dbGaP via the `ngc_file` input, you are bound by the data use agreement attached to that study, which typically restricts where the downloaded FASTQs (and any data derived from them) may be stored. **Make sure to run any workflow that imports this module in a location approved for regulated data storage** and avoid persisting outputs to general-purpose or shared filesystems.

**Fred Hutch users:** Use [PROOF Regulated](https://sciwiki.fredhutch.org/datademos/proof-regulated/) to submit workflows that use this module. PROOF Regulated stages analysis data under `/fh/regulated`, which is set up for compliance with dbGaP and similar data use agreements.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `fastqdump`
- **Test workflow**: `testrun.wdl` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/sra-tools:3.1.1`

## Tasks

### `fastqdump`
Downloads FASTQ files from SRA accessions with automatic read structure detection. Handles paired-end data including datasets with index fastqs (e.g., 10x Chromium with barcode, index, and cDNA reads) by including technical reads and then classifying them by length so that R1 and R2 outputs always contain the biological reads (see [Output File Classification](#output-file-classification) below). Supports controlled-access dbGaP data via optional NGC repository key file.

**Inputs:**
- `sra_id` (String): SRA accession ID
- `ncpu` (Int): Number of CPUs for parallel download (default: 8)
- `max_reads` (Int, optional): Maximum number of reads to download for testing/downsampling
- `ngc_file` (File, optional): NGC repository key file for downloading controlled-access dbGaP data
- `docker_image` (String): Docker image to use for this task (default: `getwilds/sra-tools:3.1.1`)

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

### Example with dbGaP Controlled-Access Data

To download controlled-access data from dbGaP, provide your NGC repository key file:

```wdl
call sra_tasks.fastqdump {
  input:
    sra_id = "SRR12345678",
    ncpu = 4,
    ngc_file = ngc_key  # Your dbGaP NGC repository key file
}
```

> **Note:** Controlled-access dbGaP data must be handled in a regulated computing environment. At Fred Hutch, this means using [PROOF Regulated](https://sciwiki.fredhutch.org/datademos/proof-regulated/) to ensure compliance with data use agreements and institutional security requirements.

### Integration Examples

This module pairs well with other WILDS modules:
- **ww-star**: For RNA-seq alignment after SRA download
- **ww-bwa**: For DNA-seq alignment after SRA download
- **ww-sra-star pipeline**: Complete pipeline from SRA to alignment
- **ww-sra-cellranger pipeline**: Complete pipeline from SRA to single-cell gene expression quantification

## Testing the Module

The `sra_example` test workflow in `testrun.wdl` requires no input parameters and automatically downloads data with hardcoded settings:

- **SRA ID**: ERR1258306 (for efficient testing)
- **CPU cores**: 2 (balanced performance)

### Running the Test Workflow

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl

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

- **Automatic read structure detection**: Handles single-end and paired-end data automatically
- **Parallel downloading**: Multi-threaded downloads for improved performance
- **Standardized output**: Consistent naming for downstream processing
- **Cross-platform**: Works with SRA accessions from NCBI, ENA, and DDBJ
- **Optional downsampling**: Limit read count for testing and development via `max_reads` parameter
- **dbGaP support**: Download controlled-access data using NGC repository key files via optional `ngc_file` parameter

## Performance Considerations

- **Download speed**: Performance scales with CPU count and network bandwidth
- **Storage requirements**: FASTQ files can be large; ensure adequate disk space
- **Memory usage**: Scales automatically with CPU allocation
- **Downsampling**: Use `max_reads` parameter to limit data download for testing or when working with constrained resources (e.g., CI runners, local development)

## Output File Classification

`fasterq-dump --split-files --include-technical` can produce 2–4 output files per SRA accession depending on what the submitter uploaded. The file ordering is not consistent across submissions: a 10x Chromium dataset, for example, may land as `{R1=28bp, R2=98bp, I1=8bp}` in any file-index order. To return biological reads reliably, `fastqdump` inspects each output file and assigns roles by read length:

- Any file whose median read length is `<= 15bp` is treated as a technical/index read and dropped.
- Of the remaining biological-read files, the one with the longest median is `R2` and the next is `R1`.
- Files of equal median length break ties by original file index, preserving the conventional `R1=_1, R2=_2` mapping for plain paired-end submissions.
- If only one biological-read file survives, the task treats the sample as single-end (overriding `is_paired_end` to `false`) and writes an empty `R2`.

The classification table is echoed to stdout for every run, so the file→role mapping is auditable from task logs.

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

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to SRA toolkit usage, please refer to the [NCBI SRA toolkit documentation](https://github.com/ncbi/sra-tools).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.

