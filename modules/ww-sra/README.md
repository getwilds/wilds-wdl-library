# ww-sra
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for downloading genomic data from the NCBI Sequence Read Archive (SRA).

## Overview

This module provides reusable WDL tasks for downloading sequencing data from the SRA using the SRA toolkit. It handles both single-end and paired-end reads, automatically detecting the read type and processing accordingly. The module includes built-in validation to ensure successful downloads.

The module uses `parallel-fastq-dump` for efficient, multi-threaded downloading of FASTQ files from SRA accessions.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `fastqdump`, `validate_outputs`
- **Workflow**: `sra_download` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/sra-tools:3.1.1`

## Tasks

### `fastqdump`
Downloads FASTQ files from SRA accessions with automatic paired-end detection.

**Inputs:**
- `sra_id` (String): SRA accession ID
- `ncpu` (Int): Number of CPUs for parallel download

**Outputs:**
- `r1_end` (File): R1 FASTQ file
- `r2_end` (File): R2 FASTQ file (empty for single-end)
- `is_paired_end` (Boolean): Paired-end detection flag

### `validate_outputs`
Validates downloaded files and generates a report.

**Inputs:**
- Arrays of files and metadata from `fastqdump`

**Outputs:**
- `report` (File): Validation summary

## Usage as a Module

### Importing into Your Workflow

```wdl
import "path/to/ww-sra.wdl" as sra_tasks

workflow my_workflow {
  input {
    String sra_accession = "SRR12345678"
  }
  
  call sra_tasks.fastqdump {
    input: 
      sra_id = sra_accession,
      ncpu = 4
  }
  
  output {
    File r1_fastq = fastqdump.r1_end
    File r2_fastq = fastqdump.r2_end
  }
}
```

### Integration Examples

This module pairs well with other WILDS modules:
- **ww-star**: For RNA-seq alignment after SRA download
- **ww-sra-star vignette**: Complete pipeline from SRA to alignment

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-sra.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-sra.wdl -i inputs.json

# Using Sprocket
sprocket run ww-sra.wdl inputs.json
```

### Test Input Format

```json
{
  "sra_download.sra_id_list": ["SRR13191702"],
  "sra_download.n_cpu": 1
}
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for SRA downloads

## Features

- **Automatic paired-end detection**: Determines read structure automatically
- **Parallel downloading**: Multi-threaded downloads for improved performance
- **Validation**: Built-in file validation and reporting
- **Standardized output**: Consistent naming for downstream processing
- **Error handling**: Robust error detection and reporting

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline. Tests run on multiple WDL executors to ensure compatibility and reliability.

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.

