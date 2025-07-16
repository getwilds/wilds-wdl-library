# ww-samtools
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for processing genomic files with [Samtools](http://www.htslib.org/).

## Overview

This module provides reusable WDL tasks for merging and converting genomic data (CRAM/BAM/SAM) to FASTQ format using **Samtools**. It includes built-in validation to ensure that the resulting FASTQ files are successfully generated.

Designed to be a modular component in the WILDS ecosystem, this workflow is suitable for both standalone execution and integration into larger bioinformatics pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Task**: `crams_to_fastq`, `validate_outputs`
- **Workflow**: `samtools_example` (demonstration workflow executing all tasks)
- **Container**: `getwilds/samtools:1.11`

## Tasks

### `crams_to_fastq`

Merges one or more CRAM/BAM/SAM files for a sample, sorts by read name, and converts the result to a compressed FASTQ.

**Inputs:**
- `cram_files` (Array[String]): List of CRAM/BAM/SAM files to merge and convert
- `name` (String): Sample name used for output filenames
- `cpu_cores` (Int): Number of CPU cores to use (default: 23)
- `memory_gb` (Int): Memory allocation in GB (default: 36)

**Outputs:**
- `fastq_file` (File): FASTQ output file (`.fastq.gz`)
- `sample_name` (String): Sample name that was processed

### `validate_outputs`

Checks that all FASTQ files were created and are non-empty. Also generates a human-readable summary report.

**Inputs:**
- `fastq_files` (Array[File]): FASTQ files to validate
- `sample_names` (Array[String]): Names of corresponding samples

**Outputs:**
- `report` (File): Validation report with status and basic stats

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-samtools/ww-samtools.wdl" as samtools_tasks

struct SampleInfo {
    String name
    Array[String] cram_files
}

workflow my_preprocessing_pipeline {
  input {
    Array[SampleInfo] samples
  }

  scatter (sample in samples) {
    call samtools_tasks.crams_to_fastq {
      input:
        cram_files = sample.cram_files,
        name = sample.name
    }
  }

  call samtools_tasks.validate_outputs {
    input:
      fastq_files = crams_to_fastq.fastq_file,
      sample_names = crams_to_fastq.sample_name
  }

  output {
    Array[File] fastqs = crams_to_fastq.fastq_file
    File validation_report = validate_outputs.report
  }
}
```

### Integration Examples

This module pairs well with other WILDS modules:
- **ww-sra**: For downloading genomic data
- **ww-bwa**: For FASTQ sequence alignment.

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-samtools.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-samtools.wdl -i inputs.json

# Using Sprocket
sprocket run ww-samtools.wdl inputs.json
```

### Test Input Format

```json
{
  "samtools_example.samples": [
    {
      "name": "ERR1258306",
      "cram_files": ["path/to/first.cram", "pat/to/second.cram"]
    }
  ],
  "samtools_example.cpus": 2,
  "samtools_example.memory_gb": 4
}
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Sufficient computational resources (can be memory-intensive for large genomes)

## Features

- **Parallel processing**: Multi-threaded execution for improved performance
- **Quality validation**: Built-in output validation and statistics
- **Standardized output**: Uncompressed FASTQ format compatible with downstream tools
- **Detailed reporting**: Validation that output files exist and are non-empty

## Performance Considerations

- **Memory usage**: Whole-genome analysis typically requires 16-32GB RAM
- **CPU scaling**: Performance improves with additional cores (recommend 8-16 CPUs)
- **Region exclusion**: Using exclude regions BED file significantly improves runtime and accuracy
- **SV type filtering**: Calling specific SV types reduces runtime compared to calling all types

## Output Description
- **FASTQ files**: FASTQ output files for each sample
- **Validation report**: Validation report confirming all expected outputs were generated

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (chromosome 1 subset for efficiency)
- Comprehensive validation of all outputs and statistics

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to Delly usage or configuration, please refer to the documentation present in the [Delly GitHub repository](https://github.com/dellytools/delly). Please make sure to cite their work if you use Delly in your analyses:

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.
DELLY: structural variant discovery by integrated paired-end and split-read analysis.
Bioinformatics. 2012 Sep 15;28(18):i333-i339.
https://doi.org/10.1093/bioinformatics/bts378

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.