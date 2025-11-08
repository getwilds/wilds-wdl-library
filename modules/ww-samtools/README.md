# ww-samtools
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for processing genomic files with [Samtools](http://www.htslib.org/).

## Overview

This module provides reusable WDL tasks for processing genomic data with **Samtools**, including converting CRAM/BAM/SAM files to FASTQ format and merging BAM files to CRAM format. Designed to be a modular component in the WILDS ecosystem, this module is suitable for integration into larger bioinformatics pipelines and is automatically validated with real sequencing data via its test workflow.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `crams_to_fastq`, `merge_bams_to_cram`
- **Test workflow**: `testrun.wdl` (demonstration workflow executing all tasks)
- **Container**: `getwilds/samtools:1.19`

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

### `merge_bams_to_cram`

Merges multiple BAM files into a single queryname-sorted CRAM file using samtools merge and collate. The output CRAM is sorted by query name (read name), which is appropriate for unmapped reads and does not support indexing.

**Inputs:**
- `bams_to_merge` (Array[File]): Array of BAM files to merge into a single CRAM file
- `base_file_name` (String): Base name for output CRAM file
- `cpu_cores` (Int): Number of CPU cores to use (threads = cpu_cores - 1) (default: 6)
- `memory_gb` (Int): Memory allocation in GB (default: 12)

**Outputs:**
- `cram` (File): Merged queryname-sorted CRAM file containing all reads from input BAMs

**Note:** The output CRAM is sorted by query name and does not have an index file. CRAM indexes are only applicable to coordinate-sorted files.

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

  output {
    Array[File] r1_fastqs = crams_to_fastq.r1_fastq
    Array[File] r2_fastqs = crams_to_fastq.r2_fastq
  }
}
```

### Integration Examples

This module pairs well with other WILDS modules:
- **ww-sra**: For downloading genomic data
- **ww-bwa**: For FASTQ sequence alignment.

## Testing the Module

The module includes a test workflow that automatically downloads test data and runs without requiring input files:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint samtools_example
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Sufficient computational resources (can be memory-intensive for large genomes)

## Features

- **Parallel processing**: Multi-threaded execution for improved performance
- **Standardized output**: Uncompressed FASTQ format compatible with downstream tools

## Performance Considerations

- **Memory usage**: Whole-genome analysis typically requires 16-32GB RAM
- **CPU scaling**: Performance improves with additional cores (recommend 8-16 CPUs)
- **Region exclusion**: Using exclude regions BED file significantly improves runtime and accuracy
- **SV type filtering**: Calling specific SV types reduces runtime compared to calling all types

## Output Description
- **FASTQ files**: FASTQ output files for each sample

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (chromosome 1 subset for efficiency)
- Comprehensive validation of all outputs and statistics

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to Samtools usage or configuration, please refer to the documentation present in the [Samtools website](http://www.htslib.org/). Please make sure to cite their work if you use Samtools in your analyses:

Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PMID: 19505943; PMCID: PMC2723002.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.