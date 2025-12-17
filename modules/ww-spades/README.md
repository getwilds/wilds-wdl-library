# ww-spades
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for de novo metagenomic assembly using [SPAdes](https://github.com/ablab/spades).

## Overview

This module provides a reusable WDL task for performing de novo metagenomic assembly with **metaSPAdes**. The module runs the assembler only (without read error correction) and requires a single paired-end library as input.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `metaspades`
- **Test workflow**: `testrun.wdl` (demonstration workflow executing the task)
- **Container**: `staphb/spades:4.2.0`

## Tasks

### `metaspades`

Performs de novo metagenomic assembly using metaSPAdes in assembler-only mode (skips read error correction).

**Inputs:**
- `r1_fastq` (File, optional): Read 1 FASTQ file
- `r2_fastq` (File, optional): Read 2 FASTQ file
- `interleaved_fastq` (File, optional): Interleaved FASTQ file of paired-end reads
- `sample_name` (String): Sample name for output file naming
- `cpu_cores` (Int): Number of CPU cores allocated for the task (default: 16)
- `memory_gb` (Int): Memory allocated for the task in GB (default: 64)

**Note**: You must provide either `interleaved_fastq` OR both `r1_fastq` and `r2_fastq`. Input must be from a single, paired-end library.

**Outputs:**
- `scaffolds` (File): Assembled scaffolds in compressed FASTA format (`.scaffolds.fasta.gz`)
- `contigs` (File): Assembled contigs in compressed FASTA format (`.contigs.fasta.gz`)
- `log` (File): SPAdes log file

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-spades/ww-spades.wdl" as spades_tasks

workflow my_assembly_pipeline {
  input {
    File read1_fastq
    File read2_fastq
    String sample_id
  }

  call spades_tasks.metaspades {
    input:
      r1_fastq = read1_fastq,
      r2_fastq = read2_fastq,
      sample_name = sample_id,
      cpu_cores = 4,
      memory_gb = 8
  }

  output {
    File assembled_scaffolds = metaspades.scaffolds_fasta
    File assembled_contigs = metaspades.contigs_fasta
    File spades_log = metaspades.log
  }
}
```

### Integration Examples

This module pairs well with other WILDS modules:
- **ww-sra**: For downloading metagenomic data from SRA

## Testing the Module

The module includes a test workflow that automatically downloads test data and runs without requiring input files:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint spades_example
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Sufficient computational resources (metagenomic assembly is computationally intensive)

## Features

- **Flexible input formats**: Accepts either separate R1/R2 files or interleaved FASTQ
- **Compressed output**: Scaffolds are automatically gzip-compressed to save storage space

## Performance Considerations

- **Memory usage**: Metagenomic assembly is memory-intensive; recommend 64GB+ RAM for typical datasets
- **CPU scaling**: Performance improves with additional cores (recommend 16-32 CPUs)
- **Disk space**: Ensure adequate temporary disk space (SPAdes intermediate files are cleaned up automatically)
- **Runtime**: Assembly time varies significantly based on dataset complexity and size

## Output Description

- **Scaffolds**: Assembled contiguous sequences (scaffolds) in FASTA format, compressed with gzip. Reecommend for use as resulting sequences.
- **Contigs**: Assembled contiguous sequences (contigs) in FASTA format, compressed with gzip.
- **Log file**: Detailed log of the SPAdes assembly process for troubleshooting and record-keeping.

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data
- Comprehensive validation of all outputs

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to SPAdes usage or configuration, please refer to the [SPAdes manual](https://ablab.github.io/spades/). Please make sure to cite their work if you use metaSPAdes in your analyses:

Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. metaSPAdes: a new versatile metagenomic assembler. Genome Res. 2017 May;27(5):824-834. doi: 10.1101/gr.213959.116. Epub 2017 Mar 15. PMID: 28298430; PMCID: PMC5411777.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
