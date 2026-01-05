# ww-megahit
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for de novo metagenomic assembly using [MEGAHIT](https://github.com/voutcn/megahit).

## Overview

This module provides a reusable WDL task for performing fast and memory-efficient de novo metagenomic assembly with MEGAHIT.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `megahit`
- **Test workflow**: `testrun.wdl` (demonstration workflow executing the task)
- **Container**: `getwilds/megahit:1.2.9` (WILDS Docker image with MEGAHIT installed)

## Tasks

### `megahit`

Performs de novo metagenomic assembly using MEGAHIT.

**Inputs:**
- `input_fastq` (File): Interleaved FASTQ file of paired-end reads
- `sample_name` (String): Sample name for output file naming (default: basename of input FASTQ)
- `cpu_cores` (Int): Number of CPU cores allocated for the task (default: 2)
- `memory_gb` (Int): Memory allocated for the task in GB (default: 4)

**Outputs:**
- `contigs` (File): Assembled contigs in compressed FASTA format (`.contigs.fasta.gz`)

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-megahit/ww-megahit.wdl" as megahit_tasks

workflow my_assembly_pipeline {
  input {
    File interleaved_reads
    String sample_id
  }

  call megahit_tasks.megahit {
    input:
      input_fastq = interleaved_reads,
      sample_name = sample_id,
      cpu_cores = 4,
      memory_gb = 8
  }

  output {
    File assembled_contigs = megahit.contigs
  }
}
```

### Integration Examples

This module pairs well with other WILDS modules:
- **ww-sra**: For downloading FASTQ files from SRA
- **ww-testdata**: For interleaving paired-end FASTQs

## Testing the Module

The module includes a test workflow that automatically downloads test data and runs without requiring input files:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

The test workflow (`megahit_example`) automatically:
1. Downloads metagenomic FASTQ data from SRA using `ww-sra`
2. Interleaves paired-end reads using `ww-testdata`
3. Performs assembly using MEGAHIT
4. Validates all outputs

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Sufficient computational resources (MEGAHIT is optimized for modest memory requirements)

## Performance Considerations

- **Memory usage**: MEGAHIT is designed to be memory-efficient; typically requires 4-16GB RAM for most datasets
- **CPU scaling**: Performance improves with additional cores (recommend 4-8 CPUs)
- **Disk space**: Ensure adequate temporary disk space (MEGAHIT intermediate files are cleaned up automatically)
- **Runtime**: Assembly time varies based on dataset complexity and size

## Output Description

- **Contigs**: Assembled contiguous sequences in FASTA format, compressed with gzip. These represent the assembled metagenomic sequences.

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data
- Comprehensive validation of all outputs

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to MEGAHIT usage or configuration, please refer to the [MEGAHIT documentation](https://github.com/voutcn/megahit). Please make sure to cite their work if you use MEGAHIT in your analyses:

```
Li D, Liu CM, Luo R, Sadakane K, Lam TW. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics. 2015 May 15;31(10):1674-6. doi: 10.1093/bioinformatics/btv033. Epub 2015 Jan 20. PMID: 25609793.
```

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
