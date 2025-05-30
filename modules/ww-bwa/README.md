# ww-bwa
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)  
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for sequence alignment using the Burrows-Wheeler Aligner (BWA-MEM).

## Overview

This module provides reusable WDL tasks for aligning sequencing reads to a reference genome using **BWA-MEM**, a fast and accurate aligner suitable for reads ranging from 70bp to 1Mbp. The module includes a demonstration workflow that performs alignment and outputs sorted BAM files, and is designed to integrate seamlessly into broader pipelines such as germline or somatic variant calling.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Task**: `bwa_mem`
- **Workflow**: `bwa_example` (demonstration workflow for aligning paired-end reads)
- **Container**: `getwilds/bwa:0.7.17`

## Tasks

### `bwa_mem`

Aligns paired-end reads to a reference genome using BWA-MEM and outputs a sorted BAM file.

**Inputs:**
- `sample_data` (`SampleInfo`): Struct containing sample name and R1/R2 FASTQ files
- `reference_fasta` (`File`): Reference genome FASTA file (indexed and preprocessed externally)
- `cpu_cores` (`Int`): Number of CPU cores allocated (default: 8)
- `memory_gb` (`Int`): Memory allocated in GB (default: 62)

**Outputs:**
- `sorted_bam` (`File`): Sorted alignment file in BAM format

## Workflow

### `bwa_example`

A demonstration workflow that aligns an array of samples to a provided reference genome.

**Inputs:**
- `samples` (`Array[SampleInfo]`): Array of sample structs
- `reference_genome` (`RefGenome`): Struct containing the reference FASTA and associated index files
- `cpus` (`Int`): Number of CPUs per task (default: 8)
- `memory_gb` (`Int`): Memory per task in GB (default: 64)

**Outputs:**
- `bwa_bam` (`Array[File]`): Sorted BAM files for each sample

## Usage as a Module

### Importing into Your Workflow

```wdl
import "path/to/ww-bwa.wdl" as bwa_tasks

struct SampleInfo {
    String name
    File r1
    File r2
}

struct RefGenome {
    File ref_dict
    File ref_fasta
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_fasta_index
    File ref_pac
    File ref_sa
}

workflow my_alignment_pipeline {
  input {
    Array[SampleInfo] samples
    RefGenome reference_genome
  }

  scatter (sample in samples) {
    call bwa_tasks.bwa_mem {
      input:
        sample_data = sample,
        reference_fasta = reference_genome.ref_fasta
    }
  }

  output {
    Array[File] aligned_bams = bwa_mem.sorted_bam
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-sra**: Download sequencing data prior to alignment
- **ww-gatk**: Use the output BAMs for germline variant calling with GATK
- **ww-qc**: Validate and summarize BAM file quality

## Testing the Module

The module includes a demonstration workflow (`bwa_example`) with support for execution on multiple WDL backends:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-bwa.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-bwa.wdl -i inputs.json

# Using Sprocket
sprocket run ww-bwa.wdl inputs.json
```

### Test Input Format

```json
{
  "bwa_example.samples": [
    {
      "name": "sample1",
      "r1": "/path/to/sample1_R1.fastq.gz",
      "r2": "/path/to/sample1_R2.fastq.gz"
    }
  ],
  "bwa_example.reference_genome": {
    "ref_dict": "/path/to/genome.dict",
    "ref_fasta": "/path/to/genome.fasta",
    "ref_amb": "/path/to/genome.amb",
    "ref_ann": "/path/to/genome.ann",
    "ref_bwt": "/path/to/genome.bwt",
    "ref_fasta_index": "/path/to/genome.fasta.fai",
    "ref_pac": "/path/to/genome.pac",
    "ref_sa": "/path/to/genome.sa"
  },
  "bwa_example.cpus": 8,
  "bwa_example.memory_gb": 64
}
```

## Configuration Guidelines

### Resource Allocation

- **Memory**: 62–64 GB recommended for human genomes; can be tuned for smaller references
- **CPUs**: 8 cores typically sufficient for most samples; BWA-MEM and samtools scale with available threads

### Advanced Parameters

- Ensure the reference genome is pre-indexed with BWA and has an associated `.fai` and `.dict`
- Set `cpu_cores` and `memory_gb` based on your environment and dataset size

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker or Apptainer support
- Indexed reference genome files (e.g., `.bwt`, `.amb`, `.sa`, `.fai`, `.dict`)

## Features

- **BWA-MEM alignment**: Fast and accurate for reads 70bp to 1Mbp
- **Sorted BAM output**: Ready for downstream analysis (variant calling, QC, etc.)
- **Modular design**: Integrates with other WILDS workflows and tools
- **Scalable**: Supports batch alignment across many samples
- **Flexible**: Customizable resource settings per task
- **Robust**: Includes error handling and reproducible output filenames

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Subsets of real human genome data (e.g., chr22)
- Automated output verification

For development, contributions, and collaboration, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, contact the Fred Hutch Data Science Lab (DaSL) at [wilds@fredhutch.org](mailto:wilds@fredhutch.org) or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

To contribute to this module, please review the [WILDS Contributor Guide](https://getwilds.org/guide/) and the [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md).

## License

Distributed under the MIT License. See `LICENSE` for full details.