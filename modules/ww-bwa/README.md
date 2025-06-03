# ww-bwa
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)  
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for sequence alignment using the Burrows-Wheeler Aligner (BWA).

## Overview

This module provides reusable WDL tasks for aligning sequencing reads to a reference genome using **BWA-MEM** (Burrows-Wheeler aligner maximal exact matches algorithm). 

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger sequence analysis pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Task**: `bwa_index`, `bwa_mem`
- **Workflow**: `bwa_example` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/bwa:0.7.17`

## Tasks

### `bwa_index`
Builds BWA index files from reference FASTA.

**Inputs:**
- `reference_fasta` (File): Reference genome FASTA file
- `memory_gb` (Int): Memory allocation in GB (default: 32)
- `cpu_cores` (Int): Number of CPU cores (default: 8)

**Outputs:**
- `fasta` (File): "Reference genome FASTA file"
- `amb` (File): Text file of ambiguous bases
- `ann` (File): Text file of reference sequence information, such as name and length
- `bwt` (File): Binary file of Burrows-Wheeler transformed reference sequence
- `pac` (File): Binary file of compressed reference sequence
- `sa` (File): Binary file of the suffix array

### `bwa_mem`

Aligns paired-end reads to a reference using BWA-MEM.

**Inputs:**
- `sample_data` (SampleInfo): Sample information struct containing sample name and R1/R2 FASTQ file paths
- `reference_fasta` (File): Reference genome FASTA file (indexed and preprocessed externally)
- `cpu_cores` (Int): Number of CPU cores (default: 8)
- `memory_gb` (Int): Memory allocation in GB (default: 16)

**Outputs:**
- `sorted_bam` (File): Sorted BAM alignment file
- `sorted_bai` (File): BAM index file

## Usage as a Module

### Importing into Your Workflow

```wdl
import "path/to/ww-bwa.wdl" as bwa_tasks

struct SampleInfo {
    String name
    File r1
    File r2
}

workflow my_alignment_pipeline {
  input {
    Array[SampleInfo] samples
    File reference_fasta
  }

  call bwa_tasks.bwa_index {
    reference_fasta = reference_fasta
  }

  scatter (sample in samples) {
    call bwa_tasks.bwa_mem {
      input:
        sample_data = sample,
        reference_fasta = bwa_index.fasta
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
    "ref_fasta": "/path/to/genome.fasta",
  },
  "bwa_example.cpus": 8,
  "bwa_example.memory_gb": 32
}
```

## Configuration Guidelines

### Resource Allocation

- **Memory**: 16-32 GB recommended for human genomes; can be tuned for smaller references
- **CPUs**: 8 cores typically sufficient for most samples; BWA-MEM and samtools scale with available threads

### Advanced Parameters

- Ensure the reference genome is pre-indexed with BWA and has an associated `.fai` and `.dict`
- Set `cpu_cores` and `memory_gb` based on your environment and dataset size

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker or Apptainer support
- Sufficient memory for genome indexing (varies by genome size)

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
- Real RNA-seq data (chromosome 22 subset for efficiency)

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, contact the Fred Hutch Data Science Lab (DaSL) at [wilds@fredhutch.org](mailto:wilds@fredhutch.org) or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

To contribute to this module, please review the [WILDS Contributor Guide](https://getwilds.org/guide/) and the [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md).

## License

Distributed under the MIT License. See `LICENSE` for full details.