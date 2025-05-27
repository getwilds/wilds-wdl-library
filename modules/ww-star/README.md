# ww-star
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for RNA-seq alignment using STAR's two-pass methodology.

## Overview

This module provides reusable WDL tasks for high-quality RNA-seq alignment using STAR (Spliced Transcripts Alignment to a Reference). It implements STAR's two-pass approach for improved splice junction detection and includes comprehensive validation of outputs.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger RNA-seq analysis pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `build_index`, `align_two_pass`, `validate_outputs`
- **Workflow**: `star_example` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/star:2.7.6a`

## Tasks

### `build_index`
Builds STAR genome index from reference FASTA and GTF files.

**Inputs:**
- `reference_fasta` (File): Reference genome FASTA file
- `reference_gtf` (File): Reference genome GTF annotation file
- `sjdb_overhang` (Int): Splice junction database overhang (default: 100)
- `genome_sa_index_nbases` (Int): SA index string length (default: 14)
- `memory_gb` (Int): Memory allocation in GB (default: 64)
- `cpu_cores` (Int): Number of CPU cores (default: 8)

**Outputs:**
- `star_index_tar` (File): Compressed STAR genome index

### `align_two_pass`
Performs two-pass STAR alignment with gene counting.

**Inputs:**
- `star_genome_tar` (File): STAR genome index from `build_index`
- `sample_data` (SampleInfo): Sample information struct
- `ref_genome_name` (String): Reference genome name for output files
- Various resource and parameter settings

**Outputs:**
- `bam` (File): Sorted BAM alignment file
- `bai` (File): BAM index file
- `gene_counts` (File): Gene-level read counts
- `log_final`, `log_progress`, `log` (Files): STAR log files
- `sj_out` (File): Splice junction file

### `validate_outputs`
Validates alignment outputs and generates a comprehensive report.

**Inputs:**
- Arrays of output files from alignment tasks

**Outputs:**
- `report` (File): Validation summary with alignment statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "path/to/ww-star.wdl" as star_tasks

struct SampleInfo {
    String name
    File r1
    File r2
}

workflow my_rna_seq_pipeline {
  input {
    Array[SampleInfo] samples
    File reference_fasta
    File reference_gtf
  }
  
  call star_tasks.build_index {
    input:
      reference_fasta = reference_fasta,
      reference_gtf = reference_gtf
  }
  
  scatter (sample in samples) {
    call star_tasks.align_two_pass {
      input:
        sample_data = sample,
        star_genome_tar = build_index.star_index_tar,
        ref_genome_name = "hg38"
    }
  }
  
  output {
    Array[File] aligned_bams = align_two_pass.bam
    Array[File] gene_counts = align_two_pass.gene_counts
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-sra**: Download sequencing data before alignment
- **ww-sra-star vignette**: Complete SRA-to-alignment pipeline
- **ww-star-deseq2 vignette**: Extended pipeline with differential expression

## Testing the Module

The module includes a demonstration workflow with comprehensive testing:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-star.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-star.wdl -i inputs.json

# Using Sprocket
sprocket run ww-star.wdl inputs.json
```

### Test Input Format

```json
{
  "star_example.samples": [
    {
      "name": "sample1",
      "r1": "/path/to/sample1_R1.fastq.gz",
      "r2": "/path/to/sample1_R2.fastq.gz"
    }
  ],
  "star_example.reference_genome": {
    "name": "hg38",
    "fasta": "/path/to/genome.fasta",
    "gtf": "/path/to/annotation.gtf"
  },
  "star_example.genome_sa_index_nbases": 14,
  "star_example.cpus": 8,
  "star_example.memory_gb": 64
}
```

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: Scales with genome size (64GB recommended for human genome)
- **CPUs**: Adjustable based on available resources (8 cores recommended)
- **Index Parameters**: `genome_sa_index_nbases` should be set based on genome size:
  - Human genome (~3GB): 14 (default)
  - Mouse genome (~2.7GB): 14
  - C. elegans (~100MB): 11
  - Small test genomes: 8-11

### Advanced Parameters

- `sjdb_overhang`: Set to (read_length - 1) for optimal junction detection
- Resource allocation can be tuned for different compute environments

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Sufficient memory for genome indexing (varies by genome size)

## Features

- **Two-pass alignment**: Improved splice junction detection
- **Comprehensive outputs**: BAM files, gene counts, QC metrics
- **Validation**: Built-in output validation and reporting
- **Scalable**: Configurable resource allocation
- **Robust**: Extensive error handling and cleanup
- **Compatible**: Works with multiple WDL executors

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real RNA-seq data (chromosome 22 subset for efficiency)
- Comprehensive validation of all outputs

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
