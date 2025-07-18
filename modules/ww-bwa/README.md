# ww-bwa
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for sequence alignment using the Burrows-Wheeler Aligner (BWA).

## Overview

This module provides reusable WDL tasks for aligning sequencing reads to a reference genome using **BWA-MEM** (Burrows-Wheeler aligner maximal exact matches algorithm). It includes comprehensive indexing, alignment, and validation capabilities.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger sequence analysis pipelines. It can run completely standalone with automatic test data download, or integrate with existing FASTQ files.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `bwa_index`, `bwa_mem`, `validate_outputs`
- **Workflow**: `bwa_example` (demonstration workflow with automatic test data support)
- **Container**: `getwilds/bwa:0.7.17`
- **Test Data**: Automatically downloads reference genome and FASTQ data when not provided using `ww-testdata` module

## Tasks

### `bwa_index`
Builds BWA index files from reference FASTA and packages them in a compressed tarball.

**Inputs:**
- `reference_fasta` (File): Reference genome FASTA file
- `cpu_cores` (Int): Number of CPU cores (default: 8)
- `memory_gb` (Int): Memory allocation in GB (default: 32)

**Outputs:**
- `bwa_index_tar` (File): Compressed tarball containing BWA genome index files

### `bwa_mem`
Aligns paired-end reads to a reference using BWA-MEM with automatic read group addition and BAM sorting.

**Inputs:**
- `bwa_genome_tar` (File): Compressed tarball containing BWA genome index
- `reference_fasta` (File): Reference genome FASTA file
- `reads` (File): FASTQ file for forward (R1) reads or interleaved reads
- `mates` (File): Optional FASTQ file for reverse (R2) reads
- `name` (String): Sample name for read group information and output files
- `paired_end` (Boolean): Optional, indicating if reads are paired end (default: true)
- `cpu_cores` (Int): Number of CPU cores (default: 8)
- `memory_gb` (Int): Memory allocation in GB (default: 16)

**Outputs:**
- `sorted_bam` (File): Sorted BAM alignment file
- `sorted_bai` (File): BAM index file

### `validate_outputs`
Validates alignment outputs and generates a comprehensive report.

**Inputs:**
- `bam_files` (Array[File]): Array of BAM files to validate
- `bai_files` (Array[File]): Array of BAM index files to validate

**Outputs:**
- `report` (File): Validation summary with alignment statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bwa/ww-bwa.wdl" as bwa_tasks

struct BwaSample {
    String name
    File reads
    File mates
}

workflow my_alignment_pipeline {
  input {
    Array[BwaSample] samples
    File reference_fasta
  }

  call bwa_tasks.bwa_index {
    input:
      reference_fasta = reference_fasta,
      cpu_cores = 8,
      memory_gb = 32
  }

  scatter (sample in samples) {
    call bwa_tasks.bwa_mem {
      input:
        bwa_genome_tar = bwa_index.bwa_index_tar,
        reference_fasta = reference_fasta,
        reads = sample.reads,
        mates = sample.mates,
        name = sample.name,
        cpu_cores = 8,
        memory_gb = 16
    }
  }

  call bwa_tasks.validate_outputs {
    input:
      bam_files = bwa_mem.sorted_bam,
      bai_files = bwa_mem.sorted_bai
  }

  output {
    Array[File] aligned_bams = bwa_mem.sorted_bam
    Array[File] aligned_bais = bwa_mem.sorted_bai
    File validation_report = validate_outputs.report
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Provides reference genomes and test FASTQ files automatically
- **ww-bcftools**: Use aligned BAM files for variant calling
- **ww-bedtools**: Use aligned BAM files for genomic interval analysis
- **ww-delly**: Use aligned BAM files for structural variant calling
- **Custom workflows**: Foundation for any analysis requiring aligned reads

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

### Automatic Demo Mode

When no samples or reference files are provided, the workflow automatically:
1. Downloads reference genome data using `ww-testdata`
2. Downloads FASTQ test data using `ww-testdata`
3. Builds BWA index from reference FASTA
4. Aligns reads using BWA-MEM with proper read groups
5. Sorts and indexes BAM output
6. Validates all outputs with detailed statistics

### Test Input Format

**Minimal input (uses automatic demo data):**
```json
{
  "bwa_example.demo_sra_id": "ERR1258306",
  "bwa_example.cpus": 2,
  "bwa_example.memory_gb": 8
}
```

**Full input (provide your own data):**
```json
{
  "bwa_example.samples": [
    {
      "name": "sample1",
      "reads": "/path/to/sample1_R1.fastq.gz",
      "mates": "/path/to/sample1_R2.fastq.gz"
    }
  ],
  "bwa_example.reference_fasta": "/path/to/genome.fasta",
  "bwa_example.cpus": 8,
  "bwa_example.memory_gb": 32
}
```

**Note**: If no `samples` or `reference_fasta` are provided, the workflow will automatically download test data using the `ww-testdata` module.

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:
- **Memory**: 16-32 GB recommended for human genomes; can be reduced for smaller references
- **CPUs**: 8 cores optimal for most samples; BWA-MEM and samtools benefit from multi-threading
- **Storage**: Sufficient space for FASTQ files, reference genome, BWA index, and output BAM files

### Advanced Parameters

- **BWA index building**: Memory-intensive step; ensure adequate memory allocation (32GB recommended for human genome)
- **CPU thread optimization**: BWA-MEM uses `cpu_cores - 1` threads automatically for optimal performance
- **Read group information**: Automatically added using the sample name for downstream compatibility
- **SAM to BAM conversion**: Includes sorting and indexing in a single step

### Demo Configuration

- `demo_sra_id`: Currently a placeholder parameter; actual demo uses test FASTQ data from `ww-testdata`
- Resource parameters apply to both demo and user-provided data modes

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker or Apptainer support
- Sufficient memory for genome indexing (varies by genome size; 32GB for human genome)
- Paired-end FASTQ files (when providing your own data)

## Features

- **Standalone execution**: Complete workflow with automatic test data download
- **Optional input**: Use your own data or automatic demo data from `ww-testdata`
- **Paired-end or single-end input**: Use separate or interleaved paired-end FASTQs, or single-end data.
- **BWA-MEM alignment**: Fast and accurate for reads 70bp to 1Mbp
- **Automatic indexing**: Builds BWA index from reference FASTA with tarball packaging
- **Sorted BAM output**: Ready for downstream analysis (variant calling, QC, etc.)
- **Read group addition**: Proper read group tags for downstream compatibility
- **Comprehensive validation**: Built-in output validation with alignment statistics
- **Module integration**: Seamlessly integrates with `ww-testdata` and other WILDS modules
- **Scalable**: Supports batch alignment across many samples
- **Robust**: Includes error handling and thread optimization

## Advanced Usage

### Resource Optimization for Large Genomes

For human genome alignment:

```json
{
  "bwa_example.cpus": 16,
  "bwa_example.memory_gb": 64
}
```

### Multiple Sample Processing

```json
{
  "bwa_example.samples": [
    {
      "name": "control_1",
      "reads": "/data/control_1_R1.fastq.gz",
      "mates": "/data/control_1_R2.fastq.gz"
    },
    {
      "name": "treatment_1",
      "reads": "/data/treatment_1_interleaved.fastq.gz",
      "mates": ""
    }
  ],
  "bwa_example.cpus": 12,
  "bwa_example.memory_gb": 48
}
```

### Custom Reference Genome

```json
{
  "bwa_example.reference_fasta": "/data/custom_genome.fasta",
  "bwa_example.samples": [
    {
      "name": "sample1",
      "reads": "/data/sample1_R1.fastq.gz",
      "mates": "/data/sample1_R2.fastq.gz"
    }
  ]
}
```

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Test data from `ww-testdata` module (chromosome 1 subset for efficiency)
- Comprehensive validation of all outputs including alignment statistics
- Integration testing with `ww-testdata` module
- Cross-platform compatibility testing

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, contact the Fred Hutch Data Science Lab (DaSL) at [wilds@fredhutch.org](mailto:wilds@fredhutch.org) or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

To contribute to this module, please review the [WILDS Contributor Guide](https://getwilds.org/guide/) and the [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md).

## License

Distributed under the MIT License. See `LICENSE` for full details.
