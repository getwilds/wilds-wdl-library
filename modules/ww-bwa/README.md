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
- **Dependencies**: Integrates with `ww-sra` and `ww-testdata` modules for complete workflows
- **Test Data**: Automatically downloads reference genome and SRA data when not provided

## Tasks

### `bwa_index`
Builds BWA index files from reference FASTA.

**Inputs:**
- `reference_fasta` (File): Reference genome FASTA file
- `cpu_cores` (Int): Number of CPU cores (default: 8)
- `memory_gb` (Int): Memory allocation in GB (default: 32)

**Outputs:**
- `bwa_index_tar` (File): Compressed tarball containing BWA genome index

### `bwa_mem`

Aligns paired-end reads to a reference using BWA-MEM.

**Inputs:**
- `bwa_genome_tar` (File): Compressed tarball containing BWA genome index
- `reference_fasta` (File): Reference genome FASTA file
- `r1` (File): FASTQ file for read 1
- `r2` (File): FASTQ file for read 2
- `name` (String): Sample name for read group information
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
    File r1
    File r2
}

workflow my_alignment_pipeline {
  input {
    Array[BwaSample] samples
    File reference_fasta
  }

  call bwa_tasks.bwa_index {
    input: 
      reference_fasta = reference_fasta
  }

  scatter (sample in samples) {
    call bwa_tasks.bwa_mem {
      input:
        bwa_genome_tar = bwa_index.bwa_index_tar,
        reference_fasta = reference_fasta,
        r1 = sample.r1,
        r2 = sample.r2,
        name = sample.name
    }
  }

  output {
    Array[File] aligned_bams = bwa_mem.sorted_bam
    Array[File] aligned_bais = bwa_mem.sorted_bai
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-sra**: Download sequencing data prior to alignment (built into demo workflow)
- **ww-bcftools**: Use aligned BAM files for variant calling
- **ww-bedtools**: Use aligned BAM files for genomic interval analysis
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
2. Downloads SRA data (default: ERR1258306) using `ww-sra`
3. Builds BWA index
4. Aligns reads using BWA-MEM
5. Validates outputs

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
      "r1": "/path/to/sample1_R1.fastq.gz",
      "r2": "/path/to/sample1_R2.fastq.gz"
    }
  ],
  "bwa_example.reference_fasta": "/path/to/genome.fasta",
  "bwa_example.cpus": 8,
  "bwa_example.memory_gb": 32
}
```

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:
- **Memory**: 16-32 GB recommended for human genomes; can be tuned for smaller references
- **CPUs**: 8 cores typically sufficient for most samples; BWA-MEM and samtools scale with available threads
- **Storage**: Sufficient space for FASTQ files, reference genome, and output BAM files

### Advanced Parameters

- BWA index building is memory-intensive; ensure adequate memory allocation
- CPU threads are automatically calculated as cpu_cores - 1 for optimal performance
- Read group information is automatically added using the sample name
- Set `cpu_cores` and `memory_gb` based on your environment and dataset size

### Demo Configuration

- `demo_sra_id`: Change to use different SRA sample for testing
- Resource parameters apply to both demo and user-provided data modes

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker or Apptainer support
- Sufficient memory for genome indexing (varies by genome size)
- Paired-end FASTQ files (when providing your own data)

## Features

- **Standalone execution**: Complete workflow with automatic test data download
- **Flexible input**: Use your own data or automatic demo data
- **BWA-MEM alignment**: Fast and accurate for reads 70bp to 1Mbp
- **Automatic indexing**: Builds BWA index from reference FASTA
- **Sorted BAM output**: Ready for downstream analysis (variant calling, QC, etc.)
- **Read group addition**: Proper read group tags for downstream compatibility
- **Validation**: Built-in output validation and reporting with alignment statistics
- **Module integration**: Seamlessly combines with ww-sra and ww-testdata
- **Modular design**: Integrates with other WILDS workflows and tools
- **Scalable**: Supports batch alignment across many samples
- **Flexible**: Customizable resource settings per task
- **Robust**: Includes error handling and reproducible output filenames

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
      "r1": "/data/control_1_R1.fastq.gz",
      "r2": "/data/control_1_R2.fastq.gz"
    },
    {
      "name": "treatment_1", 
      "r1": "/data/treatment_1_R1.fastq.gz",
      "r2": "/data/treatment_1_R2.fastq.gz"
    }
  ]
}
```

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (SRA sample ERR1258306 for integration testing)
- Comprehensive validation of all outputs including alignment statistics
- Integration testing with ww-sra and ww-testdata modules
- Chromosome 22 subset for efficiency during CI/CD

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, contact the Fred Hutch Data Science Lab (DaSL) at [wilds@fredhutch.org](mailto:wilds@fredhutch.org) or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

To contribute to this module, please review the [WILDS Contributor Guide](https://getwilds.org/guide/) and the [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md).

## License

Distributed under the MIT License. See `LICENSE` for full details.
