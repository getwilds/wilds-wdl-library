# ww-bowtie Module

[![Project Status: Prototype -- Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for performing short-read alignment using [Bowtie](https://bowtie-bio.sourceforge.net/index.shtml), an ultrafast, memory-efficient short-read aligner.

## Overview

This module provides WDL tasks for building Bowtie genome indexes and aligning short sequencing reads to a reference genome. Bowtie is particularly well-suited for aligning short reads (up to ~50 bp) with high speed and low memory usage, making it a standard tool for ChIP-seq, small RNA-seq, and other short-read applications.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-bowtie.wdl` - Contains task definitions for Bowtie indexing and alignment
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `bowtie_build`

Builds Bowtie index files from a reference FASTA file.

**Inputs:**
- `reference_fasta` (File): Reference genome FASTA file
- `index_prefix` (String, default="bowtie_index"): Prefix for the Bowtie index files
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=16): Memory allocated for the task in GB

**Outputs:**
- `bowtie_index_tar` (File): Compressed tarball containing Bowtie genome index files

### `bowtie_align`

Aligns short reads to a reference genome using Bowtie.

**Inputs:**
- `bowtie_index_tar` (File): Compressed tarball containing Bowtie genome index files
- `index_prefix` (String, default="bowtie_index"): Prefix used when building the Bowtie index
- `reads` (File): FASTQ file for forward (R1) reads
- `name` (String): Sample name for output file naming and read group information
- `mates` (File, optional): FASTQ file for reverse (R2) reads for paired-end alignment
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `sorted_bam` (File): Sorted Bowtie alignment output BAM file
- `sorted_bai` (File): Index file for the sorted Bowtie alignment BAM file

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bowtie/ww-bowtie.wdl" as bowtie_tasks

struct BowtieSample {
    String name
    File reads
    File? mates
}

workflow my_alignment_pipeline {
  input {
    File reference_fasta
    Array[BowtieSample] samples
  }

  call bowtie_tasks.bowtie_build { input:
      reference_fasta = reference_fasta
  }

  scatter (sample in samples) {
    call bowtie_tasks.bowtie_align { input:
        bowtie_index_tar = bowtie_build.bowtie_index_tar,
        reads = sample.reads,
        mates = sample.mates,
        name = sample.name
    }
  }

  output {
    Array[File] aligned_bams = bowtie_align.sorted_bam
    Array[File] aligned_bais = bowtie_align.sorted_bai
  }
}
```

### Advanced Usage Examples

**Custom resource allocation:**
```wdl
call bowtie_tasks.bowtie_align {
  input:
    bowtie_index_tar = bowtie_build.bowtie_index_tar,
    reads = sample_reads,
    name = "large_sample",
    cpu_cores = 8,
    memory_gb = 16
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of test data for demonstrations
- **ww-samtools**: Post-alignment BAM processing and statistics
- **ww-fastqc**: Pre-alignment quality control of FASTQ files

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint bowtie_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads reference genome data using `ww-testdata`
2. Downloads test FASTQ data using `ww-testdata`
3. Builds a Bowtie index from the reference
4. Aligns both paired-end and single-end test samples
5. Produces sorted BAM files and their indexes

## Docker Container

This module uses the `getwilds/bowtie:1.3.1` container image, which includes:
- Bowtie aligner (v1.3.1)
- Samtools for BAM conversion and sorting
- All necessary system dependencies

## Citation

> Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology 10:R25 (2009).
> DOI: [10.1186/gb-2009-10-3-r25](https://doi.org/10.1186/gb-2009-10-3-r25)

## Parameters and Resource Requirements

### Default Resources

| Task | CPU | Memory | Typical Runtime |
|------|-----|--------|-----------------|
| `bowtie_build` | 4 cores | 16 GB | 5-30 min (varies by genome size) |
| `bowtie_align` | 4 cores | 8 GB | 5-30 min per sample |

### Resource Scaling
- `cpu_cores`: Bowtie supports multithreaded alignment; increase for faster processing
- `memory_gb`: Index building requires more memory for larger genomes; alignment is generally memory-efficient

## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for detailed guidelines

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
- **[Bowtie Manual](https://bowtie-bio.sourceforge.net/manual.shtml)**: Official Bowtie documentation
