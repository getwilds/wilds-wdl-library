# ww-bowtie2 Module

[![Project Status: Prototype -- Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for performing sequence alignment using [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), a fast and sensitive read aligner.

## Overview

This module provides WDL tasks for building Bowtie 2 genome indexes and aligning sequencing reads to a reference genome. Bowtie 2 supports gapped, local, and paired-end alignment modes, making it suitable for a wide range of applications including DNA-seq, RNA-seq (when used with a splice-unaware approach), ChIP-seq, and other short-read sequencing experiments. It handles reads of varying lengths and is the successor to Bowtie.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-bowtie2.wdl` - Contains task definitions for Bowtie 2 indexing and alignment
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `bowtie2_build`

Builds Bowtie 2 index files from a reference FASTA file.

**Inputs:**
- `reference_fasta` (File): Reference genome FASTA file
- `index_prefix` (String, default="bowtie2_index"): Prefix for the Bowtie 2 index files
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=16): Memory allocated for the task in GB

**Outputs:**
- `bowtie2_index_tar` (File): Compressed tarball containing Bowtie 2 genome index files

### `bowtie2_align`

Aligns reads to a reference using Bowtie 2, optionally filtering the output using Samtools

**Inputs:**
- `bowtie2_index_tar` (File): Compressed tarball containing Bowtie 2 genome index files
- `index_prefix` (String, default="bowtie2_index"): Prefix used when building the Bowtie 2 index
- `reads` (File): FASTQ file for forward (R1) reads
- `name` (String): Sample name for output file naming and read group information
- `mates` (File, optional): FASTQ file for reverse (R2) reads for paired-end alignment
- `preset` (String, optional): Bowtie 2 sensitivity preset. One of `fast`, `sensitive`, `very-sensitive`, `fast-local`, `sensitive-local`, `very-sensitive-local`. Unset means bowtie2's default end-to-end `--sensitive` preset.
- `capture_unaligned` (Boolean, default=false): If true, write reads that fail to align concordantly to gzipped FASTQ outputs (`--un-gz` / `--un-conc-gz`). Used for e.g. rRNA depletion upstream of a second alignment step.
- `min_mapq` (Int, default=0): Minimum MAPQ score for `samtools view -q` post-alignment filter. 0 disables.
- `samtools_filter_flags` (String, default=""): Extra flags passed to `samtools view` for post-alignment filtering (e.g. `-f 2` to keep only proper pairs).
- `extra_bowtie2_args` (String, default=""): Additional arguments forwarded verbatim to bowtie2 (e.g. `--no-mixed --no-discordant`).
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `sorted_bam` (File): Sorted Bowtie 2 alignment output BAM file
- `sorted_bai` (File): Index file for the sorted Bowtie 2 alignment BAM file
- `unaligned_se` (File?, optional): FASTQ of unaligned reads. Only produced when `capture_unaligned` is true and `mates` is unset (single-end input).
- `unaligned_r1` (File?, optional): FASTQ of R1 reads that failed to align concordantly. Only produced when `capture_unaligned` is true and `mates` is set (paired-end input).
- `unaligned_r2` (File?, optional): FASTQ of R2 reads that failed to align concordantly. Only produced when `capture_unaligned` is true and `mates` is set (paired-end input).

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bowtie2/ww-bowtie2.wdl" as bowtie2_tasks

struct Bowtie2Sample {
    String name
    File reads
    File? mates
}

workflow my_alignment_pipeline {
  input {
    File reference_fasta
    Array[Bowtie2Sample] samples
  }

  call bowtie2_tasks.bowtie2_build { input:
      reference_fasta = reference_fasta
  }

  scatter (sample in samples) {
    call bowtie2_tasks.bowtie2_align { input:
        bowtie2_index_tar = bowtie2_build.bowtie2_index_tar,
        reads = sample.reads,
        mates = sample.mates,
        name = sample.name
    }
  }

  output {
    Array[File] aligned_bams = bowtie2_align.sorted_bam
    Array[File] aligned_bais = bowtie2_align.sorted_bai
  }
}
```

### Advanced Usage Examples

**Custom resource allocation:**
```wdl
call bowtie2_tasks.bowtie2_align {
  input:
    bowtie2_index_tar = bowtie2_build.bowtie2_index_tar,
    reads = sample_reads,
    name = "large_sample",
    cpu_cores = 8,
    memory_gb = 16
}
```

**rRNA depletion (capture unaligned reads for a second alignment step):**
```wdl
call bowtie2_tasks.bowtie2_align as deplete_rrna {
  input:
    bowtie2_index_tar = rrna_index_tar,
    reads = r1, mates = r2,
    name = "sample",
    preset = "fast-local",
    capture_unaligned = true
}
# Pass deplete_rrna.unaligned_r1 / unaligned_r2 into the next bowtie2_align call
# against the experimental genome.
```

**PRO-seq spike-in alignment (sensitive local + proper-pair + MAPQ filter):**
```wdl
call bowtie2_tasks.bowtie2_align as align_spikein {
  input:
    bowtie2_index_tar = spikein_index_tar,
    reads = trimmed_r1, mates = trimmed_r2,
    name = "sample_spikein",
    preset = "very-sensitive-local",
    min_mapq = 10,
    samtools_filter_flags = "-f 2",
    extra_bowtie2_args = "--no-mixed --no-discordant"
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of test data for demonstrations
- **ww-samtools**: Post-alignment BAM processing and statistics
- **ww-fastqc**: Pre-alignment quality control of FASTQ files
- **ww-fastp**: Adapter trimming and (optionally) UMI extraction; trimmed reads feed directly into `bowtie2_align`
- **ww-umi-tools**: Downstream PCR-duplicate removal when fastp encoded UMIs into read names

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint bowtie2_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads reference genome data using `ww-testdata`
2. Downloads test FASTQ data using `ww-testdata`
3. Builds a Bowtie 2 index from the reference
4. Aligns both paired-end and single-end test samples
5. Produces sorted BAM files and their indexes

## Docker Container

This module uses the `getwilds/bowtie2:2.5.4` container image, which includes:
- Bowtie 2 aligner (v2.5.4)
- Samtools for BAM conversion and sorting
- All necessary system dependencies

## Citation

> Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012;9(4):357-359.
> DOI: [10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)

## Parameters and Resource Requirements

### Default Resources

| Task | CPU | Memory | Typical Runtime |
|------|-----|--------|-----------------|
| `bowtie2_build` | 4 cores | 16 GB | 5-30 min (varies by genome size) |
| `bowtie2_align` | 4 cores | 8 GB | 5-30 min per sample |

### Resource Scaling
- `cpu_cores`: Bowtie 2 supports multithreaded alignment; increase for faster processing
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
- **[Bowtie 2 Manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)**: Official Bowtie 2 documentation
