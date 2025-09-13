# ww-star
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for RNA-seq alignment using STAR's two-pass methodology.

## Overview

This module provides reusable WDL tasks for high-quality RNA-seq alignment using STAR (Spliced Transcripts Alignment to a Reference). It implements STAR's two-pass methodology for optimal splice junction detection and includes comprehensive validation of outputs. The module supports both single-sample analysis and batch processing.

The module can run completely standalone with automatic test data download, or integrate with existing FASTQ files for production analyses.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `build_index`, `align_two_pass`, `validate_outputs`
- **Workflow**: `star_example` (demonstration workflow with automatic test data support)
- **Container**: `getwilds/star:2.7.6a`
- **Dependencies**: Integrates with `ww-testdata` module for complete workflows
- **Test Data**: Automatically downloads reference genome, GTF annotation, and FASTQ data when not provided

## Tasks

### `build_index`
Builds STAR genome index from reference FASTA and GTF files.

**Inputs:**
- `reference_fasta` (File): Reference genome FASTA file
- `reference_gtf` (File): Reference genome GTF annotation file
- `sjdb_overhang` (Int): Length of genomic sequence around junctions (default: 100)
- `genome_sa_index_nbases` (Int): SA pre-indexing string length (default: 14)
- `memory_gb` (Int): Memory allocation in GB (default: 64)
- `cpu_cores` (Int): Number of CPU cores (default: 8)

**Outputs:**
- `star_index_tar` (File): Compressed tarball containing STAR genome index

### `align_two_pass`
Performs RNA-seq alignment using STAR's two-pass methodology.

**Inputs:**
- `star_genome_tar` (File): STAR genome index from `build_index`
- `r1` (File): R1 FASTQ file
- `r2` (File): R2 FASTQ file
- `name` (String): Sample name for output files
- `sjdb_overhang` (Int): Length of genomic sequence around junctions (default: 100)
- `memory_gb` (Int): Memory allocation in GB (default: 62)
- `cpu_cores` (Int): Total CPU cores (default: 8)
- `star_threads` (Int): STAR-specific thread count (default: 6)

**Outputs:**
- `bam` (File): Sorted BAM alignment file
- `bai` (File): BAM index file
- `gene_counts` (File): Gene-level read counts
- `log_final`, `log_progress`, `log` (Files): STAR log files
- `sj_out` (File): Splice junction file

### `validate_outputs`
Validates alignment outputs and generates a comprehensive report.

**Inputs:**
- `bam_files` (Array[File]): Array of BAM files to validate
- `bai_files` (Array[File]): Array of BAM index files to validate
- `gene_count_files` (Array[File]): Array of gene count files to validate

**Outputs:**
- `report` (File): Validation summary with alignment statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-star/ww-star.wdl" as star_tasks

struct StarSample {
    String name
    File r1
    File r2
}

workflow my_rna_seq_pipeline {
  input {
    Array[StarSample] samples
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
        star_genome_tar = build_index.star_index_tar,
        r1 = sample.r1,
        r2 = sample.r2,
        name = sample.name
    }
  }
  
  output {
    Array[File] aligned_bams = align_two_pass.bam
    Array[File] gene_counts = align_two_pass.gene_counts
  }
}
```

### Advanced Usage Examples

**Custom splice junction parameters:**
```wdl
call star_tasks.build_index {
  input:
    sjdb_overhang = 149,  # read_length - 1 for optimal performance
    genome_sa_index_nbases = 14  # Adjust for genome size
}
```

**Resource optimization:**
```wdl
call star_tasks.align_two_pass {
  input:
    memory_gb = 64,
    cpu_cores = 16,
    star_threads = 12  # Leave some cores for I/O
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-sra**: Download sequencing data before alignment
- **ww-testdata**: Automatic provisioning of reference data and test samples
- **ww-sra-star vignette**: Complete SRA-to-alignment pipeline
- **Custom workflows**: Foundation for RNA-seq analysis pipelines

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

The demonstration workflow automatically downloads test data and runs without requiring input files:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-star.wdl

# Using miniWDL
miniwdl run ww-star.wdl

# Using Sprocket
sprocket run ww-star.wdl
```

The demonstration workflow (`star_example`) automatically:
1. Downloads reference genome data using `ww-testdata`
2. Downloads demonstration FASTQ data using `ww-testdata`
3. Builds STAR genome index
4. Performs RNA-seq alignment using STAR two-pass methodology
5. Validates all outputs

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 8-64 GB recommended (scales with genome size and sample complexity)
- **CPUs**: 8 cores recommended; STAR benefits from multi-threading
- **Index Parameters**: `genome_sa_index_nbases` should be set based on genome size:
  - Human genome (~3GB): 14 (default)
  - Mouse genome (~2.7GB): 14
  - C. elegans (~100MB): 11
  - Small test genomes: 8-11

### Advanced Parameters

- `sjdb_overhang`: Set to (read_length - 1) for optimal junction detection
- `star_threads`: Usually set slightly less than `cpu_cores` to leave resources for I/O
- Resource allocation can be tuned for different compute environments


## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Docker/Apptainer support for containerized execution
- Sufficient computational resources (STAR can be memory-intensive for large genomes)

## Features

- **Standalone execution**: Complete workflow with automatic test data download
- **Automatic test data**: Uses test data from `ww-testdata` module for demonstration
- **Two-pass methodology**: Optimal splice junction detection using STAR's two-pass approach
- **Comprehensive outputs**: BAM files, gene counts, splice junctions, and detailed logs
- **Multi-sample support**: Process multiple samples in parallel
- **Validation**: Built-in output validation and reporting
- **Module integration**: Seamlessly combines with ww-sra and ww-testdata
- **Scalable**: Configurable resource allocation
- **Compatible**: Works with multiple WDL executors

## Performance Considerations

- **Memory usage**: Index building requires significant RAM (32-64GB for human genome)
- **CPU scaling**: Both index building and alignment benefit from multiple cores
- **Storage requirements**: Ensure sufficient space for index files and BAM outputs
- **Two-pass optimization**: Second pass uses splice junctions from first pass for improved accuracy

## Output Description

- **BAM files**: Sorted alignment files ready for downstream analysis
- **BAI files**: Index files for rapid BAM access
- **Gene count files**: Tab-delimited files with read counts per gene
- **Log files**: Detailed alignment statistics and runtime information
- **Splice junction files**: Coordinates and support for detected splice junctions
- **Validation report**: Comprehensive validation with alignment statistics and file integrity checks

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (demonstration FASTQ for integration testing)
- Comprehensive validation of all outputs including BAM format validation
- Integration testing with ww-testdata modules

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to STAR usage or configuration, please refer to the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). Please make sure to cite their work if you use STAR in your analyses:

Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
