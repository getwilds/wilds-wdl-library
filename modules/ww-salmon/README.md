# ww-salmon
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for RNA-seq transcript quantification using Salmon.

## Overview

This module provides reusable WDL tasks for fast and accurate transcript-level quantification using Salmon. Salmon uses a quasi-mapping approach for rapid quantification of transcript abundance from RNA-seq data without requiring alignment to a reference genome. The module supports both single-sample and batch processing, and includes utilities for merging quantification results across multiple samples.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `build_index`, `quantify`, `merge_results`
- **Test workflow**: `testrun.wdl` (demonstration workflow with automatic test data support)
- **Container**: `getwilds/salmon:1.10.3` (WILDS Docker image with Salmon installed)

## Tasks

### `build_index`
Builds Salmon index from a reference transcriptome FASTA file.

**Inputs:**
- `transcriptome_fasta` (File): Reference transcriptome in FASTA format
- `cpu_cores` (Int): Number of CPU cores (default: 8)
- `memory_gb` (Int): Memory allocation in GB (default: 16)

**Outputs:**
- `salmon_index` (File): Compressed tarball containing the Salmon index for quantification

### `quantify`
Quantifies transcript expression from paired-end RNA-seq reads using Salmon.

**Inputs:**
- `salmon_index_dir` (File): Compressed tarball containing Salmon index from `build_index`
- `sample_name` (String): Sample name identifier for output files
- `fastq_r1` (File): FASTQ file for read 1
- `fastq_r2` (File): FASTQ file for read 2
- `cpu_cores` (Int): Number of CPU cores (default: 8)
- `memory_gb` (Int): Memory allocation in GB (default: 16)

**Outputs:**
- `salmon_quant_dir` (File): Compressed tarball containing Salmon quantification results including abundance estimates
- `output_sample_name` (String): Sample name used for output files

### `merge_results`
Merges Salmon quantification results from multiple samples into TPM and count matrices.

**Inputs:**
- `salmon_quant_dirs` (Array[File]): Array of compressed Salmon quantification directories from multiple samples
- `sample_names` (Array[String]): Array of sample names corresponding to the quantification directories
- `cpu_cores` (Int): Number of CPU cores (default: 4)
- `memory_gb` (Int): Memory allocation in GB (default: 8)

**Outputs:**
- `tpm_matrix` (File): Tab-separated matrix of TPM (Transcripts Per Million) values for all samples
- `counts_matrix` (File): Tab-separated matrix of estimated read counts for all samples
- `sample_list` (File): Text file listing all sample names included in the matrices

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-salmon/ww-salmon.wdl" as salmon_tasks

struct SalmonSample {
    String name
    File r1_fastq
    File r2_fastq
}

workflow my_rnaseq_quantification {
  input {
    Array[SalmonSample] samples
    File transcriptome_fasta
  }

  # Build Salmon index
  call salmon_tasks.build_index {
    input:
      transcriptome_fasta = transcriptome_fasta
  }

  # Quantify each sample
  scatter (sample in samples) {
    call salmon_tasks.quantify {
      input:
        salmon_index_dir = build_index.salmon_index,
        sample_name = sample.name,
        fastq_r1 = sample.r1_fastq,
        fastq_r2 = sample.r2_fastq
    }
  }

  # Merge results
  call salmon_tasks.merge_results {
    input:
      salmon_quant_dirs = quantify.salmon_quant_dir,
      sample_names = quantify.output_sample_name
  }

  output {
    Array[File] individual_quants = quantify.salmon_quant_dir
    File tpm_matrix = merge_results.tpm_matrix
    File counts_matrix = merge_results.counts_matrix
  }
}
```

### Advanced Usage Examples

**Custom resource allocation:**
```wdl
call salmon_tasks.build_index {
  input:
    transcriptome_fasta = my_transcriptome,
    cpu_cores = 16,
    memory_gb = 32
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-sra**: Download sequencing data before quantification
- **ww-testdata**: Automatic provisioning of reference data and test samples
- **ww-deseq2**: Perform differential expression analysis on merged count matrices
- **Custom workflows**: Foundation for RNA-seq quantification pipelines

## Testing the Module

The module includes a demonstration workflow that can be tested independently. The workflow in `testrun.wdl` automatically downloads test data and runs without requiring input files:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

The test workflow (`salmon_example`) automatically:
1. Downloads protein-coding transcriptome from GENCODE using `ww-testdata`
2. Downloads demonstration FASTQ data using `ww-testdata`
3. Builds Salmon index from the transcriptome
4. Performs transcript quantification using Salmon
5. Merges results into TPM and count matrices
6. Validates all outputs

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 8-32 GB recommended (scales with transcriptome size and sample complexity)
- **CPUs**: 4-8 cores recommended; Salmon benefits from multi-threading
- **Index building**: Memory requirements depend on transcriptome size
  - Human transcriptome: 16 GB recommended
  - Mouse transcriptome: 16 GB recommended
  - Smaller organisms: 8 GB may suffice

### Advanced Considerations

- **Library type**: The `quantify` task uses `--libType A` for automatic library type detection
- **Bias correction**: The task includes `--gcBias` and `--seqBias` flags for improved accuracy
- **Mapping validation**: Uses `--validateMappings` for accurate alignment-based selective alignment
- Resource allocation can be tuned for different compute environments

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support for containerized execution
- Sufficient computational resources (Salmon is generally lightweight compared to alignment-based methods)

## Features

- **Multi-sample support**: Process multiple samples in parallel using scatter-gather patterns
- **Result aggregation**: Automatically merge quantification results across samples into TPM and count matrices
- **Module integration**: Seamlessly combines with ww-sra, ww-testdata, and ww-deseq2 modules
- **Best practices**: Implements recommended Salmon parameters (--validateMappings, --gcBias, --seqBias)
- **Scalable**: Configurable resource allocation for varying dataset sizes
- **Comprehensive metadata**: Extensive WDL annotations for workflow documentation and validation
- **Cross-platform compatible**: Works with multiple WDL executors (Cromwell, miniWDL, Sprocket)

## Performance Considerations

- **Speed**: Salmon provides fast transcript-level quantification through quasi-mapping, making it suitable for large-scale RNA-seq studies
- **Memory usage**: Index building and quantification require moderate RAM (8-32GB typically)
- **CPU scaling**: Both index building and quantification benefit from multiple cores
- **Storage requirements**: Moderate disk space needed for index files and quantification results

## Output Description

- **Salmon index**: Compressed index files for rapid quantification
- **Quantification directories**: Per-sample abundance estimates including:
  - `quant.sf`: Transcript-level abundance, TPM, and counts
  - `quant.genes.sf`: Gene-level abundance estimates (if gene mapping provided)
  - Auxiliary files and metadata
- **TPM matrix**: Sample-by-transcript matrix of TPM values
- **Counts matrix**: Sample-by-transcript matrix of estimated counts
- **Sample list**: Text file documenting which samples were included

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (demonstration FASTQ for integration testing)
- Comprehensive validation of all outputs
- Integration testing with ww-testdata modules

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to Salmon usage or configuration, please refer to the [Salmon documentation](https://salmon.readthedocs.io/). Please make sure to cite their work if you use Salmon in your analyses:

Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14(4), 417-419.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
