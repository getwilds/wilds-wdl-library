# ww-sra-salmon Vignette
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL vignette demonstrating RNA-seq quantification from SRA download to transcript-level abundance estimation using modular WILDS components.

## Overview

This vignette showcases how to combine WILDS WDL modules to create a complete RNA-seq quantification pipeline. It demonstrates the integration of the `ww-sra` and `ww-salmon` modules to download sequencing data from NCBI's Sequence Read Archive and perform fast, accurate transcript quantification using Salmon's quasi-mapping approach.

This vignette serves as both a functional workflow and a demonstration of modular WDL design patterns within the WILDS ecosystem.

## Vignette Structure

This vignette is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and demonstrates:

- **Module Integration**: Combining `ww-sra` and `ww-salmon` modules
- **Data Flow**: Seamless passing of outputs between modules
- **Best Practices**: Modular workflow design patterns

## Pipeline Steps

1. **SRA Download** (using `ww-sra` module):
   - Downloads FASTQ files from SRA accessions
   - Optimizes download for paired-end data by fetching it in parallel
   - Validates successful downloads

2. **Salmon Index Building** (using `ww-salmon` module):
   - Builds Salmon index from reference transcriptome
   - Optimizes k-mer indexing for fast quasi-mapping

3. **Transcript Quantification** (using `ww-salmon` module):
   - Performs Salmon quasi-mapping and quantification for each sample
   - Generates transcript-level abundance estimates (TPM and counts)
   - Produces per-sample quantification directories

4. **Result Merging** (using `ww-salmon` module):
   - Combines quantification results across all samples
   - Creates TPM and counts matrices for downstream analysis

## Module Dependencies

This vignette imports and uses:
- **ww-sra module**: For SRA data download (`fastqdump` task)
- **ww-salmon module**: For transcriptome indexing and quantification (`build_index`, `quantify`, `merge_results` tasks)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for SRA downloads
- Sufficient compute resources for Salmon quantification
- **Paired-end RNA-seq data**: This workflow is designed for paired-end sequencing data. Ensure your SRA accessions contain paired-end reads (R1 and R2 FASTQ files)

### Input Configuration

Create an inputs JSON file with your SRA accessions and reference transcriptome:

```json
{
  "sra_salmon.sra_id_list": ["SRR3589956"],
  "sra_salmon.transcriptome_fasta": "/path/to/transcriptome.fa",
  "sra_salmon.ncpu": 8,
  "sra_salmon.memory_gb": 16,
  "sra_salmon.max_reads": 100000
}
```

**Note**: The `max_reads` parameter is optional and useful for testing with a subset of reads. Omit this parameter or set it to a higher value for production runs.

**Important**: Salmon requires a **transcriptome FASTA** file (mature mRNA sequences), not a genome FASTA. The transcriptome contains only the exonic sequences for each transcript, which Salmon uses for quasi-mapping. See [Salmon documentation](https://salmon.readthedocs.io/) for details on obtaining or building transcriptome references.

### Running the Vignette

```bash
# Using Cromwell
java -jar cromwell.jar run ww-sra-salmon.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-sra-salmon.wdl -i inputs.json

# Using Sprocket
sprocket run ww-sra-salmon.wdl inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this workflow directly to the on-premise HPC cluster:

1. Update the inputs JSON with your desired SRA accessions
2. Ensure transcriptome reference file is accessible to the cluster
3. Submit through the PROOF interface

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `sra_id_list` | List of SRA accession IDs to download | Array[String] | Yes | - |
| `transcriptome_fasta` | Reference transcriptome FASTA file | File | Yes | - |
| `ncpu` | Number of CPU cores | Int | No | 8 |
| `memory_gb` | Memory allocation in GB | Int | No | 16 |
| `max_reads` | Maximum reads to download per sample (for testing) | Int | No | all reads |

### Transcriptome Reference

The transcriptome FASTA should contain mature transcript sequences (not genomic sequences). Common sources include:

- **GENCODE**: [https://www.gencodegenes.org/](https://www.gencodegenes.org/) - Download `gencode.vXX.transcripts.fa.gz` or `gencode.vXX.pc_transcripts.fa.gz` (protein-coding only)
- **Ensembl**: [https://www.ensembl.org/](https://www.ensembl.org/) - Download cDNA FASTA files
- **RefSeq**: Extract from NCBI genome annotations using tools like `gffread`

**Note**: Ensure your transcriptome matches the annotation source you plan to use for downstream analysis.

## Output Files

The vignette produces comprehensive quantification outputs:

| Output | Description | Source Module |
|--------|-------------|---------------|
| `salmon_quant_dirs` | Per-sample quantification directories | ww-salmon |
| `tpm_matrix` | Matrix of TPM values across samples | ww-salmon |
| `counts_matrix` | Matrix of estimated counts across samples | ww-salmon |

### Quantification Directory Contents

Each Salmon quantification directory contains:
- `quant.sf`: Transcript-level abundance estimates
- `cmd_info.json`: Command and parameter information
- `lib_format_counts.json`: Library format detection results
- `aux_info/`: Auxiliary quantification information
- `logs/`: Salmon log files

## Resource Considerations

### Compute Requirements
- **Memory**: 16GB recommended for human transcriptome (8GB for testing purposes)
- **CPUs**: 8+ cores recommended for efficient processing (2 cores for testing purposes)
- **Storage**: Sufficient space for transcriptome, FASTQ files, and outputs
- **Network**: Stable internet connection for SRA downloads

### Optimization Tips
- Use call caching to avoid re-downloading SRA data
- Index building is a one-time cost - reuse indices when possible
- Consider pre-built indices for common organisms

## Vignette Development

This vignette is automatically tested as part of the WILDS WDL Library CI/CD pipeline:
- Tests run on multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Uses real but small RNA-seq datasets for efficiency
- Ensures module integration works correctly

## Integration Patterns

This vignette demonstrates several key WDL patterns:
- **Module Composition**: How to combine multiple modules effectively
- **Data Passing**: Seamless transfer of outputs between modules
- **Scatter-Gather**: Parallel processing of multiple SRA samples
- **Resource Management**: Coordinated resource allocation across modules

## Extending the Vignette

This vignette can be extended by:
- Adding quality control (e.g., FastQC, MultiQC)
- Integrating differential expression analysis (e.g., DESeq2, edgeR)
- Including multiple quantification tools for comparison

## Related WILDS Components

- **ww-sra module**: SRA download functionality
- **ww-salmon module**: Salmon quantification functionality
- **ww-sra-star vignette**: Alternative pipeline using STAR alignment
- **Other vignettes**: Additional integration examples

## Testing the Vignette

The vignette includes a test workflow with support for execution on multiple WDL backends that automatically runs with minimal test data:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint sra_salmon_example
```

The test workflow automatically:
1. Downloads a small test transcriptome from GENCODE
2. Downloads a small SRA test dataset
3. Builds Salmon index from transcriptome
4. Performs quantification
5. Merges results across samples
6. Validates all outputs

The vignette is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Small test datasets (protein-coding transcripts) for efficiency
- Comprehensive output validation
- Performance benchmarking

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL vignette, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
