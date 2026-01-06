# ww-sra-star Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline demonstrating RNA-seq analysis from SRA download to alignment using modular WILDS components.

## Overview

This pipeline showcases how to combine WILDS WDL modules to create a complete RNA-seq processing pipeline. It demonstrates the integration of the `ww-sra` and `ww-star` modules to download sequencing data from NCBI's Sequence Read Archive and perform high-quality alignment using STAR's two-pass methodology.

This pipeline serves as both a functional workflow and a demonstration of modular WDL design patterns within the WILDS ecosystem.

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and demonstrates:

- **Module Integration**: Combining `ww-sra` and `ww-star` modules
- **Data Flow**: Seamless passing of outputs between modules
- **Best Practices**: Modular workflow design patterns

## Pipeline Steps

1. **SRA Download** (using `ww-sra` module):
   - Downloads FASTQ files from SRA accessions
   - Automatically detects paired-end vs single-end reads
   - Validates successful downloads

2. **STAR Index Building** (using `ww-star` module):
   - Builds STAR genome index from reference files
   - Optimizes parameters for the reference genome

3. **Two-Pass Alignment** (using `ww-star` module):
   - Performs STAR two-pass alignment for each sample
   - Generates BAM files, gene counts, and QC metrics

## Module Dependencies

This pipeline imports and uses:
- **ww-sra module**: For SRA data download (`fastqdump` task)
- **ww-star module**: For genome indexing and alignment (`build_index`, `align_two_pass` tasks)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for SRA downloads
- Sufficient compute resources for STAR alignment

### Input Configuration

Create an inputs JSON file with your SRA accessions and reference genome:

```json
{
  "sra_star.sra_id_list": ["SRR3589956"],
  "sra_star.ref_genome": {
    "name": "hg38",
    "fasta": "/path/to/genome.fasta",
    "gtf": "/path/to/annotation.gtf"
  },
  "sra_star.ncpu": 12,
  "sra_star.memory_gb": 64
}
```

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-sra-star.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-sra-star.wdl -i inputs.json

# Using Sprocket
sprocket run ww-sra-star.wdl inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this pipeline directly to the on-premise HPC cluster:

1. Update the inputs JSON with your desired SRA accessions
2. Ensure reference genome files are accessible to the cluster
3. Submit through the PROOF interface

### Running on Cirro

This pipeline includes [Cirro](https://cirro.bio/) platform configuration files in the [.cirro/](.cirro/) directory for cloud execution. To run this pipeline on Cirro:

1. Ensure you have access to a Cirro instance
2. Follow the [Cirro pipeline documentation](https://docs.cirro.bio/pipelines/adding-pipelines/) to add this workflow to your Cirro instance
3. When launching the workflow, provide:
   - A comma-separated list of SRA accession IDs (e.g., `SRR1039508,SRR1039509,SRR1039512`)
   - Reference genome files (FASTA and GTF) from your Cirro references library
   - Optional parameters for compute resources and STAR settings
4. The configuration files in `.cirro/` define:
   - Input form fields for the Cirro web interface
   - Compute resource requirements
   - Output file handling and organization
   - Data preprocessing steps

For detailed information on configuring and using Cirro pipelines, see the [official Cirro documentation](https://docs.cirro.bio/).

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `sra_id_list` | List of SRA accession IDs to download | Array[String] | Yes | - |
| `ref_genome` | Reference genome information | RefGenome | Yes | - |
| `sjdb_overhang` | STAR splice junction overhang | Int | No | 100 |
| `genome_sa_index_nbases` | STAR index parameter | Int | No | 14 |
| `ncpu` | Number of CPU cores | Int | No | 12 |
| `memory_gb` | Memory allocation in GB | Int | No | 64 |
| `max_reads` | Maximum reads to download per sample (for testing) | Int | No | all reads |

### RefGenome Structure

```json
{
  "name": "genome_name",
  "fasta": "/path/to/genome.fasta",
  "gtf": "/path/to/annotation.gtf"
}
```

## Output Files

The pipeline produces comprehensive outputs from both modules:

| Output | Description | Source Module |
|--------|-------------|---------------|
| `star_bam` | Aligned BAM files | ww-star |
| `star_bai` | BAM index files | ww-star |
| `star_gene_counts` | Gene-level read counts | ww-star |
| `star_log_final` | STAR final logs | ww-star |
| `star_log_progress` | STAR progress logs | ww-star |
| `star_log` | STAR main logs | ww-star |
| `star_sj` | Splice junction files | ww-star |

## Resource Considerations

### Compute Requirements
- **Memory**: 64GB recommended for human genome alignment (8GB for testing purposes)
- **CPUs**: 8+ cores recommended for efficient processing (2 cores for testing purposes)
- **Storage**: Sufficient space for reference genome, FASTQ files, and outputs
- **Network**: Stable internet connection for SRA downloads

### Optimization Tips
- Use call caching to avoid re-downloading SRA data
- Consider smaller test datasets (e.g., chromosome subsets) for development
- Adjust `genome_sa_index_nbases` based on genome size

## Pipeline Development

This pipeline is automatically tested as part of the WILDS WDL Library CI/CD pipeline:
- Tests run on multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Uses real but small RNA-seq datasets for efficiency
- Ensures module integration works correctly

## Integration Patterns

This pipeline demonstrates several key WDL patterns:
- **Module Composition**: How to combine multiple modules effectively
- **Data Passing**: Seamless transfer of outputs between modules
- **Struct Usage**: Efficient organization of complex data types
- **Resource Management**: Coordinated resource allocation across modules

## Extending the Pipeline

This pipeline can be extended by:
- Adding quality control modules (e.g., FastQC, MultiQC)
- Integrating differential expression analysis
- Including additional alignment tools for comparison
- Adding downstream analysis steps

## Related WILDS Components

- **ww-sra module**: SRA download functionality
- **ww-star module**: STAR alignment functionality  
- **ww-star-deseq2 pipeline**: Extended pipeline with differential expression
- **Other pipelines**: Additional integration examples

## Testing the Pipeline

The pipeline includes a test workflow with support for execution on multiple WDL backends that automatically runs with minimal test data:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint sra_star_example
```

The test workflow automatically:
1. Downloads a small SRA test dataset
2. Builds STAR genome index from reference files
3. Performs two-pass alignment
4. Validates all outputs

The pipeline is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Small test datasets (chromosome 22 subset) for efficiency
- Comprehensive output validation
- Performance benchmarking

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
