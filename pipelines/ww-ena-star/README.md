# ww-ena-star Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline demonstrating RNA-seq analysis from ENA download to alignment using modular WILDS components.

## Overview

This pipeline showcases how to combine WILDS WDL modules to create a complete RNA-seq processing pipeline. It demonstrates the integration of the `ww-ena` and `ww-star` modules to download sequencing data from the European Nucleotide Archive and perform high-quality alignment using STAR's two-pass methodology.

This pipeline serves as both a functional workflow and a demonstration of modular WDL design patterns within the WILDS ecosystem.

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and demonstrates:

- **Module Integration**: Combining `ww-ena` and `ww-star` modules
- **Data Flow**: Seamless passing of outputs between modules
- **Best Practices**: Modular workflow design patterns

## Pipeline Steps

1. **ENA Download** (using `ww-ena` module):
   - Downloads FASTQ files from ENA accessions using `download_files` task
   - Supports FTP and Aspera transfer protocols

2. **FASTQ Pair Extraction** (using `ww-ena` module):
   - Extracts paired-end R1 and R2 files using `extract_fastq_pairs` task
   - Identifies files by common naming patterns (_1/_2, _R1/_R2)
   - Creates standardized outputs for downstream processing

3. **STAR Index Building** (using `ww-star` module):
   - Builds STAR genome index from reference files
   - Optimizes parameters for the reference genome

4. **Two-Pass Alignment** (using `ww-star` module):
   - Performs STAR two-pass alignment for each sample
   - Generates BAM files, gene counts, and QC metrics

## Module Dependencies

This pipeline imports and uses:
- **ww-ena module**: For ENA data download and FASTQ extraction (`download_files`, `extract_fastq_pairs` tasks)
- **ww-star module**: For genome indexing and alignment (`build_index`, `align_two_pass` tasks)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for ENA downloads
- Sufficient compute resources for STAR alignment

### Input Configuration

Create an inputs JSON file with your ENA accessions and reference genome:

```json
{
  "ena_star.ena_accession_list": ["ERR000001"],
  "ena_star.ref_genome": {
    "name": "hg38",
    "fasta": "/path/to/genome.fasta",
    "gtf": "/path/to/annotation.gtf"
  },
  "ena_star.ncpu": 12,
  "ena_star.memory_gb": 64
}
```

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-ena-star.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-ena-star.wdl -i inputs.json

# Using Sprocket
sprocket run ww-ena-star.wdl inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this pipeline directly to the on-premise HPC cluster:

1. Update the inputs JSON with your desired ENA accessions
2. Ensure reference genome files are accessible to the cluster
3. Submit through the PROOF interface

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `ena_accession_list` | List of ENA accession IDs to download | Array[String] | Yes | - |
| `ref_genome` | Reference genome information | RefGenome | Yes | - |
| `sjdb_overhang` | STAR splice junction overhang | Int | No | 100 |
| `genome_sa_index_nbases` | STAR index parameter | Int | No | 14 |
| `ncpu` | Number of CPU cores | Int | No | 12 |
| `memory_gb` | Memory allocation in GB | Int | No | 64 |
| `ena_file_format` | Format of files to download from ENA | String | No | READS_FASTQ |
| `ena_protocol` | Transfer protocol (FTP or ASPERA) | String | No | FTP |

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
- **Network**: Stable internet connection for ENA downloads

### Optimization Tips
- Use call caching to avoid re-downloading ENA data
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

- **ww-ena module**: ENA download functionality
- **ww-star module**: STAR alignment functionality
- **ww-sra-star pipeline**: Similar pipeline using SRA instead of ENA
- **Other pipelines**: Additional integration examples

## Testing the Pipeline

The pipeline includes a test workflow with support for execution on multiple WDL backends that automatically runs with minimal test data:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint ena_star_example
```

The test workflow automatically:
1. Downloads a small ENA test dataset
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
