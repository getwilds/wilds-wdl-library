# ww-bwa-gatk Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline demonstrating read alignment and initial quality control using modular WILDS components.

## Overview

This pipeline combines WILDS WDL modules to create the foundation of a DNA-seq processing pipeline. It integrates the `ww-bwa` and `ww-gatk` modules to perform sequence alignment, mark duplicate reads, and recalibrate base quality scores.

This pipeline serves as both a functional workflow and a demonstration of modular WDL design patterns within the WILDS ecosystem.

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library).

## Pipeline Steps

1. **BWA Index Building** (using `ww-bwa` module):
   - Builds BWA genome index from reference FASTA
   - Generates a tarball of reference index files

2. **Alignment** (using `ww-bwa` module):
   - Performs BWA-MEM alignment for each sample
   - Generates aligned BAM files and their .BAI index files

3. **Quality Control** (using `ww-gatk` module):
   - Marks duplicate reads and recalibrates base quality scores using GATK
   - Generates BAM files with marked duplicate reads and updated base quality scores and their .BAI index files

## Module Dependencies

This pipeline imports and uses:
- **ww-bwa module**: For genome indexing and read alignment (`bwa_index`, `bwa_mem` tasks)
- **ww-gatk module**: For marking duplicate reads and recalibrating base quality scores (`mark_duplicates`, `create_sequence_dictionary`, `base_recalibrator` tasks)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for module imports
- Sufficient compute resources for BWA-MEM alignment

### Input Configuration

Create an inputs JSON file with your sample(s) and reference genome. Modify the CPU and memory requirements as needed.

```json
{
  "bwa_gatk.samples": [
    {
      "name": "sample_name",
      "reads": "/path/to/R1.fastq.gz",
      "mates": "/path/to/R2.fastq.gz"
    }
  ],
  "bwa_gatk.reference_fasta": "/path/to/reference.fasta",
  "bwa_gatk.reference_fasta_index": "/path/to/reference.fasta.fai",
  "bwa_gatk.dbsnp_vcf": "/path/to/dbsnp.vcf",
  "bwa_gatk.known_indels_vcf": ["/path/to/known1.vcf", "/path/to/known2.vcf"],
  "bwa_gatk.cpu_cores": 6,
  "bwa_gatk.memory_gb": 12
}
```

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-bwa-gatk.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-bwa-gatk.wdl -i inputs.json

# Using Sprocket
sprocket run ww-bwa-gatk.wdl inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this pipeline directly to the on-premise HPC cluster:

1. Ensure your input files are accessible by the cluster
2. Update the inputs JSON
3. Submit through the PROOF interface

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `samples` | List of BwaSample objects | Array[File] | Yes | - |
| `reference_fasta` | Reference genome FASTA file | File | Yes | - |
| `reference_fasta_index` | Index for reference genome FASTA file | File | Yes | - |
| `dbsnp_vcf` | dbSNP VCF file for known variants | File | Yes | - |
| `known_indels_vcf` | List of VCFs with known indels | Array[File] | Yes | - |
| `cpu_cores` | Number of CPU cores | Int | No | 6 |
| `memory_gb` | Memory allocation in GB | Int | No | 12 |

### BwaSample Structure

```json
{
  "name": "sample_name",
  "reads": "/path/to/R1.fastq.gz",
  "mates": "/path/to/R2.fastq.gz"
}
```

## Output Files

| Output | Description | Source Module |
|--------|-------------|---------------|
| `recalibrated_bam` | Bam files with duplicate reads marked and base qualities recalibrated for each sample | ww-gatk |
| `recalibrated_bai` | Index files for each BAM | ww-gatk |
| `duplicate_metrics` | Duplicate marking statistics for each sample | ww-gatk |
| `recalibration_report` | Base recalibration report table | ww-gatk |

## Resource Considerations

### Compute Requirements
- **Memory**: 64GB recommended for human genome alignment (16 GB for testing purposes)
- **CPUs**: 8+ cores recommended for efficient processing (8 cores for testing purposes)
- **Storage**: Sufficient space for indexed reference genome and BAM files
- **Network**: Stable internet connection for module imports

### Optimization Tips
- Use call caching to save intermediate files if the pipeline crashes partway

## Pipeline Testing

This pipeline is manually tested using a small input dataset.

The modules used in this pipeline are automatically tested as part of the WILDS WDL Library CI/CD pipeline:
- Tests run on multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Uses real but small RNA-seq datasets for efficiency
- Validates complete end-to-end functionality

## Extending the Pipeline

This pipeline can be extended by:
- Adding downstream variant calling (e.g., via the modules `ww-gatk`, `ww-delly`, `ww-manta`, `ww-smoove`)

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
