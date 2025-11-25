# ww-saturation Vignette
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL vignette for analyzing saturation mutagenesis experiments using modular WILDS components.

## Overview

This vignette combines WILDS WDL modules to create a complete saturation mutagenesis analysis pipeline. It integrates the `ww-bwa` and `ww-gatk` modules to perform sequence alignment and comprehensive saturation mutagenesis analysis using GATK's AnalyzeSaturationMutagenesis tool.

Saturation mutagenesis is a powerful technique used to systematically introduce mutations across a target gene or region to understand the functional impact of every possible amino acid substitution. This workflow processes sequencing data from saturation mutagenesis experiments to generate detailed reports on variant frequencies, codon usage, amino acid changes, and coverage metrics.

This vignette serves as both a functional workflow and a demonstration of modular WDL design patterns within the WILDS ecosystem.

## Acknowledgments

This vignette was inspired by the original shell script implementation developed by [Siobhan O'Brien](https://github.com/sobrien29). Special thanks to Siobhan for pioneering this analysis approach and providing the foundation for this WDL implementation.

## Vignette Structure

This vignette is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library).

## Pipeline Steps

1. **BWA Index Building** (using `ww-bwa` module):
   - Builds BWA genome index from reference FASTA
   - Generates a tarball of reference index files
   - Performed once and reused across all samples

2. **Alignment** (using `ww-bwa` module):
   - Performs BWA-MEM alignment for each sample using paired-end Nextera sequencing reads
   - Generates sorted, indexed BAM files

3. **Saturation Mutagenesis Analysis** (using `ww-gatk` module):
   - Analyzes aligned reads using GATK AnalyzeSaturationMutagenesis
   - Generates comprehensive output tables for downstream analysis:
     - Variant counts at each position
     - Amino acid counts and fractions
     - Codon counts and fractions
     - Coverage metrics
     - Read counts

## Module Dependencies

This vignette imports and uses:
- **ww-bwa module**: For genome indexing and read alignment (`bwa_index`, `bwa_mem` tasks)
- **ww-gatk module**: For saturation mutagenesis analysis (`analyze_saturation_mutagenesis` task)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for module imports
- Sufficient compute resources for BWA-MEM alignment and GATK analysis

### Input Configuration

Create an inputs JSON file with your sample(s) and reference genome. Modify the CPU and memory requirements as needed.

```json
{
  "saturation_mutagenesis.samples": [
    {
      "name": "sample1",
      "reads": "/path/to/sample1_R1.fastq.gz",
      "mates": "/path/to/sample1_R2.fastq.gz"
    },
    {
      "name": "sample2",
      "reads": "/path/to/sample2_R1.fastq.gz",
      "mates": "/path/to/sample2_R2.fastq.gz"
    }
  ],
  "saturation_mutagenesis.reference_fasta": "/path/to/reference.fasta",
  "saturation_mutagenesis.reference_fasta_index": "/path/to/reference.fasta.fai",
  "saturation_mutagenesis.reference_dict": "/path/to/reference.dict",
  "saturation_mutagenesis.orf_range": "1-300",
  "saturation_mutagenesis.cpu_cores": 4,
  "saturation_mutagenesis.memory_gb": 16
}
```

### Running the Vignette

```bash
# Using Cromwell
java -jar cromwell.jar run ww-saturation.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-saturation.wdl -i inputs.json

# Using Sprocket
sprocket run ww-saturation.wdl inputs.json
```

### Running the Test Workflow

To test the vignette with automatically downloaded test data:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this workflow directly to the on-premise HPC cluster:

1. Ensure your input files are accessible by the cluster
2. Update the inputs JSON
3. Submit through the PROOF interface

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `samples` | List of SaturationSample objects | Array[SaturationSample] | Yes | - |
| `reference_fasta` | Reference genome FASTA file | File | Yes | - |
| `reference_fasta_index` | Index for reference genome FASTA file | File | Yes | - |
| `reference_dict` | Reference genome sequence dictionary | File | Yes | - |
| `orf_range` | Open reading frame range to analyze (e.g., '1-300') | String | Yes | - |
| `cpu_cores` | Number of CPU cores | Int | No | 4 |
| `memory_gb` | Memory allocation in GB | Int | No | 16 |

### SaturationSample Structure

```json
{
  "name": "sample_name",
  "reads": "/path/to/R1.fastq.gz",
  "mates": "/path/to/R2.fastq.gz"
}
```

**Note**: The `mates` field is optional. If omitted, the workflow will process single-end reads.

## Output Files

| Output | Description | Source Module |
|--------|-------------|---------------|
| `variant_counts` | Variant count table showing mutation frequencies at each position | ww-gatk |
| `aa_counts` | Amino acid count table showing raw counts for each amino acid change | ww-gatk |
| `aa_fractions` | Amino acid fraction table showing relative frequencies of each amino acid change | ww-gatk |
| `codon_counts` | Codon count table showing raw counts for each codon variant | ww-gatk |
| `codon_fractions` | Codon fraction table showing relative frequencies of each codon variant | ww-gatk |
| `cov_length_counts` | Coverage length count table providing read coverage statistics | ww-gatk |
| `read_counts` | Read count table summarizing total read counts | ww-gatk |
| `ref_coverage` | Reference coverage table showing per-position coverage depth | ww-gatk |

All outputs are generated as arrays, with one set of files per sample.

## Resource Considerations

### Compute Requirements
- **Memory**: 16-32GB recommended depending on reference genome size and sequencing depth
- **CPUs**: 4-8 cores recommended for efficient processing
- **Storage**: Sufficient space for:
  - BWA index files
  - Aligned BAM files
  - Analysis output tables
- **Network**: Stable internet connection for module imports

### Optimization Tips
- Use call caching to save intermediate files if the pipeline crashes partway
- For large sample sets, consider adjusting scatter parallelization settings in your workflow executor
- The BWA index is built once and reused across all samples for efficiency

## Understanding ORF Range

The `orf_range` parameter defines the open reading frame region to analyze. This should correspond to the coordinates of your gene or region of interest in the reference sequence.

**Format**: `"start-end"` (e.g., `"1-300"` for nucleotide positions 1 through 300)

**Important**: Ensure the ORF range matches your experimental design and reference sequence coordinates.

## Vignette Testing

This vignette includes a `testrun.wdl` workflow that automatically downloads test data and demonstrates the complete pipeline functionality.

The modules used in this vignette are automatically tested as part of the WILDS WDL Library CI/CD pipeline:
- Tests run on multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Uses real but small datasets for efficiency
- Validates complete end-to-end functionality

## Extending the Vignette

This vignette can be extended by:
- Adding quality control steps before alignment (e.g., using `ww-fastqc`)
- Incorporating additional filtering or normalization steps
- Adding custom visualization or statistical analysis downstream
- Integrating with variant annotation tools for functional predictions

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL vignette, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
