# ww-leukemia
[![Project Status: Work in Progress – Initial development in progress, not yet recommended for general use.](https://getwilds.org/badges/badges/wip.svg)](https://getwilds.org/badges/#wip)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Work in Progress**: This workflow is currently under active development. While functional, the modular components are still being organized for optimal reusability within the WILDS ecosystem. The current implementation provides a complete analysis pipeline but will be refactored into importable modules in future releases.

A WILDS WDL workflow for comprehensive genomic analysis using consensus variant calling from multiple algorithms on targeted panel sequencing data, with planned support for structural variant detection.

## Overview

This workflow performs end-to-end analysis of human panel/PCR-based sequencing data with a focus on leukemia analysis. It implements a consensus variant calling approach that combines results from three different variant callers (GATK HaplotypeCaller, GATK Mutect2, and samtools/bcftools) to provide high-confidence detection of single nucleotide variants (SNVs) and small insertions/deletions (indels).

The workflow is designed for targeted sequencing panels commonly used in leukemia research and clinical diagnostics, providing comprehensive variant annotation and quality control metrics. **Planned enhancements include integration of structural variant (SV) callers** to detect larger genomic rearrangements, translocations, and copy number variations that are particularly relevant in leukemia analysis.

## Features

- **Multi-caller consensus approach**: Combines GATK HaplotypeCaller, Mutect2, and bcftools for robust SNV/indel detection
- **Planned structural variant detection**: Future integration of SV callers for translocations, CNVs, and large rearrangements
- **Comprehensive annotation**: Uses Annovar with customizable protocols for variant annotation
- **Quality control**: Includes bedtools coverage analysis and Picard hybrid selection metrics
- **Leukemia-focused**: Optimized for targeted panel sequencing in leukemia analysis with clinically relevant variant types
- **Base quality recalibration**: GATK BQSR for improved variant calling accuracy
- **Flexible reference support**: Compatible with hg19 and hg38 reference genomes

## Workflow Steps

1. **Read Alignment**: BWA-MEM alignment to reference genome
2. **Duplicate Marking**: Picard MarkDuplicates for PCR duplicate identification
3. **Base Quality Recalibration**: GATK BQSR using known variant sites
4. **Quality Control**: 
   - Bedtools coverage analysis over target regions
   - Picard hybrid selection metrics for enrichment assessment
5. **Variant Calling**: Three independent variant callers for SNVs/indels:
   - GATK HaplotypeCaller (germline variants)
   - GATK Mutect2 (somatic variants, tumor-only mode)
   - samtools/bcftools mpileup (general variant calling)
6. **Variant Annotation**: Annovar annotation with customizable protocols
7. **Consensus Analysis**: R-based consensus calling combining all three callers
8. **Planned: Structural Variant Detection**: Integration of SV callers for comprehensive genomic analysis

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Target panel BED file defining regions of interest
- Reference genome files (FASTA, indices, known variant sites)
- Paired-end FASTQ files from targeted sequencing

## Usage

### Basic Usage

1. Prepare your input JSON file with sample information and reference files:

```json
{
  "ww_leukemia.samples": [
    {
      "omics_sample_name": "SAMPLE001",
      "molecular_id": "LIB001",
      "r1_fastq": "/path/to/sample_R1.fastq.gz",
      "r2_fastq": "/path/to/sample_R2.fastq.gz"
    }
  ],
  "ww_leukemia.bed_location": "/path/to/target_regions.bed",
  "ww_leukemia.ref_name": "hg38",
  "ww_leukemia.ref_fasta": "/path/to/reference.fasta",
  "ww_leukemia.annovar_protocols": "refGene,gnomad211_exome,cosmic70",
  "ww_leukemia.annovar_operation": "g,f,f"
}
```

2. Run the workflow:

```bash
# Cromwell
java -jar cromwell.jar run ww-leukemia.wdl --inputs inputs.json --options options.json

# miniWDL
miniwdl run ww-leukemia.wdl -i inputs.json

# Sprocket
sprocket run ww-leukemia.wdl inputs.json
```

### Input Parameters

| Parameter | Description | Type | Required? |
|-----------|-------------|------|-----------|
| `samples` | Array of sample information objects | Array[SampleInfo] | Yes |
| `bed_location` | BED file defining target regions | File | Yes |
| `ref_fasta` | Reference genome FASTA file | File | Yes |
| `ref_fasta_index` | Reference FASTA index (.fai) | File | Yes |
| `ref_dict` | Reference sequence dictionary | File | Yes |
| `ref_*` | BWA index files (.amb, .ann, .bwt, .pac, .sa) | Files | Yes |
| `dbsnp_vcf` | dbSNP VCF for BQSR and annotation | File | Yes |
| `dbsnp_vcf_index` | dbSNP VCF index | File | Yes |
| `known_indels_sites_vcfs` | Known indel sites for BQSR | Array[File] | Yes |
| `known_indels_sites_indices` | Indices for known indels | Array[File] | Yes |
| `af_only_gnomad` | gnomAD allele frequencies | File | Yes |
| `af_only_gnomad_index` | gnomAD VCF index | File | Yes |
| `ref_name` | Reference genome name (hg19/hg38) | String | Yes |
| `annovar_protocols` | Annovar annotation protocols | String | Yes |
| `annovar_operation` | Annovar operations | String | Yes |

#### SampleInfo Structure

Each sample in the `samples` array should contain:
- `omics_sample_name`: Primary sample identifier
- `molecular_id`: Library/molecular identifier  
- `r1_fastq`: Path to R1 FASTQ file
- `r2_fastq`: Path to R2 FASTQ file

### Output Files

The workflow produces comprehensive outputs from each analysis step:

| Output | Description |
|--------|-------------|
| `analysis_ready_bam` | Recalibrated BAM files ready for analysis |
| `analysis_ready_bai` | Index files for recalibrated BAMs |
| `gatk_vcf` | GATK HaplotypeCaller variant calls |
| `sam_vcf` | samtools/bcftools variant calls |
| `mutect_vcf` | GATK Mutect2 variant calls |
| `*_annotated_vcf` | Annovar-annotated VCF files |
| `*_annotated_table` | Tabular variant annotations |
| `consensus_variants` | **Consensus variant calls combining all callers** |
| `panel_qc` | Coverage QC metrics over target regions |
| `picard_qc` | Hybrid selection metrics |
| `picard_qc_per_target` | Per-target coverage metrics |

## For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this workflow to the on-premise HPC cluster:

1. Clone or download this repository
2. Update `inputs.json` with your sample paths and reference files
3. Update `options.json` with your preferred output location
4. Submit through PROOF interface

**Important Notes:**
- All file paths must be accessible to the cluster (e.g., `/fh/fast/`, AWS S3)
- Reference files are available on the cluster at standard locations
- Enable call caching (`write_to_cache`, `read_from_cache`) to avoid re-processing

## Docker Containers

This workflow uses the following WILDS Docker containers:

- `getwilds/gatk:4.3.0.0` - GATK variant calling and processing
- `getwilds/bwa:0.7.17` - BWA alignment
- `getwilds/bedtools:2.31.1` - Coverage analysis
- `getwilds/bcftools:1.19` - Variant calling with samtools/bcftools
- `getwilds/annovar:{ref_name}` - Variant annotation (reference-specific)
- `getwilds/consensus:0.1.1` - Consensus variant calling in R

All containers are available on both DockerHub and GitHub Container Registry.

## Configuration Guidelines

### Reference Genome Setup

The workflow supports both hg19 and hg38. Ensure all reference files are from the same genome build:
- Use `ref_name: "hg19"` or `ref_name: "hg38"`
- All reference files must match the specified build
- Annovar container is automatically selected based on `ref_name`

### Annovar Configuration

Customize annotation protocols based on your analysis needs:
```json
{
  "annovar_protocols": "refGene,knownGene,cosmic70,gnomad211_exome,clinvar_20180603",
  "annovar_operation": "g,f,f,f,f"
}
```

Common protocols for leukemia analysis:
- `refGene`, `knownGene`: Gene-based annotation
- `cosmic70`: COSMIC cancer variants
- `gnomad211_exome`: Population frequencies
- `clinvar_20180603`: Clinical significance

## Development Status

**Current Status**: The workflow is functional and produces reliable results for leukemia panel analysis.

**Upcoming Changes**:
- Refactoring into modular components for the WILDS WDL Library
- **Integration of structural variant callers** for detection of translocations, copy number variations, and large genomic rearrangements
- Individual tasks will be separated into importable modules
- Enhanced testing and validation framework
- Additional quality control and visualization components

**Using This Workflow**: While the workflow is stable for analysis, expect structural changes as it transitions to a modular architecture. Pin to specific commits for reproducible analyses.

## Support

For questions, bugs, and feature requests:
- **Issues**: [GitHub Issues](https://github.com/getwilds/wilds-wdl-library/issues)
- **General Questions**: Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- **Documentation**: [WILDS Guide](https://getwilds.org/guide/)

## Contributing

As this workflow transitions to modular architecture, we welcome:
- Bug reports and fixes
- Performance optimizations
- Documentation improvements
- Testing with additional panel designs

See our [Contributing Guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) and the [WILDS Contributor Guide](https://getwilds.org/guide/) for details.

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by this workflow
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **Other WILDS Workflows**: Additional analysis pipelines in the WILDS ecosystem

## License

Distributed under the MIT License. See `LICENSE` for details.
