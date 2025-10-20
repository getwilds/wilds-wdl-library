# ww-leukemia
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive WILDS WDL workflow for genomic analysis using consensus variant calling, structural variant detection, and tumor fraction estimation on targeted panel sequencing data with a focus on leukemia analysis.

## Overview

This workflow performs end-to-end analysis of human sequencing data with a focus on leukemia analysis. It implements a consensus variant calling approach that combines results from three different variant callers (GATK HaplotypeCaller, GATK Mutect2, and samtools/bcftools) to provide high-confidence detection of single nucleotide variants (SNVs) and small insertions/deletions (indels).

The workflow is designed for whole genome sequencing and clinical diagnostics, providing comprehensive variant annotation, structural variant detection, tumor fraction estimation, and quality control metrics.

## Features

### Core Analysis Capabilities
- **Multi-caller consensus approach**: Combines GATK HaplotypeCaller, Mutect2, and bcftools for robust SNV/indel detection
- **Comprehensive structural variant detection**: Integrated analysis using Manta, Smoove, and Delly with AnnotSV annotation
- **Tumor fraction estimation**: ichorCNA analysis for cfDNA tumor fraction and copy number profiling
- **Advanced parallelization**: Scatter-gather processing with interval-based BAM splitting for optimal performance
- **Comprehensive annotation**: Uses Annovar with customizable protocols for variant annotation
- **Quality control**: GATK best practices including base quality recalibration and comprehensive QC metrics

### Modular Architecture
- **Built with 10 WILDS modules**: Maximum reusability and maintainability through modular design
- **GATK best practices**: Complete preprocessing pipeline with duplicate marking and base recalibration
- **Flexible reference support**: Compatible with hg19 and hg38 reference genomes (hg38 preferred)
- **CRAM input support**: Directly processes CRAM files with automatic FASTQ conversion

## Workflow Steps

### Data Preprocessing
1. **CRAM to FASTQ Conversion**: samtools extraction of paired-end reads
2. **Read Alignment**: BWA-MEM alignment with proper read group assignment
3. **Duplicate Marking**: Picard MarkDuplicates for PCR duplicate identification
4. **Base Quality Recalibration**: GATK BQSR using known variant sites
5. **Quality Control**: Comprehensive QC metrics and validation

### Variant Detection Pipeline
6. **Interval-based Parallelization**: Dynamic scatter-gather across genomic intervals
7. **Multi-caller Variant Calling**:
   - GATK HaplotypeCaller (germline variants)
   - GATK Mutect2 (somatic variants, tumor-only mode)
   - samtools/bcftools mpileup (general variant calling)
8. **Variant Annotation**: Annovar annotation with customizable protocols
9. **Consensus Analysis**: R-based consensus calling combining evidence from all three callers

### Structural Variant Analysis
10. **Multi-tool SV Detection**:
    - Manta: Breakpoint assembly and local realignment
    - Smoove: Lumpy-based detection with optimization
    - Delly: Split-read and paired-end analysis
11. **SV Annotation**: AnnotSV comprehensive structural variant annotation

### Copy Number and Tumor Fraction Analysis
12. **Read Count Generation**: Windowed read counting for copy number analysis
13. **ichorCNA Analysis**: Tumor fraction estimation and copy number profiling for cfDNA

## WILDS Modules Used

This workflow integrates the following WILDS WDL modules:

| Module | Purpose | Container |
|--------|---------|-----------|
| `ww-annotsv` | Structural variant annotation | `getwilds/annotsv` |
| `ww-annovar` | Variant annotation | `getwilds/annovar` |
| `ww-bcftools` | Variant calling and processing | `getwilds/bcftools` |
| `ww-bwa` | Read alignment and indexing | `getwilds/bwa` |
| `ww-delly` | Structural variant detection | `getwilds/delly` |
| `ww-gatk` | Variant calling and preprocessing | `getwilds/gatk` |
| `ww-ichorcna` | Tumor fraction estimation | `getwilds/ichorcna` |
| `ww-manta` | Structural variant detection | `getwilds/manta` |
| `ww-samtools` | SAM/BAM/CRAM processing | `getwilds/samtools` |
| `ww-smoove` | Structural variant detection | `getwilds/smoove` |

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Reference genome files (FASTA, indices, known variant sites)
- CRAM files from targeted sequencing
- ichorCNA reference files for tumor fraction analysis

## Usage

### Basic Usage

1. Prepare your input JSON file with sample information and reference files. We recommend using the `inputs.json` file provided in this directory as a starting point and modifying the parameters/file paths accordingly.

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
| `samples` | Array of sample information objects | Array[SampleDetails] | Yes |
| `ref_fasta` | Reference genome FASTA file | File | Yes |
| `ref_fasta_index` | Reference FASTA index (.fai) | File | Yes |
| `ref_dict` | Reference sequence dictionary | File | Yes |
| `dbsnp_vcf` | dbSNP VCF for BQSR and annotation | File | Yes |
| `dbsnp_vcf_index` | dbSNP VCF index | File | Yes |
| `known_indels_sites_vcfs` | Known indel sites for BQSR | Array[File] | Yes |
| `known_indels_sites_indices` | Indices for known indels | Array[File] | Yes |
| `af_only_gnomad` | gnomAD allele frequencies | File | Yes |
| `af_only_gnomad_index` | gnomAD VCF index | File | Yes |
| `wig_gc` | GC-content WIG file for ichorCNA | File | Yes |
| `wig_map` | Mappability WIG file for ichorCNA | File | Yes |
| `panel_of_norm_rds` | Panel of normals RDS file for ichorCNA | File | Yes |
| `centromeres` | Centromere locations file for ichorCNA | File | Yes |
| `ref_name` | Reference genome name (hg19/hg38) | String | Yes |
| `annovar_protocols` | Annovar annotation protocols | String | Yes |
| `annovar_operation` | Annovar operations | String | Yes |
| `ichorcna_chromosomes` | Array of chromosomes for ichorCNA read counting | Array[String] | Yes |
| `ichorcna_chrs_string` | R-style chromosome string for ichorCNA analysis | String | Yes |
| `scatter_count` | Number of intervals for parallelization | Int | No (default: 32) |
| `high_intensity_cpus` | CPU cores for high-intensity tasks | Int | No (default: 8) |
| `high_intensity_memory_gb` | Memory (GB) for high-intensity tasks | Int | No (default: 16) |
| `standard_cpus` | CPU cores for standard-intensity tasks | Int | No (default: 4) |
| `standard_memory_gb` | Memory (GB) for standard-intensity tasks | Int | No (default: 8) |

#### SampleDetails Structure

Each sample in the `samples` array should contain:
- `name`: Primary sample identifier
- `cramfiles`: Array of CRAM files for the sample

### Output Files

The workflow produces comprehensive outputs from each analysis step:

#### Variant Calling Outputs
| Output | Description |
|--------|-------------|
| `haplotype_vcf` | GATK HaplotypeCaller variant calls |
| `mpileup_vcf` | samtools/bcftools variant calls |
| `mutect_vcf` | GATK Mutect2 variant calls |
| `mutect_vcf_index` | Index files for Mutect2 VCF files |
| `mutect_stats` | Mutect2 statistics files |
| `*_annotated_vcf` | Annovar-annotated VCF files |
| `*_annotated_table` | Tabular variant annotations |
| `consensus_variants` | **Consensus variant calls combining all callers** |
| `gatk_wgs_metrics` | Comprehensive QC metrics |

#### Structural Variant Outputs
| Output | Description |
|--------|-------------|
| `manta_sv_vcf` | Manta structural variant calls |
| `manta_sv_vcf_index` | Index files for Manta VCF files |
| `manta_sv_annotated_tsv` | Manta SV calls annotated with AnnotSV |
| `smoove_sv_vcf` | Smoove structural variant calls |
| `smoove_sv_vcf_index` | Index files for Smoove VCF files |
| `smoove_sv_annotated_tsv` | Smoove SV calls annotated with AnnotSV |
| `delly_sv_bcf` | Delly structural variant calls |
| `delly_sv_bcf_index` | Index files for Delly BCF files |
| `delly_sv_annotated_tsv` | Delly SV calls annotated with AnnotSV |

#### Tumor Fraction and Copy Number Outputs
| Output | Description |
|--------|-------------|
| `ichorcna_params` | Final converged parameters and solutions |
| `ichorcna_seg` | Segments with subclonal status annotations |
| `ichorcna_genomewide_pdf` | Genome-wide plots with copy number annotations |
| `ichorcna_allgenomewide_pdf` | Combined PDF of all solutions |
| `ichorcna_correct_pdf` | Genome-wide correction comparisons |
| `ichorcna_rdata` | Complete R workspace with all solutions |
| `ichorcna_wig` | Binned read count WIG files |

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

### ichorCNA Chromosome Configuration

The workflow requires explicit chromosome specification for ichorCNA analysis to ensure compatibility with your sequencing data:

**For whole-genome data** (all chromosomes):
```json
{
  "ichorcna_chromosomes": ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"],
  "ichorcna_chrs_string": "c(1:22, 'X', 'Y')"
}
```

**For targeted/subset data** (e.g., chr1 only for testing):
```json
{
  "ichorcna_chromosomes": ["chr1"],
  "ichorcna_chrs_string": "c('1')"
}
```

**Notes**:
- `ichorcna_chromosomes`: Array used by `readcounter_wig` for read counting
- `ichorcna_chrs_string`: R-style string used by `ichorcna_call` for analysis
- Both must match the chromosomes present in your BAM/CRAM files
- For targeted sequencing, specify only the chromosomes covered by your panel

### Performance Optimization

The workflow includes several optimization features:
- **Scatter Count**: Adjust `scatter_count` parameter based on available resources (default: 32)
- **Parallel Processing**: Internal parallelization within GATK tasks
- **Memory Scaling**: Automatic memory allocation based on scatter count

### Resource Configuration

The workflow provides fine-grained control over computational resources through task-category-specific parameters. This allows you to optimize resource usage for different environments (production HPC vs. local testing).

#### Resource Parameters

**High-Intensity Tasks** (`high_intensity_cpus`, `high_intensity_memory_gb`):
- **BWA alignment** (`bwa_mem`): Read alignment to reference genome
- **GATK preprocessing** (`markdup_recal_metrics`): Duplicate marking and base recalibration
- **Variant calling** (`mpileup_call`): bcftools variant calling
- **Structural variant calling** (`manta_call`, `smoove_call`, `delly_call`): SV detection

**Standard-Intensity Tasks** (`standard_cpus`, `standard_memory_gb`):
- **Reference indexing** (`bwa_index`): One-time BWA index generation
- **Data preparation** (`split_intervals`, `crams_to_fastq`): Preprocessing tasks
- **Variant annotation** (`annovar_annotate`): Functional annotation (3 calls)
- **SV annotation** (`annotsv_annotate`): Structural variant annotation (3 calls)
- **Copy number analysis** (`readcounter_wig`, `ichorcna_call`): ichorCNA tasks

#### Resource Configuration Examples

**Production Environment** (default values):
```json
{
  "scatter_count": 32,
  "high_intensity_cpus": 8,
  "high_intensity_memory_gb": 16,
  "standard_cpus": 4,
  "standard_memory_gb": 8
}
```

**Testing Environment** (resource-constrained):
```json
{
  "scatter_count": 2,
  "high_intensity_cpus": 2,
  "high_intensity_memory_gb": 4,
  "standard_cpus": 2,
  "standard_memory_gb": 4
}
```

**Large-scale HPC** (maximum performance):
```json
{
  "scatter_count": 64,
  "high_intensity_cpus": 16,
  "high_intensity_memory_gb": 32,
  "standard_cpus": 8,
  "standard_memory_gb": 16
}
```

**Notes**:
- Resource parameters are **optional** - the workflow uses production defaults if not specified
- Adjust based on your compute environment's available resources
- For testing on GitHub Actions or local machines, use the testing configuration
- Higher scatter counts improve parallelization but require more concurrent resources

## Development Status

**Current Status**: The workflow is production-ready and actively maintained as part of the WILDS WDL Library.

**Recent Achievements**:
- **Complete modularization**: Built using 10 WILDS WDL modules for maximum reusability
- **Structural variant detection**: Integrated Manta, Smoove, and Delly with AnnotSV annotation
- **Tumor fraction analysis**: Complete ichorCNA integration for cfDNA analysis
- **Advanced parallelization**: Scatter-gather optimization with interval-based processing
- **CRAM support**: Direct processing of CRAM files with automatic conversion
- **Flexible resource management**: Task-category-specific CPU and memory controls for optimization across different compute environments

**Ongoing Enhancements**:
- Enhanced validation and testing framework
- Additional visualization components
- Performance optimizations for large-scale studies

## Support

For questions, bugs, and feature requests:
- **Issues**: [GitHub Issues](https://github.com/getwilds/wilds-wdl-library/issues)
- **General Questions**: Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- **Documentation**: [WILDS Guide](https://getwilds.org/guide/)

## Contributing

We welcome contributions to improve the workflow:
- Bug reports and fixes
- Performance optimizations
- Documentation improvements
- Testing with additional panel designs

See our [Contributing Guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) and the [WILDS Contributor Guide](https://getwilds.org/guide/) for details.

## Related Resources

- **[WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library)**: Complete collection of WILDS WDL modules and workflows
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by this workflow
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices

## License

Distributed under the MIT License. See `LICENSE` for details.
