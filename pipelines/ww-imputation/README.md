# ww-imputation Pipeline
[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline for genotype imputation from low-coverage whole genome sequencing data using GLIMPSE2.

## Overview

This pipeline provides a unified workflow for genotype imputation, combining the functionality of multiple WILDS modules to process CRAM/BAM files against a reference panel and produce imputed VCF files. It is inspired by the [Broad Institute's GLIMPSE Imputation Pipeline](https://github.com/broadinstitute/palantir-workflows/tree/main/GlimpseImputationPipeline) but simplifies the workflow into a single, easy-to-use pipeline.

**Key Features:**
- Single unified workflow (vs. 9 separate workflows in the Broad implementation)
- Direct CRAM/BAM input support
- Parallel processing across samples and chromosomes
- One output file per sample (all chromosomes concatenated)
- Configurable imputation parameters

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library).

## Pipeline Steps

1. **Reference Preparation** (per chromosome):
   - Chunk genomic regions for parallel processing (`glimpse2_chunk`)
   - Convert reference panel to binary format (`glimpse2_split_reference`)

2. **Imputation** (per sample × chromosome × chunk):
   - Phase and impute genotypes from CRAM/BAM (`glimpse2_phase_cram`)

3. **Ligation** (per sample × chromosome):
   - Combine imputed chunks into chromosome-level files (`glimpse2_ligate`)

4. **Concatenation** (per sample):
   - Merge all chromosomes into a single output file (`bcftools.concat`)

5. **Concordance** (optional, per sample):
   - Evaluate imputation accuracy against truth genotypes (`glimpse2_concordance`)
   - Only runs when `truth_vcf` is provided in the sample

## Module Dependencies

This pipeline imports and uses:
- **ww-glimpse2 module**: For chunking, reference splitting, phasing, and ligation
- **ww-bcftools module**: For concatenating chromosome-level files

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for module imports
- Sufficient compute resources for imputation

### Input Configuration

Create an inputs JSON file with your samples and reference data:

```json
{
  "imputation.samples": [
    {
      "sample_id": "SAMPLE001",
      "cram": "/path/to/sample001.cram",
      "cram_index": "/path/to/sample001.cram.crai",
      "truth_vcf": "/path/to/sample001.truth.vcf.gz",
      "truth_vcf_index": "/path/to/sample001.truth.vcf.gz.tbi"
    }
  ],
  "imputation.chromosomes": [
    {
      "chromosome": "chr1",
      "reference_vcf": "/path/to/1000GP_chr1.bcf",
      "reference_vcf_index": "/path/to/1000GP_chr1.bcf.csi",
      "genetic_map": "/path/to/genetic_map_chr1.txt"
    },
    {
      "chromosome": "chr2",
      "reference_vcf": "/path/to/1000GP_chr2.bcf",
      "reference_vcf_index": "/path/to/1000GP_chr2.bcf.csi",
      "genetic_map": "/path/to/genetic_map_chr2.txt"
    }
  ],
  "imputation.reference_fasta": "/path/to/GRCh38.fa",
  "imputation.reference_fasta_index": "/path/to/GRCh38.fa.fai"
}
```

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-imputation.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-imputation.wdl -i inputs.json

# Using Sprocket
sprocket run ww-imputation.wdl inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this pipeline directly to the cluster:

1. Ensure your input files are accessible by the cluster
2. Update the inputs JSON
3. Submit through the PROOF interface

## Input Parameters

### Required Inputs

| Parameter | Description | Type |
|-----------|-------------|------|
| `samples` | Array of ImputationSample objects | Array[ImputationSample] |
| `chromosomes` | Array of ChromosomeData objects | Array[ChromosomeData] |
| `reference_fasta` | Reference genome FASTA file | File |
| `reference_fasta_index` | Reference genome FASTA index (.fai) | File |

### Optional Inputs

| Parameter | Description | Type | Default |
|-----------|-------------|------|---------|
| `output_format` | Output format (bcf, vcf, vcf.gz) | String | "bcf" |
| `impute_reference_only_variants` | Only impute variants in reference panel | Boolean | false |
| `window_size_cm` | Chunk window size in centiMorgans | Float | 2.0 |
| `buffer_size_cm` | Chunk buffer size in centiMorgans | Float | 0.2 |
| `n_burnin` | MCMC burn-in iterations | Int | 5 |
| `n_main` | MCMC main iterations | Int | 15 |
| `effective_population_size` | Effective population size for HMM | Int | 15000 |
| `chunk_cpu_cores` | CPU cores for chunking tasks | Int | 4 |
| `chunk_memory_gb` | Memory (GB) for chunking tasks | Int | 8 |
| `phase_cpu_cores` | CPU cores for phasing tasks | Int | 4 |
| `phase_memory_gb` | Memory (GB) for phasing tasks | Int | 8 |
| `ligate_cpu_cores` | CPU cores for ligation tasks | Int | 4 |
| `ligate_memory_gb` | Memory (GB) for ligation tasks | Int | 16 |
| `concordance_allele_frequencies` | Allele frequencies file for concordance binning | File | (optional) |
| `concordance_min_val_dp` | Minimum depth in validation/truth data | Int | 0 |
| `concordance_min_val_gq` | Minimum genotype quality in validation/truth data | Int | 0 |
| `concordance_cpu_cores` | CPU cores for concordance tasks | Int | 4 |
| `concordance_memory_gb` | Memory (GB) for concordance tasks | Int | 8 |

### Data Structures

**ImputationSample:**
```json
{
  "sample_id": "SAMPLE001",
  "cram": "/path/to/sample.cram",
  "cram_index": "/path/to/sample.cram.crai",
  "truth_vcf": "/path/to/sample.truth.vcf.gz",
  "truth_vcf_index": "/path/to/sample.truth.vcf.gz.tbi"
}
```

Note: `truth_vcf` and `truth_vcf_index` are optional. When provided, the pipeline will compute concordance metrics comparing imputed genotypes against the truth data. This is useful for validation studies where high-confidence genotypes are available (e.g., from high-coverage sequencing or genotyping arrays).

**ChromosomeData:**
```json
{
  "chromosome": "chr1",
  "reference_vcf": "/path/to/reference_panel.bcf",
  "reference_vcf_index": "/path/to/reference_panel.bcf.csi",
  "genetic_map": "/path/to/genetic_map.txt"
}
```

## Output Files

| Output | Description |
|--------|-------------|
| `imputed_vcfs` | Array of imputed VCF/BCF files (one per sample, all chromosomes combined) |
| `imputed_vcf_indices` | Array of index files for imputed VCFs |
| `concordance_outputs` | Array of concordance metrics files (null for samples without truth VCFs) |

## Reference Data

### Reference Panels

GLIMPSE2 works best with large phased reference panels. Recommended sources:
- [1000 Genomes Project](https://www.internationalgenome.org/) - High-coverage phased data
- [TOPMed](https://topmed.nhlbi.nih.gov/) - Large diverse reference panel
- [Haplotype Reference Consortium](http://www.haplotype-reference-consortium.org/)

### Genetic Maps

Genetic maps can be downloaded from:
- [GLIMPSE GitHub repository](https://github.com/odelaneau/GLIMPSE/tree/master/maps)
- [Eagle genetic maps](https://alkesgroup.broadinstitute.org/Eagle/)

## Resource Considerations

### Compute Requirements
- **Memory**: 8-16GB per phasing task (scale with reference panel size)
- **CPUs**: 4+ cores recommended for efficient processing
- **Storage**: Sufficient space for CRAM files and output VCFs
- **Network**: Stable internet connection for module imports

### Scaling
- Samples are processed in parallel
- Chromosomes within each sample are processed in parallel
- Chunks within each chromosome are processed in parallel

## Pipeline Testing

The pipeline includes a zero-configuration test workflow (`testrun.wdl`) that automatically downloads test data and runs imputation on a small chr1 region.

```bash
# Run the test workflow
miniwdl run testrun.wdl
```

## Citation

If you use GLIMPSE2 in your research, please cite:

> Rubinacci S, Ribeiro DM, Hofmeister RJ, Delaneau O. Efficient phasing and imputation of low-coverage sequencing data using large reference panels. Nat Genet. 2021;53(1):120-126. doi:10.1038/s41588-020-00756-0

> Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O. Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes. Nat Genet. 2023;55(7):1088-1090. doi:10.1038/s41588-023-01438-3

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
