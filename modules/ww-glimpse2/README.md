# ww-glimpse2 Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for GLIMPSE2 genotype imputation from low-coverage whole genome sequencing data. This module provides tasks for reference panel preparation, phasing, and imputation.

## Overview

GLIMPSE2 is a set of tools for phasing and imputation of low-coverage whole genome sequencing (WGS) data. It enables accurate genotype imputation from sequencing data with coverage as low as 0.5x, making it highly cost-effective for large-scale genomic studies.

This module wraps the core GLIMPSE2 tools and is inspired by the [Broad Institute's GLIMPSE Imputation Pipeline](https://github.com/broadinstitute/palantir-workflows/tree/main/GlimpseImputationPipeline).

**Key Features:**
- Efficient chunking of genomic regions for parallel processing
- Support for both VCF (with genotype likelihoods) and CRAM/BAM input formats
- Binary reference panel conversion for optimized imputation
- Chunk ligation for seamless chromosome-wide results
- Concordance analysis for quality assessment

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-glimpse2.wdl` - Contains task definitions for GLIMPSE2 operations
- **Test workflow**: `testrun.wdl` - Example workflows demonstrating the imputation pipeline
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `glimpse2_chunk`

Split a genomic region into chunks for parallel imputation.

**Inputs:**
- `reference_vcf` (File): Reference panel VCF/BCF file
- `reference_vcf_index` (File): Index file for reference panel
- `genetic_map` (File): Genetic map file for the chromosome
- `region` (String): Genomic region to process (e.g., chr22 or chr22:16000000-20000000)
- `output_prefix` (String): Prefix for output files
- `window_size_cm` (Float, default=2.0): Minimal window size in centiMorgans
- `buffer_size_cm` (Float, default=0.2): Buffer size in centiMorgans
- `uniform_number_variants` (Boolean, default=false): Use uniform variants per chunk
- `cpu_cores` (Int, default=4): Number of CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `chunks_file` (File): Text file containing chunk definitions

### `glimpse2_split_reference`

Convert reference panel VCF to binary format for a specific chunk.

**Inputs:**
- `reference_vcf` (File): Reference panel VCF/BCF file
- `reference_vcf_index` (File): Index file for reference panel
- `genetic_map` (File): Genetic map file for the chromosome
- `input_region` (String): Input region from chunks file
- `output_region` (String): Output region from chunks file
- `output_prefix` (String): Prefix for output files
- `keep_monomorphic_ref_sites` (Boolean, default=false): Keep monomorphic sites
- `cpu_cores` (Int, default=4): Number of CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `reference_chunk` (File): Binary reference chunk file

### `glimpse2_phase`

Perform imputation and phasing from VCF input.

**Inputs:**
- `input_vcf` (File): Input VCF/BCF with genotype likelihoods (GL/PL field)
- `input_vcf_index` (File): Index file for input VCF
- `reference_chunk` (File): Binary reference chunk from glimpse2_split_reference
- `output_prefix` (String): Prefix for output files
- `impute_reference_only_variants` (Boolean, default=false): Impute only reference variants
- `n_burnin` (Int, default=5): Number of burn-in iterations
- `n_main` (Int, default=15): Number of main iterations
- `effective_population_size` (Int, default=15000): Effective population size
- `cpu_cores` (Int, default=4): Number of CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `imputed_chunk` (File): Imputed BCF file
- `imputed_chunk_index` (File): Index for imputed BCF

### `glimpse2_phase_cram`

Perform imputation directly from CRAM/BAM files.

**Inputs:**
- `input_cram` (File): Input CRAM or BAM file
- `input_cram_index` (File): Index file for input CRAM/BAM
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Reference genome FASTA index
- `reference_chunk` (File): Binary reference chunk
- `output_prefix` (String): Prefix for output files
- Additional parameters same as `glimpse2_phase`

**Outputs:**
- `imputed_chunk` (File): Imputed BCF file
- `imputed_chunk_index` (File): Index for imputed BCF

### `glimpse2_ligate`

Ligate multiple imputed chunks into a single file.

**Inputs:**
- `imputed_chunks` (Array[File]): Array of imputed chunk BCF files
- `imputed_chunks_indices` (Array[File]): Array of index files
- `output_prefix` (String): Prefix for output files
- `output_format` (String, default="bcf"): Output format (bcf, vcf, or vcf.gz)
- `cpu_cores` (Int, default=4): Number of CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `ligated_vcf` (File): Ligated VCF/BCF file
- `ligated_vcf_index` (File): Index file

### `glimpse2_concordance`

Compute concordance metrics between imputed and truth genotypes.

**Inputs:**
- `imputed_vcf` (File): Imputed VCF/BCF file
- `imputed_vcf_index` (File): Index file
- `truth_vcf` (File): Truth/validation VCF/BCF file
- `truth_vcf_index` (File): Index file for truth VCF
- `allele_frequencies` (File, optional): Allele frequencies for binning
- `output_prefix` (String): Prefix for output files
- `region` (String, optional): Genomic region to evaluate
- `min_val_dp` (Int, default=0): Minimum depth in validation data
- `min_val_gq` (Int, default=0): Minimum genotype quality
- `cpu_cores` (Int, default=4): Number of CPU cores
- `memory_gb` (Int, default=8): Memory in GB

**Outputs:**
- `concordance_output` (Array[File]): Concordance metrics files

### `parse_chunks_file`

Utility task to parse chunks file for parallel processing.

**Inputs:**
- `chunks_file` (File): Chunks file from glimpse2_chunk

**Outputs:**
- `input_regions` (Array[String]): Input regions
- `output_regions` (Array[String]): Output regions
- `chunk_ids` (Array[String]): Chunk identifiers

## Usage as a Module

### Basic Import

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-glimpse2/ww-glimpse2.wdl" as ww_glimpse2
```

### Complete Imputation Pipeline Example

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-glimpse2/ww-glimpse2.wdl" as ww_glimpse2

workflow imputation_pipeline {
  input {
    File reference_vcf
    File reference_vcf_index
    File genetic_map
    File input_vcf
    File input_vcf_index
    String region
    String output_prefix
  }

  # Step 1: Create chunks
  call ww_glimpse2.glimpse2_chunk {
    input:
      reference_vcf = reference_vcf,
      reference_vcf_index = reference_vcf_index,
      genetic_map = genetic_map,
      region = region,
      output_prefix = output_prefix
  }

  # Step 2: Parse chunks
  call ww_glimpse2.parse_chunks_file {
    input:
      chunks_file = glimpse2_chunk.chunks_file
  }

  # Step 3 & 4: Split reference and phase in parallel
  scatter (idx in range(length(parse_chunks_file.input_regions))) {
    call ww_glimpse2.glimpse2_split_reference {
      input:
        reference_vcf = reference_vcf,
        reference_vcf_index = reference_vcf_index,
        genetic_map = genetic_map,
        input_region = parse_chunks_file.input_regions[idx],
        output_region = parse_chunks_file.output_regions[idx],
        output_prefix = "~{output_prefix}_chunk_~{idx}"
    }

    call ww_glimpse2.glimpse2_phase {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        reference_chunk = glimpse2_split_reference.reference_chunk,
        output_prefix = "~{output_prefix}_imputed_~{idx}"
    }
  }

  # Step 5: Ligate chunks
  call ww_glimpse2.glimpse2_ligate {
    input:
      imputed_chunks = glimpse2_phase.imputed_chunk,
      imputed_chunks_indices = glimpse2_phase.imputed_chunk_index,
      output_prefix = "~{output_prefix}_final"
  }

  output {
    File imputed_vcf = glimpse2_ligate.ligated_vcf
    File imputed_vcf_index = glimpse2_ligate.ligated_vcf_index
  }
}
```

## Testing the Module

The module includes a zero-configuration test workflow (`testrun.wdl`) that automatically downloads all required test data using the `ww-testdata` module.

### Automatic Demo Mode

The test workflow automatically:
1. Downloads a reference genome region (chr22:20000000-21000000) from UCSC
2. Downloads genetic map files from the GLIMPSE repository
3. Downloads and prepares a 1000 Genomes reference panel subset
4. Downloads test BAM data (NA12878)
5. Generates genotype likelihoods from BAM at reference panel sites
6. Runs the complete GLIMPSE2 imputation pipeline
7. Outputs the final imputed VCF

### Running the Test Workflow

No inputs required! Simply run:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint glimpse2_example
```

### Custom Testing with Inputs

For testing with your own data, you can provide inputs.json:

```json
{
  "glimpse2_example.test_region": "chr22:16000000-17000000",
  "glimpse2_example.test_chromosome": "chr22",
  "glimpse2_example.output_prefix": "my_test"
}
```

## Docker Container

This module uses the `getwilds/glimpse2:2.0.0` container image, which includes:
- GLIMPSE2 tools (chunk, split_reference, phase, ligate, concordance)
- bcftools for VCF manipulation and indexing
- HTSlib for file format support

## Reference Data

### Reference Panels

GLIMPSE2 works best with large reference panels. Recommended sources:
- [1000 Genomes Project](https://www.internationalgenome.org/)
- [TOPMed](https://topmed.nhlbi.nih.gov/)
- [Haplotype Reference Consortium](http://www.haplotype-reference-consortium.org/)

### Genetic Maps

Genetic maps can be downloaded from:
- [GLIMPSE GitHub repository](https://github.com/odelaneau/GLIMPSE)
- [Eagle genetic maps](https://alkesgroup.broadinstitute.org/Eagle/)

## Parameters and Resource Requirements

### Default Resources

| Task | CPU | Memory | Typical Runtime |
|------|-----|--------|-----------------|
| glimpse2_chunk | 4 | 8 GB | Minutes |
| glimpse2_split_reference | 4 | 8 GB | Minutes per chunk |
| glimpse2_phase | 4 | 8 GB | Minutes per chunk |
| glimpse2_ligate | 4 | 8 GB | Minutes |
| glimpse2_concordance | 4 | 8 GB | Minutes |

### Resource Scaling

For larger datasets or higher coverage:
- **Higher sample counts**: Increase memory for phase tasks
- **Whole genome**: Process chromosomes in parallel
- **Large reference panels**: Increase memory for split_reference

## Citation

If you use GLIMPSE2 in your research, please cite:

> Rubinacci S, Ribeiro DM, Hofmeister RJ, Delaneau O. Efficient phasing and imputation of low-coverage sequencing data using large reference panels. Nat Genet. 2021;53(1):120-126. doi:10.1038/s41588-020-00756-0

> Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O. Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes. Nat Genet. 2023;55(7):1088-1090. doi:10.1038/s41588-023-01438-3

## Related Resources

- **[GLIMPSE Documentation](https://odelaneau.github.io/GLIMPSE/)**: Official GLIMPSE documentation
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images for WDL workflows
- **[Broad Institute GLIMPSE Pipeline](https://github.com/broadinstitute/palantir-workflows/tree/main/GlimpseImputationPipeline)**: Reference implementation
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices

## Support and Feedback

For questions or issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines
