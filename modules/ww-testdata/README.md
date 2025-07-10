# ww-testdata
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

WILDS WDL module for downloading reference and test data.

## Overview

The `ww-testdata` module provides standardized tasks for downloading reference genomes, annotations, and specialized datasets used across WILDS WDL workflows. This module serves as the foundation for consistent, reproducible test data management in the WILDS ecosystem.

## Purpose

Rather than maintaining large static test datasets, `ww-testdata` enables:
- **On-demand data retrieval**: Download only the data you need, when you need it
- **Standardized datasets**: Consistent reference data across all WILDS modules
- **Modular testing**: Each module can specify its exact data requirements
- **Reproducibility**: Versioned, trackable data sources with checksums
- **Efficiency**: Avoid downloading unnecessary test data for targeted module testing

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `download_ref_data`, `download_fastq_data`, `download_bam_data`, `download_ichor_data`, `download_annotsv_vcf`, `validate_outputs`
- **Workflow**: `testdata_example` (demonstration workflow that executes all tasks)

## Tasks

### `download_ref_data`
Downloads chromosome-specific reference genome data including:
- Reference FASTA file (compressed and decompressed)
- FASTA index (.fai) file created with samtools
- Gene annotations (GTF format, chromosome-specific)
- Chromosome coverage BED file

**Supported genomes**: hg38, hg19
**Configurable**: Any chromosome (chr1, chr2, chrX, etc.)

### `download_fastq_data`
Downloads small example paired-end FASTQ files for testing sequencing analysis workflows from GATK test data.

### `download_bam_data`
Downloads and processes example BAM files for testing alignment-based workflows. This task:
- Downloads BAM data from GATK test repository
- Filters to chromosome 1 only
- Removes supplementary alignments and keeps only primary alignments
- Subsamples to 10% of reads for smaller test files
- Creates a clean, indexed BAM file suitable for testing

### `download_ichor_data`
Downloads specialized reference files for ichorCNA copy number analysis:
- GC content WIG file (500kb bins)
- Mappability WIG file (500kb bins)
- Centromere location annotations
- Panel of normals RDS file

### `download_annotsv_vcf`
Downloads test VCF files for structural variant annotation workflows from the AnnotSV repository.

### `validate_outputs`
Validates all downloaded test data files to ensure they exist and are non-empty.

**Inputs**: All 13 output files from the download tasks
**Outputs**: `report` (File): Validation summary confirming file presence and basic integrity

## Usage

### As Imported Tasks

Import specific tasks into your workflows:

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as testdata

workflow my_analysis {
  call testdata.download_ref_data { 
    input: 
      chromo = "chr22",
      version = "hg38"
  }
  
  call my_analysis_task { 
    input: 
      reference_fasta = download_ref_data.fasta,
      reference_gtf = download_ref_data.gtf
  }
}
```

### Demonstration Workflow

Run the complete demonstration workflow to download a full test dataset:

```bash
# Navigate to module directory
cd modules/ww-testdata

# Run with your preferred executor
miniwdl run ww-testdata.wdl -i inputs.json
java -jar cromwell.jar run ww-testdata.wdl --inputs inputs.json
sprocket run ww-testdata.wdl inputs.json
```

### Common Integration Patterns

**RNA-seq analysis**:
```wdl
call testdata.download_ref_data { 
  input: 
    chromo = "chr22",
    version = "hg38" 
}
call star_tasks.build_star_index { 
  input: 
    reference_fasta = download_ref_data.fasta 
}
```

**Copy number analysis**:
```wdl
call testdata.download_ichor_data { }
call ichor_tasks.run_ichor { 
  input: 
    gc_wig = download_ichor_data.wig_gc,
    map_wig = download_ichor_data.wig_map
}
```

**Variant annotation**:
```wdl
call testdata.download_annotsv_vcf { }
call annotsv_tasks.annotate_variants { 
  input: 
    test_vcf = download_annotsv_vcf.test_vcf 
}
```

**BAM processing workflows**:
```wdl
call testdata.download_bam_data { }
call my_bam_analysis { 
  input: 
    input_bam = download_bam_data.bam,
    input_bai = download_bam_data.bai
}
```

## Task Reference

### download_ref_data

**Inputs**:
- `chromo` (String): Chromosome to download (default: "chr1")
- `version` (String): Genome version (default: "hg38")
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `fasta` (File): Reference chromosome FASTA file (decompressed)
- `fasta_index` (File): Samtools FASTA index (.fai)
- `gtf` (File): Chromosome-specific gene annotations
- `bed` (File): BED file covering entire chromosome

### download_fastq_data

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `r1_fastq` (File): R1 FASTQ file for paired-end sequencing
- `r2_fastq` (File): R2 FASTQ file for paired-end sequencing

### download_bam_data

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `bam` (File): Processed BAM alignment file (chr1 only, primary alignments, 10% subsampled)
- `bai` (File): BAM index file

### download_ichor_data

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `wig_gc` (File): GC content in 500kb bins
- `wig_map` (File): Mappability in 500kb bins
- `centromeres` (File): Centromere coordinates
- `panel_of_norm_rds` (File): Panel of normals for normalization

### download_annotsv_vcf

**Inputs**:
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `test_vcf` (File): Example VCF file for testing

### validate_outputs

**Inputs**:
- `ref_fasta` (File): Reference FASTA file to validate
- `ref_fasta_index` (File): Reference FASTA index file to validate
- `ref_gtf` (File): GTF annotation file to validate
- `ref_bed` (File): BED file to validate
- `r1_fastq` (File): R1 FASTQ file to validate
- `r2_fastq` (File): R2 FASTQ file to validate
- `bam` (File): BAM file to validate
- `bai` (File): BAM index file to validate
- `ichor_gc_wig` (File): ichorCNA GC content file to validate
- `ichor_map_wig` (File): ichorCNA mapping quality file to validate
- `ichor_centromeres` (File): ichorCNA centromere locations file to validate
- `ichor_panel_of_norm_rds` (File): ichorCNA panel of normals file to validate
- `annotsv_test_vcf` (File): AnnotSV test VCF file to validate
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 2)

**Outputs**:
- `report` (File): Validation summary reporting file checks and status

## Data Sources

All reference data is downloaded from authoritative public repositories:

- **UCSC Genome Browser**: Reference genomes and annotations
- **GATK Test Data**: Example FASTQ and BAM files  
- **ichorCNA Repository**: Copy number analysis references  
- **AnnotSV Repository**: Structural variant test data

Data integrity is maintained through the use of stable URLs and version-pinned resources.

## Requirements

### Runtime Dependencies
- **Containers**: `getwilds/samtools:1.11`, `getwilds/awscli:2.27.49`
- **Tools**: samtools (for FASTA indexing and BAM processing), wget, aws CLI
- **Network**: Internet access required for data downloads

### Resource Requirements
- **CPU**: 1 core (configurable)
- **Memory**: 2-4 GB (configurable)
- **Storage**: Varies by dataset (chr22: ~50MB, chr1: ~250MB)

## Integration with WILDS Ecosystem

This module is specifically designed to support other WILDS modules:

- **ww-star**: RNA-seq alignment (requires reference FASTA + GTF)
- **ww-bwa**: DNA alignment (requires reference FASTA)
- **ww-ichorcna**: Copy number analysis (requires ichorCNA reference files)
- **ww-annotsv**: Structural variant annotation (requires test VCF)

By centralizing test data downloads, `ww-testdata` enables:
- Consistent data across all WILDS workflows
- Efficient GitHub Actions testing (download only needed data)
- Simplified module development and testing
- Clear documentation of data dependencies
- Built-in validation to ensure data integrity

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Network connectivity validation
- Cross-platform compatibility testing
- Automated output validation

## Support

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
