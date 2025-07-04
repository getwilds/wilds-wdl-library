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

This module contains three primary tasks:

### `download_ref_data`
Downloads chromosome-specific reference genome data including:
- Reference FASTA file (compressed)
- FASTA index (.fai) file
- Gene annotations (GTF format)
- Chromosome coverage BED file

**Supported genomes**: hg38, hg19
**Configurable**: Any chromosome (chr1, chr2, chrX, etc.)

### `download_ichor_data`
Downloads specialized reference files for ichorCNA copy number analysis:
- GC content WIG file (500kb bins)
- Mappability WIG file (500kb bins)
- Centromere location annotations
- Panel of normals RDS file

### `download_annotsv_vcf`
Downloads test VCF files for structural variant annotation workflows.

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

## Task Reference

### download_ref_data

**Inputs**:
- `chromo` (String): Chromosome to download (default: "chr1")
- `version` (String): Genome version (default: "hg38")
- `cpu_cores` (Int): CPU allocation (default: 1)
- `memory_gb` (Int): Memory allocation (default: 4)

**Outputs**:
- `fasta` (File): Reference chromosome FASTA file
- `fasta_index` (File): Samtools FASTA index (.fai)
- `gtf` (File): Chromosome-specific gene annotations
- `bed` (File): BED file covering entire chromosome

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

## Data Sources

All reference data is downloaded from authoritative public repositories:

- **UCSC Genome Browser**: Reference genomes and annotations
- **ichorCNA Repository**: Copy number analysis references  
- **AnnotSV Repository**: Structural variant test data

Data integrity is maintained through the use of stable URLs and version-pinned resources.

## Requirements

### Runtime Dependencies
- **Container**: `getwilds/samtools:1.11`
- **Tools**: samtools (for FASTA indexing), wget
- **Network**: Internet access required for data downloads

### Resource Requirements
- **CPU**: 1 core (configurable)
- **Memory**: 4 GB (configurable)
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

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Network connectivity validation
- Cross-platform compatibility testing

## Support

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
