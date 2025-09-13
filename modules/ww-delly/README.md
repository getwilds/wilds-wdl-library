# ww-delly
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for structural variant calling using Delly.

## Overview

This module provides reusable WDL tasks for structural variant (SV) calling from aligned BAM files using Delly. It supports detection of all major structural variant types including deletions, duplications, inversions, translocations, and insertions. The module includes built-in validation and comprehensive reporting for quality assurance.

The module uses the latest Delly algorithms for precise breakpoint detection and supports optional region targeting and exclusion for focused analysis of specific genomic areas. It can run completely standalone with automatic test data download and alignment, or integrate with existing BAM files.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `delly_call`, `validate_outputs`
- **Workflow**: `delly_example` (demonstration workflow with automatic test data support)
- **Container**: `getwilds/delly:1.2.9`
- **Dependencies**: Integrates with `ww-sra`, `ww-bwa`, and `ww-testdata` modules for complete workflows
- **Test Data**: Automatically downloads reference genome and SRA data when not provided

## Tasks

### `delly_call`
Calls structural variants from aligned BAM files using Delly.

**Inputs:**
- `aligned_bam` (File): Input aligned BAM file containing reads for variant calling
- `aligned_bam_index` (File): Index file for the aligned BAM file
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `target_regions_bed` (File?): Optional BED file to target specific regions for calling
- `exclude_regions_bed` (File?): Optional BED file to exclude problematic regions
- `sv_type` (String): Structural variant type to call (DEL, DUP, INV, TRA, INS) or empty for all types (default: "")
- `cpu_cores` (Int): Number of CPU cores to use (default: 8)
- `memory_gb` (Int): Memory allocation in GB (default: 16)

**Outputs:**
- `vcf` (File): Structural variants in compressed VCF format
- `vcf_index` (File): Index file for the VCF output
- `summary` (File): Summary statistics and run information

### `validate_outputs`
Validates Delly outputs and generates a comprehensive report.

**Inputs:**
- `delly_vcfs` (Array[File]): Array of VCF files from Delly calling
- `delly_vcf_indices` (Array[File]): Array of index files for the VCFs

**Outputs:**
- `report` (File): Comprehensive validation report with statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-delly/ww-delly.wdl" as delly_tasks

struct DellySample {
    String name
    File bam
    File bai
}

workflow my_sv_workflow {
  input {
    Array[DellySample] samples
    File reference_fasta
    File reference_fasta_index
  }
  
  scatter (sample in samples) {
    call delly_tasks.delly_call {
      input: 
        aligned_bam = sample.bam,
        aligned_bam_index = sample.bai,
        reference_fasta = reference_fasta,
        reference_fasta_index = reference_fasta_index,
        sv_type = "DEL",  # Call only deletions
        cpu_cores = 16,
        memory_gb = 32
    }
  }
  
  output {
    Array[File] sv_calls = delly_call.vcf
    Array[File] sv_indices = delly_call.vcf_index
    Array[File] summary_stats = delly_call.summary
  }
}
```

### Advanced Usage Examples

**Calling specific SV types:**
```wdl
call delly_tasks.delly_call {
  input:
    sv_type = "DEL",  # Only deletions
    # or "DUP" for duplications, "INV" for inversions, etc.
    # Leave empty ("") for all SV types
}
```

**Targeting specific regions:**
```wdl
call delly_tasks.delly_call {
  input:
    target_regions_bed = exons.bed,
    # Focus calling on specific genomic regions
}
```

**Excluding problematic regions:**
```wdl
call delly_tasks.delly_call {
  input:
    exclude_regions_bed = centromeres_telomeres.bed,
    # Recommended for whole-genome analysis
}
```

### Integration Examples

This module pairs seamlessly with other WILDS modules:
- **ww-sra**: Download sequencing data prior to alignment (built into demo workflow)
- **ww-bwa**: Sequence alignment (built into demo workflow)  
- **ww-bcftools**: For downstream variant filtering and annotation
- **ww-annotsv**: For comprehensive structural variant annotation
- **Custom workflows**: Foundation for any structural variant analysis pipeline

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-delly.wdl

# Using miniWDL
miniwdl run ww-delly.wdl

# Using Sprocket
sprocket run ww-delly.wdl
```

### Automatic Test Data

The demonstration workflow automatically:
1. Downloads reference genome data using `ww-testdata`
2. Downloads demonstration BAM data using `ww-testdata`
3. Calls structural variants using Delly
4. Validates all outputs


## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:
- **Memory**: 8-32 GB recommended depending on genome size and analysis scope
- **CPUs**: 8-16 cores typically optimal; Delly benefits from multi-threading via OpenMP
- **Storage**: Sufficient space for BAM files, reference genome, and output VCF files

### Structural Variant Type Selection

- **DEL**: Deletions (most common, fastest to call)
- **DUP**: Duplications/tandem repeats
- **INV**: Inversions  
- **TRA**: Translocations (most computationally intensive)
- **INS**: Insertions
- **""** (empty): All SV types (comprehensive but slower)

### Region Targeting and Exclusion

- **target_regions_bed**: Focus analysis on specific regions (e.g., exons, known SV hotspots)
- **exclude_regions_bed**: Skip problematic regions (highly recommended for whole-genome analysis)
- Common exclusions: centromeres, telomeres, repetitive elements, assembly gaps


## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Input BAM files must be sorted and indexed (when providing your own data)
- Reference genome with FASTA index (when providing your own data)
- Sufficient computational resources (Delly can be memory-intensive for large genomes)

## Features

- **Standalone execution**: Complete workflow with automatic test data download
- **Comprehensive SV detection**: Supports all major structural variant types
- **Flexible filtering**: Call specific SV types or all types
- **Region targeting**: Focus calling on specific genomic regions of interest
- **Region exclusion**: Skip problematic genomic regions
- **Parallel processing**: Multi-threaded execution for improved performance via OpenMP
- **Quality validation**: Built-in output validation and statistics
- **Module integration**: Seamlessly combines with ww-sra, ww-bwa, and ww-testdata
- **Standardized output**: Compressed VCF format compatible with downstream tools
- **Detailed reporting**: Comprehensive summaries with variant counts and run parameters

## Performance Considerations

- **Memory usage**: Whole-genome analysis typically requires 16-32GB RAM
- **CPU scaling**: Performance improves with additional cores (recommend 8-16 CPUs)
- **Region exclusion**: Using exclude regions BED file significantly improves runtime and accuracy
- **SV type filtering**: Calling specific SV types reduces runtime compared to calling all types
- **OpenMP threading**: Delly automatically uses all available CPU cores via OpenMP

## Output Description

- **VCF files**: Contain structural variant calls in compressed VCF format with proper INFO and FORMAT fields
- **VCF indices**: Enable rapid random access to variant regions (.csi format)
- **Summary files**: Include variant counts, parameters used, SV type filters, and basic statistics
- **Validation report**: Comprehensive validation with file integrity checks, format validation, and overall statistics

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (demonstration BAM for integration testing)
- Comprehensive validation of all outputs including VCF format validation
- Integration testing with ww-testdata modules
- Chromosome 22 subset for efficiency during CI/CD

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to Delly usage or configuration, please refer to the documentation present in the [Delly GitHub repository](https://github.com/dellytools/delly). Please make sure to cite their work if you use Delly in your analyses:

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.
DELLY: structural variant discovery by integrated paired-end and split-read analysis.
Bioinformatics. 2012 Sep 15;28(18):i333-i339.
https://doi.org/10.1093/bioinformatics/bts378

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
