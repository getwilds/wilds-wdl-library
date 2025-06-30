# ww-delly
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for structural variant calling using Delly.

## Overview

This module provides reusable WDL tasks for structural variant (SV) calling from aligned BAM files using Delly. It supports detection of all major structural variant types including deletions, duplications, inversions, translocations, and insertions. The module includes built-in validation and comprehensive reporting for quality assurance.

The module uses the latest Delly algorithms for precise breakpoint detection and supports optional region exclusion for problematic genomic areas like centromeres and telomeres.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `delly_call`, `validate_outputs`
- **Workflow**: `delly_example` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/delly:1.2.9`

## Tasks

### `delly_call`
Calls structural variants from aligned BAM files using Delly.

**Inputs:**
- `aligned_bam` (File): Input aligned BAM file containing reads for variant calling
- `aligned_bam_index` (File): Index file for the aligned BAM file
- `reference_fasta` (File): Reference genome FASTA file
- `reference_fasta_index` (File): Index file for the reference FASTA
- `exclude_regions_bed` (File, optional): BED file to exclude problematic regions
- `sv_type` (String): Structural variant type to call (DEL, DUP, INV, TRA, INS) or empty for all types
- `cpu_cores` (Int): Number of CPU cores to use (default: 8)
- `memory_gb` (Int): Memory allocation in GB (default: 16)

**Outputs:**
- `bcf` (File): Structural variants in BCF format
- `bcf_index` (File): Index file for the BCF output
- `summary` (File): Summary statistics and run information
- `sample_name` (String): Extracted sample name from input BAM

### `validate_outputs`
Validates Delly outputs and generates a comprehensive report.

**Inputs:**
- `delly_bcfs` (Array[File]): Array of BCF files from Delly calling
- `delly_bcf_indices` (Array[File]): Array of index files for the BCFs
- `sample_names` (Array[String]): Array of sample names for validation

**Outputs:**
- `report` (File): Comprehensive validation report with statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-delly/ww-delly.wdl" as delly_tasks

workflow my_sv_workflow {
  input {
    File sample_bam
    File sample_bai
    File reference_fasta
    File reference_fasta_index
  }
  
  call delly_tasks.delly_call {
    input: 
      aligned_bam = sample_bam,
      aligned_bam_index = sample_bai,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      sv_type = "DEL",  # Call only deletions
      cpu_cores = 16,
      memory_gb = 32
  }
  
  output {
    File sv_calls = delly_call.bcf
    File sv_index = delly_call.bcf_index
    File summary_stats = delly_call.summary
  }
}
```

### Advanced Usage

**Calling specific SV types:**
```wdl
call delly_tasks.delly_call {
  input:
    sv_type = "DEL",  # Only deletions
    # or "DUP" for duplications, "INV" for inversions, etc.
    # Leave empty ("") for all SV types
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

**Processing multiple samples:**
```wdl
scatter (sample_info in samples) {
  call delly_tasks.delly_call {
    input:
      aligned_bam = sample_info.bam,
      aligned_bam_index = sample_info.bai,
      reference_fasta = reference_genome.fasta,
      reference_fasta_index = reference_genome.fasta_index
  }
}
```

### Integration Examples

This module pairs well with other WILDS modules:
- **ww-samtools**: For BAM preprocessing and quality control
- **ww-bcftools**: For downstream variant filtering and annotation

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-delly.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-delly.wdl -i inputs.json

# Using Sprocket
sprocket run ww-delly.wdl inputs.json
```

### Test Input Format

```json
{
  "delly_example.samples": [
    {
      "name": "ERR1258306",
      "bam": "/path/to/ERR1258306.chr1.aligned.bam",
      "bai": "/path/to/ERR1258306.chr1.aligned.bam.bai"
    }
  ],
  "delly_example.reference_genome": {
    "name": "chr1",
    "fasta": "/path/to/chr1.fa",
    "fasta_index": "/path/to/chr1.fa.fai"
  },
  "delly_example.sv_type": "",
  "delly_example.cpus": 2,
  "delly_example.memory_gb": 4
}
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Input BAM files must be sorted and indexed
- Reference genome with FASTA index
- Sufficient computational resources (Delly can be memory-intensive for large genomes)

## Features

- **Comprehensive SV detection**: Supports all major structural variant types
- **Flexible filtering**: Call specific SV types or all types
- **Region exclusion**: Skip problematic genomic regions
- **Parallel processing**: Multi-threaded execution for improved performance
- **Quality validation**: Built-in output validation and statistics
- **Standardized output**: BCF format compatible with downstream tools
- **Detailed reporting**: Comprehensive summaries with variant counts

## Performance Considerations

- **Memory usage**: Whole-genome analysis typically requires 16-32GB RAM
- **CPU scaling**: Performance improves with additional cores (recommend 8-16 CPUs)
- **Region exclusion**: Using exclude regions BED file significantly improves runtime and accuracy
- **SV type filtering**: Calling specific SV types reduces runtime compared to calling all types

## Output Description

- **BCF files**: Contain structural variant calls in compressed binary VCF format
- **BCF indices**: Enable rapid random access to variant regions
- **Summary files**: Include variant counts, parameters used, and basic statistics
- **Validation report**: Comprehensive validation with file integrity checks and overall statistics

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (chromosome 1 subset for efficiency)
- Comprehensive validation of all outputs and statistics

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
