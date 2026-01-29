# ww-jcast Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for [JCAST](https://github.com/ed-lau/jcast) (Junction Centric Alternative Splicing Translator), a tool that creates custom protein sequence databases from RNA-seq alternative splicing data for mass spectrometry proteomics analysis.

## Overview

JCAST processes alternative splicing events identified by [rMATS](http://rnaseq-mats.sourceforge.net/) and translates them into protein sequences. This enables researchers to identify unique protein isoforms in mass spectrometry experiments that arise from alternative splicing events detected in RNA-seq data.

**Key Features:**
- Translates alternative splicing events (SE, MXE, RI, A3SS, A5SS) into protein sequences
- Integrates with rMATS output for splicing event detection
- Supports Gaussian mixture model for read count cutoff determination
- Produces FASTA files suitable for proteomics database searches

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-jcast.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `jcast`

Translates alternative splicing events from rMATS output into protein sequences for proteomics analysis.

**Inputs:**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `rmats_directory` | File | Yes | - | Directory (or tarball) containing rMATS output files |
| `gtf_file` | File | Yes | - | Ensembl GTF annotation file for the reference genome |
| `genome_fasta` | File | Yes | - | Reference genome FASTA file (should be unmasked) |
| `output_name` | String | No | `jcast_output` | Prefix for output file names |
| `min_read_count` | Int | No | `1` | Minimum skipped junction read count for translation |
| `use_gmm` | Boolean | No | `false` | Use Gaussian mixture model for read count cutoff |
| `write_canonical` | Boolean | No | `false` | Write canonical protein sequences even if splice variants are untranslatable |
| `qvalue_min` | Float | No | `0` | Minimum rMATS FDR q-value threshold |
| `qvalue_max` | Float | No | `1` | Maximum rMATS FDR q-value threshold |
| `splice_types` | String | No | `""` | Comma-separated splice types to process (MXE,RI,SE,A3SS,A5SS) |
| `cpu_cores` | Int | No | `2` | Number of CPU cores allocated |
| `memory_gb` | Int | No | `8` | Memory allocated in GB |

**Outputs:**
- `output_fasta` (File): FASTA file containing translated protein sequences from alternative splicing events
- `output_directory` (File): Tarball containing all JCAST output files

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-jcast/ww-jcast.wdl" as jcast_tasks

workflow my_proteomics_pipeline {
  input {
    File rmats_results
    File gtf_file
    File genome_fasta
  }

  call jcast_tasks.jcast {
    input:
      rmats_directory = rmats_results,
      gtf_file = gtf_file,
      genome_fasta = genome_fasta,
      output_name = "splice_proteins"
  }

  output {
    File protein_database = jcast.output_fasta
  }
}
```

### Advanced Usage Examples

**Filter by splice type and q-value:**
```wdl
call jcast_tasks.jcast {
  input:
    rmats_directory = rmats_results,
    gtf_file = gtf_file,
    genome_fasta = genome_fasta,
    output_name = "filtered_proteins",
    splice_types = "SE,MXE",
    qvalue_min = 0,
    qvalue_max = 0.05,
    min_read_count = 10
}
```

**Use Gaussian mixture model for read count filtering:**
```wdl
call jcast_tasks.jcast {
  input:
    rmats_directory = rmats_results,
    gtf_file = gtf_file,
    genome_fasta = genome_fasta,
    output_name = "gmm_proteins",
    use_gmm = true,
    write_canonical = true
}
```

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint jcast_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads reference genome data (chr1 subset) using `ww-testdata`
2. Downloads example rMATS output files from the JCAST repository
3. Runs JCAST to generate protein sequences from alternative splicing events

## Docker Container

This module uses the `getwilds/jcast:0.3.5` container image from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library), which includes:
- JCAST v0.3.5
- Python 3.7+
- All required dependencies

## Input Requirements

### rMATS Directory

The rMATS directory should contain the following files (generated by rMATS):
- `SE.MATS.JC.txt` - Skipped exon events
- `MXE.MATS.JC.txt` - Mutually exclusive exon events
- `RI.MATS.JC.txt` - Retained intron events
- `A3SS.MATS.JC.txt` - Alternative 3' splice site events
- `A5SS.MATS.JC.txt` - Alternative 5' splice site events

The input can be provided as:
- A tarball (`.tar.gz` or `.tgz`)
- A tar archive (`.tar`)
- A zip file (`.zip`)
- A directory containing the files

### Reference Files

- **GTF file**: Should be an Ensembl-format GTF annotation file matching your reference genome
- **Genome FASTA**: Should be unmasked (no `N` characters in coding regions) for accurate translation

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 2 cores
- **Memory**: 8 GB

### Resource Scaling

For larger datasets with many splicing events:
- Increase `memory_gb` to 16-32 GB
- `cpu_cores` has minimal impact as JCAST is primarily memory-bound

## Citation

If you use JCAST in your research, please cite:

> Lau E, Lam MPY. JCAST: A Tool to Create Custom Protein Databases from RNA-seq Data for Proteomics Analysis.
> bioRxiv 2022.10.31.514609; doi: https://doi.org/10.1101/2022.10.31.514609

## About JCAST

- **Purpose**: Create custom protein databases from alternative splicing events
- **Input**: rMATS output files, GTF annotation, genome FASTA
- **Output**: FASTA file with translated protein sequences
- **Documentation**: https://ed-lau.github.io/jcast/
- **Source Code**: https://github.com/ed-lau/jcast

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
- **[rMATS](http://rnaseq-mats.sourceforge.net/)**: Alternative splicing analysis tool (upstream of JCAST)
