# ww-jcast Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for [JCAST](https://github.com/ed-lau/jcast) (Junction Centric Alternative Splicing Translator), a tool that creates custom protein sequence databases from RNA-seq alternative splicing data for mass spectrometry proteomics analysis.

## Overview

JCAST processes alternative splicing events identified by [rMATS](http://rnaseq-mats.sourceforge.net/) and translates them into protein sequences. This enables researchers to identify unique protein isoforms in mass spectrometry experiments that arise from alternative splicing events detected in RNA-seq data.

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
| `rmats_directory` | File | Yes | - | Directory (or tarball or zip file) containing rMATS output files (e.g., SE.MATS.JC.txt, MXE.MATS.JC.txt, etc.) |
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
- `output_fasta` (File): Combined FASTA file containing all translated protein sequences from alternative splicing events
- `output_directory` (File): Tarball containing all JCAST output files, including tier-specific FASTA files (T1, T2, T3, T4, canonical, orphan)

### Output Contents

The `output_directory` tarball contains tiered FASTA files:
- **T1**: High-confidence splice variants
- **T2-T4**: Lower confidence tiers
- **canonical**: Canonical protein sequences
- **orphan**: Untranslatable variants

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

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-rmats-turbo**: Detecting differential alternative splicing events from RNA-seq data

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using Sprocket
sprocket run testrun.wdl --entrypoint jcast_example

# Using miniWDL
miniwdl run testrun.wdl

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads Ensembl reference data (human chr15 GTF and FASTA) from the JCAST repository
2. Downloads example rMATS output files from the JCAST repository
3. Runs JCAST to generate protein sequences from alternative splicing events

Note: The test uses Ensembl-format files from the JCAST repository rather than `ww-testdata.download_ref_data` because JCAST requires Ensembl GTF format.

## Docker Container

This module uses the `getwilds/jcast:0.3.5` container image from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library), which includes:
- JCAST v0.3.5
- Python 3.11
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
- A tarball (`.tar.gz`)
- A tar archive (`.tar`)
- A zip file (`.zip`)
- A directory containing the files

### Reference Files

> **Important: Ensembl GTF Required**
>
> JCAST requires **Ensembl-format GTF** annotation files, which contain the `transcript_type` attribute. GTF files from other sources (UCSC, NCBI RefSeq) use different attribute names and will cause JCAST to fail with a `KeyError: 'transcript_type'` error.
>
> - **Use**: Ensembl GTF files (e.g., from [Ensembl FTP](https://ftp.ensembl.org/pub/))
> - **Do NOT use**: UCSC or NCBI RefSeq GTF files (e.g., from `ww-testdata.download_ref_data`)
>
> For testing, use `ww-testdata.download_jcast_test_data` which provides compatible Ensembl GTF and FASTA files.

- **GTF file**: Must be an Ensembl-format GTF annotation file (with `transcript_type` attribute)
- **Genome FASTA**: Should be unmasked (no `N` characters in coding regions) for accurate translation

## Citation

If you use JCAST in your research, please cite:

> Ludwig RW, Lau E. JCAST: Sample-Specific Protein Isoform Databases for Mass Spectrometry-based Proteomics Experiments. *Software Impacts*. 2021;10:100163.


## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
