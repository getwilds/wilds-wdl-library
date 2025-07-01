# ww-ichorcna
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for estimating tumor fraction in cell-free DNA using ichorCNA.

## Overview

This module provides reusable WDL tasks for tumor fraction estimation from shallow whole genome sequencing (WGS) of cell-free DNA (cfDNA) using ichorCNA.

Some parameters have been adjusted from the defaults to improve performance for a low expected tumor fraction (< 5%). See the relevant [wiki page](https://github.com/broadinstitute/ichorCNA/wiki/Parameter-tuning-and-settings) for this version of ichorCNA for more information.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger genomics analysis pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `ichorcna_call`, `validate_outputs`
- **Workflow**: `ichorcna_example` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/ichorcna:0.2.0`

## Tasks

### `ichorcna_call`
Estimate tumor fraction in cfDNA using ichorCNA for a single sample.

**Inputs:**
- `counts_bed` (File): Tarball of per-chromosome BED files of read counts. Used to make tumor WIG file
- `wig_gc` (File): GC-content WIG file
- `wig_map` (File): Mappability score WIG file
- `panel_of_norm_rds` (File): RDS file of median corrected depth from panel of normals
- `centromeres` (File): Text file containing Centromere locations
- `name` (String): Sample ID
- `sex` (String): User-specified: male or female
- `genome` (String): Genome build (e.g. hg19)
- `genome_style` (String): Chromosome naming convention (use UCSC if desired output is to have 'chr' string): NCBI or UCSC
- `chrs` (String): Chromosomes to analyze (default: chr 1-22, X, and Y)
- `memory_gb` (Int): Memory allocated for each task in the workflow in GB
- `cpus` (Int): Number of CPU cores allocated for each task in the workflow

**Outputs:**
- `params` (File): Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions,
- `seg` (File): Segments called by the Viterbi algorithm, including subclonal status of segments (0=clonal, 1=subclonal), and filtered fo exclude Y chromosome segments if not male,
- `genomewide_pdf` (File): Genome wide plot of data annotated for estimated copy number, tumor fraction, and ploidy for the optimal solution,
- `allgenomewide_pdf` (File): Combined PDF of all solutions,
- `correct_pdf` (File):  Genome wide correction comparisons,
- `rdata` (File): Saved R image after ichorCNA has finished. Results for all solutions will be included,
- `wig` (File): WIG file created from binned read count data within input BED files,
- `sample_name` (String): Sample ID that was processed

### `validate_outputs`
Validate that ichorCNA outputs are non-empty and generates a report.

**Inputs:**
- `params_files` (Array[File]): Array of parameter files
- `seg_files` (Array[File]): Array of segment files
- `genome_pdfs` (Array[File]): Array of genome wide plot PDFs
- `allgenome_pdfs` (Array[File]): Array of combined plot PDFs
- `correct_pdfs` (Array[File]): Array of correction comparison PDFs
- `rdata_files` (Array[File]): Array of ichorCNA RData files
- `wig_files` (Array[File]): Array of WIG files generated from BED file data
- `sample_names` (Array[String]): Array of Sample ID strings

**Outputs:**
- `report` (File): Validation report summarizing file check results

## Usage as a Module

### Importing into Your Workflow

```wdl
import "path/to/ww-ichorcna.wdl" as ichorcna_tasks

struct SampleInfo {
    String name
    File counts_bed
}

workflow my_tumor_fract_pipeline {
  input {
    File wig_gc
    File wig_map
    File panel_of_norm_rds
    File centromeres
    Array[SampleInfo] samples
    String sex
    String genome
    String genome_style
  }

  scatter (sample in samples) {
    call ichorcna_tasks.ichorcna_call {
      input:
        wig_gc = wig_gc,
        wig_map = wig_map,
        panel_of_norm_rds = panel_of_norm_rds,
        centromeres = centromeres,
        counts_bed = sample.counts_bed,
        name = sample.name,
        sex = sex,
        genome = genome,
        genome_style = genome_style,
    }
  }

  call ichorcna_tasks.validate_outputs {
    input:
      params_files = ichorcna_call.params,
      seg_files = ichorcna_call.seg,
      genome_pdfs = ichorcna_call.genomewide_pdf,
      allgenome_pdfs = ichorcna_call.allgenomewide_pdf,
      correct_pdfs = ichorcna_call.correct_pdf,
      rdata_files = ichorcna_call.rdata,
      wig_files = ichorcna_call.wig,
      sample_names = ichorcna_call.sample_name
  }

  output {
    Array[File] params = ichorcna_call.params
    Array[File] seg = ichorcna_call.seg
    Array[File] genomewide_pdf = ichorcna_call.genomewide_pdf
    Array[File] allgenomewide_pdf = ichorcna_call.allgenomewide_pdf
    Array[File] correct_pdf = ichorcna_call.correct_pdf
    Array[File] rdata = ichorcna_call.rdata
    Array[File] wig = ichorcna_call.wig
    File validation_report = validate_outputs.report
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:

- **ww-sra**: Download sequencing data prior to alignment
- **ww-bwa**: Align sequences
- **ww-bedtools**: Generate genomic interval data
- **Comprehensive genomics pipelines**: Combine with SNV/indel callers for complete variant analysis

## Testing the Module

The module includes a demonstration workflow with comprehensive testing:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-ichorcna.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-ichorcna.wdl -i inputs.json

# Using Sprocket
sprocket run ww-ichorcna.wdl inputs.json
```

### Test Input Format

```json
{
  "ichorcna_example.samples": [
    {
      "name": "sample1",
      "bam": "/path/to/sample1.bam",
      "bai": "/path/to/sample1.bam.bai"
    }
  ],
  "ichorcna_example.reference_genome": {
    "name": "hg38",
    "fasta": "/path/to/genome.fasta",
    "fasta_index": "/path/to/genome.fasta.fai"
  },
  "ichorcna_example.cpus": 8,
  "ichorcna_example.memory_gb": 16
}
```

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 16 GB recommended, may need more for very high-coverage samples
- **CPUs**: 6 cores recommended, can be adjusted based on available resources
- **Storage**: Ensure sufficient space for output files and temporary data

### Advanced Parameters

- **chrs**: Choose which chromosomes to analyze
- Resource allocation can be tuned for different compute environments

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support

## Features

- **Takes BED files as input**: Takes a tarball of BED files and genrates a WIG file from those.
- **Designed for low coverage data**: Built for ultra-low-pass whole genome sequencing (ULP-WGS, 0.1x coverage)
- **Set up for low tumor fractions**: Settings have been adjusted from default for the best performace with expected tumor fractions of < 5%

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data for validation
- Comprehensive validation of all outputs

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
