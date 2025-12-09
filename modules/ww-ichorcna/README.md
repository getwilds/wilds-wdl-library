# ww-ichorcna
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for tumor fraction estimation in cfDNA using ichorCNA.

## Overview

This module provides reusable WDL tasks for estimating tumor fraction in cell-free DNA (cfDNA) using ichorCNA. It analyzes copy number alterations and ploidy to determine the proportion of tumor-derived DNA in liquid biopsy samples.

Within the test workflow in `testrun.wdl`, the module integrates with HMMcopy's `readCounter` for generating WIG files and `ww-testdata` for automatic test data provisioning.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `readcounter_wig`, `ichorcna_call`, `validate_outputs`
- **Test workflow**: `testrun.wdl` (demonstration workflow with automatic test data support)
- **Containers**: `getwilds/hmmcopy:1.0.0`, `getwilds/ichorcna:0.2.0`

## Tasks

### `readcounter_wig`
Generates tumor WIG file from aligned BAM files using HMMcopy's readCounter.

**Inputs:**
- `bam_file` (File): Aligned BAM file containing reads to be analyzed
- `bam_index` (File): Index for the BAM file
- `sample_name` (String): Name of the sample being analyzed
- `chromosomes` (Array[String]): Chromosomes to include in WIG file
- `window_size` (Int): Window size in base pairs for WIG format (default: 500000)
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpus` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `wig_file` (File): WIG file created from binned read count data

### `ichorcna_call`
Estimates cfDNA tumor fraction using ichorCNA.

**Inputs:**
- `wig_tumor` (File): Tumor WIG file being analyzed
- `wig_gc` (File): GC-content WIG file
- `wig_map` (File): Mappability score WIG file
- `panel_of_norm_rds` (File): RDS file of median corrected depth from panel of normals
- `centromeres` (File): Text file containing centromere locations
- `name` (String): Sample ID
- `sex` (String): User-specified sex (male or female)
- `genome` (String): Genome build (default: "hg38")
- `genome_style` (String): Chromosome naming convention (default: "UCSC")
- `chrs` (String): Chromosomes to analyze as R vector (default: "c(1:22, 'X', 'Y')")
- `memory_gb` (Int): Memory allocation in GB (default: 16)
- `cpus` (Int): Number of CPU cores to use (default: 6)

**Outputs:**
- `params` (File): Final converged parameters for optimal solution
- `seg` (File): Segments called by the Viterbi algorithm
- `genomewide_pdf` (File): Genome wide plot of data annotated for estimated copy number, tumor fraction, and ploidy
- `allgenomewide_pdf` (File): Combined PDF of all solutions
- `correct_pdf` (File): Genome wide correction comparisons
- `rdata` (File): Saved R image after ichorCNA has finished

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-ichorcna/ww-ichorcna.wdl" as ichorcna_tasks

struct IchorSample {
    String name
    File bam
    File bam_index
}

workflow my_cfDNA_workflow {
  input {
    Array[IchorSample] samples
    File wig_gc
    File wig_map
    File panel_of_norm_rds
    File centromeres
  }
  
  scatter (sample in samples) {
    call ichorcna_tasks.readcounter_wig {
      input:
        bam_file = sample.bam,
        bam_index = sample.bam_index,
        sample_name = sample.name,
        chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"],
        window_size = 500000
    }
  
    call ichorcna_tasks.ichorcna_call {
      input:
        wig_tumor = readcounter_wig.wig_file,
        wig_gc = wig_gc,
        wig_map = wig_map,
        panel_of_norm_rds = panel_of_norm_rds,
        centromeres = centromeres,
        name = sample.name,
        sex = "male",
        chrs = "c(1:22, 'X', 'Y')",
        memory_gb = 16,
        cpus = 6
    }
  }
  
  output {
    Array[File] tumor_fractions = ichorcna_call.params
    Array[File] copy_number_segments = ichorcna_call.seg
    Array[File] plots = ichorcna_call.genomewide_pdf
  }
}
```

### Integration Examples

This module pairs seamlessly with other WILDS modules:
- **HMMcopy readCounter**: Generate WIG files from BAM files (built into workflow)
- **ww-testdata**: Automatic provisioning of reference data and test samples
- **Custom workflows**: Foundation for any cfDNA tumor fraction analysis pipeline

## Testing the Module

The module includes a test workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint ichorcna_example
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads ichorCNA-specific data files using `ww-testdata`
2. Downloads demonstration BAM data using `ww-testdata`
3. Generates WIG files from BAM files using HMMcopy's readCounter
4. Estimates tumor fraction using ichorCNA with hardcoded test parameters
5. Validates all outputs

The test workflow requires no input parameters and uses the following hardcoded settings:
- **Chromosomes**: Limited to chr1 only for fast testing
- **Sex**: Male
- **Genome**: hg38
- **Resources**: 8GB RAM, 2 CPUs
- **Sample**: Single demo sample from test data

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 16 GB recommended, may need more for very high-coverage samples
- **CPUs**: 2-6 cores recommended, can be adjusted based on available resources
- **Storage**: Ensure sufficient space for output files and temporary data

### Advanced Parameters

- **chrs**: Choose which chromosomes to analyze
- Resource allocation can be tuned for different compute environments

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Input BAM files must be sorted and indexed (when providing your own data)
- ichorCNA reference files (GC content, mappability, panel of normals, centromeres)
- Sufficient computational resources for copy number analysis

## Features

- **Tumor fraction estimation**: Estimates the proportion of tumor-derived cfDNA
- **Copy number analysis**: Identifies genomic regions with copy number alterations
- **Multiple solutions**: Evaluates different ploidy and normal contamination scenarios
- **Quality validation**: Built-in output validation and comprehensive reporting
- **HMMcopy integration**: Uses standard readCounter for WIG file generation
- **Standardized output**: Multiple output formats for downstream analysis
- **Sex-aware analysis**: Handles male/female samples appropriately

## Performance Considerations

- **Memory usage**: Copy number analysis typically requires 8-16GB RAM
- **CPU scaling**: Limited parallelization within ichorCNA; 2-6 cores sufficient
- **Chromosome selection**: Using fewer chromosomes reduces runtime for testing
- **Window size**: ichorCNA uses 500kb windows by default for optimal performance

## Output Description

- **Parameters file**: Contains tumor fraction estimates, ploidy, and model parameters
- **Segments file**: Genomic segments with copy number calls and subclonal status
- **PDF plots**: Visualization of copy number profiles and model fits
- **RData file**: Complete R workspace for further analysis
- **WIG file**: Read depth data in WIG format generated by HMMcopy readCounter

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (demonstration BAM for integration testing)
- Comprehensive validation of all outputs including R object validation
- Integration testing with ww-testdata module
- Chromosome 1 subset for efficiency during CI/CD

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to ichorCNA usage or configuration, please refer to the documentation present in the [ichorCNA GitHub repository](https://github.com/broadinstitute/ichorCNA). Please make sure to cite their work if you use ichorCNA in your analyses:

Adalsteinsson VA, Ha G, Freeman SS, et al. Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. Nat Commun. 2017;8(1):1324.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
