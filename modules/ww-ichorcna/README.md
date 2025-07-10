# ww-ichorcna
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for tumor fraction estimation in cfDNA using ichorCNA.

## Overview

This module provides reusable WDL tasks for estimating tumor fraction in cell-free DNA (cfDNA) using ichorCNA. It analyzes copy number alterations and ploidy to determine the proportion of tumor-derived DNA in liquid biopsy samples. The module includes built-in validation and comprehensive reporting for quality assurance.

The module integrates with `ww-bedtools` for read counting and `ww-testdata` for automatic test data provisioning when input files are not provided.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `ichorcna_call`, `validate_outputs`
- **Workflow**: `ichorcna_example` (demonstration workflow with automatic test data support)
- **Container**: `getwilds/ichorcna:0.2.0`
- **Dependencies**: Integrates with `ww-bedtools` and `ww-testdata` modules for complete workflows
- **Test Data**: Automatically downloads reference genome, ichorCNA data files, and BAM data when not provided

## Tasks

### `ichorcna_call`
Estimates cfDNA tumor fraction using ichorCNA.

**Inputs:**
- `counts_bed` (File): Tarball of per-chromosome BED files of read counts
- `wig_gc` (File): GC-content WIG file
- `wig_map` (File): Mappability score WIG file
- `panel_of_norm_rds` (File): RDS file of median corrected depth from panel of normals
- `centromeres` (File): Text file containing centromere locations
- `name` (String): Sample ID
- `sex` (String): User-specified sex (male or female)
- `genome` (String): Genome build (default: "hg38")
- `genome_style` (String): Chromosome naming convention (default: "UCSC")
- `chrs` (String): Chromosomes to analyze as R vector (default: "c(1)")
- `memory_gb` (Int): Memory allocation in GB (default: 16)
- `cpus` (Int): Number of CPU cores to use (default: 6)

**Outputs:**
- `params` (File): Final converged parameters for optimal solution
- `seg` (File): Segments called by the Viterbi algorithm
- `genomewide_pdf` (File): Genome wide plot of data annotated for estimated copy number, tumor fraction, and ploidy
- `allgenomewide_pdf` (File): Combined PDF of all solutions
- `correct_pdf` (File): Genome wide correction comparisons
- `rdata` (File): Saved R image after ichorCNA has finished
- `wig` (File): WIG file created from binned read count data
- `sample_name` (String): Sample ID that was processed

### `validate_outputs`
Validates ichorCNA outputs and generates a comprehensive report.

**Inputs:**
- `params_files` (Array[File]): Array of parameter files
- `seg_files` (Array[File]): Array of segment files
- `genome_pdfs` (Array[File]): Array of genome wide plot PDFs
- `allgenome_pdfs` (Array[File]): Array of combined plot PDFs
- `correct_pdfs` (Array[File]): Array of correction comparison PDFs
- `rdata_files` (Array[File]): Array of ichorCNA RData files
- `wig_files` (Array[File]): Array of WIG files
- `sample_names` (Array[String]): Array of sample ID strings

**Outputs:**
- `report` (File): Comprehensive validation report with statistics

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
    call ichorcna_tasks.ichorcna_call {
      input:
        counts_bed = sample.counts_bed,  # Pre-generated read counts
        wig_gc = wig_gc,
        wig_map = wig_map,
        panel_of_norm_rds = panel_of_norm_rds,
        centromeres = centromeres,
        name = sample.name,
        sex = "male",
        chrs = "c(1:22, \"X\", \"Y\")",
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
- **ww-bedtools**: Generate read count BED files from BAM files (built into demo workflow)
- **ww-testdata**: Automatic provisioning of reference data and test samples
- **Custom workflows**: Foundation for any cfDNA tumor fraction analysis pipeline

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-ichorcna.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-ichorcna.wdl -i inputs.json

# Using Sprocket
sprocket run ww-ichorcna.wdl inputs.json
```

### Automatic Demo Mode

When no samples or reference files are provided, the workflow automatically:
1. Downloads reference genome data using `ww-testdata`
2. Downloads ichorCNA-specific data files using `ww-testdata`
3. Downloads demonstration BAM data using `ww-testdata`
4. Generates read count windows using `ww-bedtools`
5. Estimates tumor fraction using ichorCNA
6. Validates all outputs

### Test Input Format

**Minimal input (uses automatic demo data):**
```json
{
  "ichorcna_example.sex": "male",
  "ichorcna_example.genome": "hg38",
  "ichorcna_example.memory_gb": 8,
  "ichorcna_example.cpus": 2,
  "ichorcna_example.chrs_list": ["chr1"],
  "ichorcna_example.chrs_vec": "c(1)"
}
```

**Full input (provide your own data):**
```json
{
  "ichorcna_example.samples": [
    {
      "name": "sample1",
      "bam": "/path/to/sample1.bam",
      "bam_index": "/path/to/sample1.bam.bai"
    }
  ],
  "ichorcna_example.bed_file": "/path/to/regions.bed",
  "ichorcna_example.reference_fasta": "/path/to/reference.fasta",
  "ichorcna_example.reference_index": "/path/to/reference.fasta.fai",
  "ichorcna_example.wig_gc": "/path/to/gc_content.wig",
  "ichorcna_example.wig_map": "/path/to/mappability.wig",
  "ichorcna_example.panel_of_norm_rds": "/path/to/panel_of_normals.rds",
  "ichorcna_example.centromeres": "/path/to/centromeres.txt",
  "ichorcna_example.sex": "male",
  "ichorcna_example.genome": "hg38",
  "ichorcna_example.memory_gb": 16,
  "ichorcna_example.cpus": 6
}
```

**Note**: You can mix and match - provide some inputs and let others use test data.

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

- **Standalone execution**: Complete workflow with automatic test data download
- **Flexible input**: Use your own data or automatic demo data
- **Tumor fraction estimation**: Estimates the proportion of tumor-derived cfDNA
- **Copy number analysis**: Identifies genomic regions with copy number alterations
- **Multiple solutions**: Evaluates different ploidy and normal contamination scenarios
- **Quality validation**: Built-in output validation and comprehensive reporting
- **Module integration**: Seamlessly combines with ww-bedtools and ww-testdata
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
- **WIG file**: Read depth data in WIG format
- **Validation report**: Comprehensive validation with file integrity checks

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real sequencing data (demonstration BAM for integration testing)
- Comprehensive validation of all outputs including R object validation
- Integration testing with ww-bedtools and ww-testdata modules
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
