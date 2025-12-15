# ww-shapemapper Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for ShapeMapper 2, a tool for analyzing RNA structure probing data to determine nucleotide flexibility and secondary structure from SHAPE-MaP and related chemical probing experiments.

## Overview

ShapeMapper analyzes data from SHAPE (Selective 2'-Hydroxyl Acylation analyzed by Primer Extension) and related RNA structure probing methods. It processes paired-end sequencing data from chemically modified and untreated control samples to generate per-nucleotide reactivity profiles that reveal RNA secondary structure.

This module provides a standardized WDL interface to ShapeMapper 2 for integration into larger bioinformatics workflows.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-shapemapper.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `run_shapemapper`

Run ShapeMapper to analyze RNA structure probing data and generate reactivity profiles.

**Inputs:**
- `sample_name` (String): Name identifier for the sample
- `target_fa` (File): FASTA file containing the target RNA sequence(s)
- `modified_r1` (File): R1 FASTQ file from modified (chemically treated) sample
- `modified_r2` (File): R2 FASTQ file from modified (chemically treated) sample
- `untreated_r1` (File): R1 FASTQ file from untreated control sample
- `untreated_r2` (File): R2 FASTQ file from untreated control sample
- `primers_fa` (File, optional): FASTA file containing primer sequences for trimming
- `min_depth` (Int, default=5000): Minimum read depth required for calculating reactivity
- `is_amplicon` (Boolean, default=false): Set to true if data is from amplicon sequencing
- `cpu_cores` (Int, default=4): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=16): Memory allocated for the task in GB

**Outputs:**
- `shape_file` (File): SHAPE reactivity profile with per-nucleotide reactivity scores
- `log_file` (File): ShapeMapper log file with processing details and quality metrics


## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-shapemapper/ww-shapemapper.wdl" as shapemapper_tasks

struct ShapeMapperSample {
    String name
    File target_fa
    File modified_r1
    File modified_r2
    File untreated_r1
    File untreated_r2
}

workflow my_rna_structure_pipeline {
  input {
    Array[ShapeMapperSample] samples
  }

  scatter (sample in samples) {
    call shapemapper_tasks.run_shapemapper {
      input:
        sample_name = sample.name,
        target_fa = sample.target_fa,
        modified_r1 = sample.modified_r1,
        modified_r2 = sample.modified_r2,
        untreated_r1 = sample.untreated_r1,
        untreated_r2 = sample.untreated_r2
    }
  }

  output {
    Array[File] shape_profiles = run_shapemapper.shape_file
    Array[File] analysis_logs = run_shapemapper.log_file
  }
}
```

### Advanced Usage Examples

**With optional primers and amplicon mode:**
```wdl
call shapemapper_tasks.run_shapemapper {
  input:
    sample_name = "amplicon_sample",
    target_fa = target_fasta,
    modified_r1 = mod_r1,
    modified_r2 = mod_r2,
    untreated_r1 = unmod_r1,
    untreated_r2 = unmod_r2,
    primers_fa = primers,
    is_amplicon = true,
    min_depth = 10000
}
```

**Custom resource allocation for large datasets:**
```wdl
call shapemapper_tasks.run_shapemapper {
  input:
    sample_name = "large_rna_sample",
    target_fa = large_target,
    modified_r1 = mod_r1,
    modified_r2 = mod_r2,
    untreated_r1 = unmod_r1,
    untreated_r2 = unmod_r2,
    cpu_cores = 8,
    memory_gb = 32,
    min_depth = 20000
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-fastqc**: Quality control of input FASTQ files before ShapeMapper analysis
- **ww-samtools**: Post-processing of alignment files if needed
- **Other WILDS modules**: Can be combined with downstream RNA analysis tools

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint shapemapper_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Demo Mode Notes

The test workflow uses placeholder data from `ww-testdata`. For actual ShapeMapper analysis, you need:
1. A target RNA FASTA file with your sequence(s) of interest
2. Paired-end FASTQ files from chemically modified samples (e.g., SHAPE-treated)
3. Paired-end FASTQ files from untreated control samples

## Docker Container

This module uses the `getwilds/shapemapper:2.3` container image, which includes:
- ShapeMapper 2.3 release package with bundled dependencies
- Python 3 runtime environment
- All necessary system dependencies

## Citation

If you use this module in your research, please cite ShapeMapper:

> Busan S, Weeks KM. (2018) Accurate detection of chemical modifications in RNA by mutational profiling (MaP) with ShapeMapper 2. RNA, 24(2):143-148.
> DOI: [10.1261/rna.061945.117](https://doi.org/10.1261/rna.061945.117)

## ShapeMapper Background

ShapeMapper 2 is designed for:
- **SHAPE-MaP**: Analysis of SHAPE (Selective 2'-Hydroxyl Acylation analyzed by Primer Extension) experiments using mutational profiling
- **DMS-MaP**: Analysis of dimethyl sulfate probing data
- **Other chemical probing methods**: Various RNA structure probing approaches

The tool provides:
- Per-nucleotide reactivity scores
- Quality metrics and read depth information
- Support for both amplicon and whole-transcript analyses
- Automatic primer trimming and alignment

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 4 cores
- **Memory**: 16 GB
- **Runtime**: Varies based on dataset size (typically 30 minutes to several hours)

### Resource Scaling
Adjust resources based on your data:
- `cpu_cores`: More cores speed up alignment and processing steps
- `memory_gb`: Increase for very large target sequences or high-coverage datasets
- `min_depth`: Higher values require more reads but provide more confident reactivity estimates

### Recommended Settings

**For amplicon data (single RNA target, high coverage):**
```wdl
cpu_cores = 4
memory_gb = 16
min_depth = 5000
is_amplicon = true
```

**For whole transcriptome data:**
```wdl
cpu_cores = 8
memory_gb = 32
min_depth = 1000
is_amplicon = false
```

## Output Files

### SHAPE File Format
The `.shape` file contains per-nucleotide reactivity data with columns:
- Nucleotide position
- Sequence
- Reactivity score
- Standard error
- Read depth

### Log File
The `.log` file contains:
- Processing parameters
- Quality metrics
- Alignment statistics
- Warnings and errors

## Troubleshooting

**Low read depth warnings:**
- Increase sequencing depth for better coverage
- Lower `min_depth` parameter if appropriate for your application

**Primer trimming issues:**
- Ensure `primers_fa` contains correct primer sequences
- Verify primer sequences match your library preparation

**Memory errors:**
- Increase `memory_gb` parameter
- Consider splitting very large target sequences

## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[ShapeMapper Documentation](https://github.com/Weeks-UNC/shapemapper2)**: Official ShapeMapper 2 documentation
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
