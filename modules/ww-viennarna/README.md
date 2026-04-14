# ww-viennarna Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for RNA secondary structure prediction and analysis using the [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) package.

## Overview

This module wraps key tools from the ViennaRNA package for predicting and analyzing RNA secondary structures. ViennaRNA provides algorithms for minimum free energy (MFE) structure prediction, partition function computation, suboptimal structure enumeration, RNA-RNA interaction prediction, and local base pair probability calculations.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-viennarna.wdl` - Contains task definitions for ViennaRNA tools
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `rnafold`

Predict minimum free energy (MFE) secondary structures and optionally compute partition function for RNA sequences.

**Inputs:**
- `input_fasta` (File): FASTA file containing RNA sequences to fold
- `partition_function` (Boolean, default=false): Compute partition function and base pair probabilities
- `no_lp` (Boolean, default=false): Disallow lonely base pairs
- `temperature` (Float, default=37.0): Temperature in degrees Celsius
- `extra_args` (String, default=""): Additional command-line arguments
- `cpu_cores` (Int, default=1): Number of CPU cores
- `memory_gb` (Int, default=4): Memory in GB

**Outputs:**
- `structure_output` (File): Predicted MFE structures in dot-bracket notation with free energies
- `postscript_plots` (Array[File]): PostScript secondary structure plots

### `rnasubopt`

Enumerate suboptimal secondary structures within an energy range of the MFE.

**Inputs:**
- `input_fasta` (File): FASTA file containing an RNA sequence
- `energy_range` (Float, default=5.0): Energy range above MFE in kcal/mol
- `no_lp` (Boolean, default=false): Disallow lonely base pairs
- `temperature` (Float, default=37.0): Temperature in degrees Celsius
- `extra_args` (String, default=""): Additional command-line arguments
- `cpu_cores` (Int, default=1): Number of CPU cores
- `memory_gb` (Int, default=4): Memory in GB

**Outputs:**
- `subopt_output` (File): Suboptimal structures in dot-bracket notation with energies

### `rnacofold`

Predict the joint secondary structure and hybridization energy of two interacting RNA molecules.

**Inputs:**
- `sequence_a` (String): First RNA sequence
- `sequence_b` (String): Second RNA sequence
- `partition_function` (Boolean, default=false): Compute partition function
- `temperature` (Float, default=37.0): Temperature in degrees Celsius
- `extra_args` (String, default=""): Additional command-line arguments
- `cpu_cores` (Int, default=1): Number of CPU cores
- `memory_gb` (Int, default=4): Memory in GB

**Outputs:**
- `cofold_output` (File): Predicted joint MFE structure with interaction energy
- `postscript_plots` (Array[File]): PostScript plots of the joint structure

### `rnaplfold`

Compute local base pair probabilities using a sliding window approach, useful for predicting accessibility of RNA regions.

**Inputs:**
- `input_fasta` (File): FASTA file containing RNA sequences
- `window_size` (Int, default=70): Window size for sliding window
- `span` (Int, default=40): Maximum base pair span
- `unpaired_length` (Int, default=4): Length of unpaired region for accessibility
- `temperature` (Float, default=37.0): Temperature in degrees Celsius
- `extra_args` (String, default=""): Additional command-line arguments
- `cpu_cores` (Int, default=1): Number of CPU cores
- `memory_gb` (Int, default=4): Memory in GB

**Outputs:**
- `accessibility_profiles` (Array[File]): Accessibility profile files with unpaired probabilities
- `dp_plots` (Array[File]): Dot plot PostScript files with local pair probabilities

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-viennarna/ww-viennarna.wdl" as viennarna_tasks

workflow rna_structure_analysis {
  input {
    File rna_sequences
  }

  call viennarna_tasks.rnafold {
    input:
      input_fasta = rna_sequences,
      partition_function = true
  }

  output {
    File structures = rnafold.structure_output
    Array[File] plots = rnafold.postscript_plots
  }
}
```

### Advanced Usage Examples

**Suboptimal structure sampling with constrained energy range:**
```wdl
call viennarna_tasks.rnasubopt {
  input:
    input_fasta = rna_sequences,
    energy_range = 2.0,
    no_lp = true
}
```

**RNA-RNA interaction prediction:**
```wdl
call viennarna_tasks.rnacofold {
  input:
    sequence_a = "GGGAAAUCCC",
    sequence_b = "GGGAUUUCCCC",
    partition_function = true
}
```

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint viennarna_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Generates test RNA sequences (tRNA, hairpin, ribozyme) using `ww-testdata`
2. Runs MFE structure prediction with partition function
3. Enumerates suboptimal structures
4. Predicts RNA-RNA interaction between two short sequences
5. Computes local base pair probabilities

## Docker Container

This module uses the `getwilds/viennarna:2.7.2` container image, which includes:
- ViennaRNA Package 2.7.2
- RNAfold, RNAsubopt, RNAcofold, RNAplfold, and other ViennaRNA tools
- All necessary dependencies for RNA structure prediction

## Citation

> Lorenz, R., Bernhart, S.H., Honer zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P.F. and Hofacker, I.L. (2011),
> "ViennaRNA Package 2.0", Algorithms for Molecular Biology, 6:26
> DOI: [10.1186/1748-7188-6-26](https://doi.org/10.1186/1748-7188-6-26)

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 1 core
- **Memory**: 4 GB
- **Runtime**: Seconds to minutes depending on sequence length and task

### Resource Scaling
- `cpu_cores`: ViennaRNA tools are primarily single-threaded; increase only for very large batch jobs
- `memory_gb`: Increase for very long sequences (>10,000 nt) or large suboptimal structure enumerations
- RNAsubopt with large energy ranges can produce very large output files

## Support and Feedback

For questions or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md)

## Related Resources

- **[ViennaRNA Documentation](https://www.tbi.univie.ac.at/RNA/documentation.html)**: Official ViennaRNA documentation
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
