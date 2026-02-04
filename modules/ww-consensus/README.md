# ww-consensus
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for consensus variant calling using annotated variants from multiple callers.

## Overview

This module provides a reusable WDL task for generating consensus variant calls by combining results from multiple variant callers (GATK HaplotypeCaller, GATK Mutect2, and bcftools/samtools mpileup) that have been annotated with Annovar. The consensus approach increases confidence in variant calls by requiring agreement across different calling algorithms.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger variant analysis pipelines such as the ww-leukemia pipeline.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `consensus_processing`
- **Test workflow**: `testrun.wdl` (demonstration workflow that executes all tasks)
- **Script**: `consensus-trio.R` (pulled from GitHub at runtime)
- **Container**: `rocker/tidyverse:4.4.2`

## Tasks

### `consensus_processing`
Generates consensus variant calls by combining annotated variant tables from three different callers.

**Inputs:**
- `gatk_vars` (File): Annotated variant table from GATK HaplotypeCaller
- `sam_vars` (File): Annotated variant table from samtools/bcftools
- `mutect_vars` (File): Annotated variant table from GATK Mutect2
- `base_file_name` (String): Base name for output files
- `cpu_cores` (Int): Number of CPU cores to use (default: 1)
- `memory_gb` (Int): Memory allocation in GB (default: 8)

**Outputs:**
- `consensus_tsv` (File): Tab-separated file containing consensus variant calls with evidence from all callers

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-consensus/ww-consensus.wdl" as consensus_tasks

workflow my_variant_pipeline {
  input {
    File haplotypecaller_annotated
    File bcftools_annotated
    File mutect2_annotated
    String sample_name
  }

  call consensus_tasks.consensus_processing {
    input:
      gatk_vars = haplotypecaller_annotated,
      sam_vars = bcftools_annotated,
      mutect_vars = mutect2_annotated,
      base_file_name = sample_name
  }

  output {
    File consensus_variants = consensus_processing.consensus_tsv
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-gatk**: HaplotypeCaller and Mutect2 variant calling
- **ww-bcftools**: mpileup variant calling
- **ww-annovar**: Variant annotation before consensus processing
- **ww-leukemia**: Full pipeline integration for leukemia analysis

## Testing the Module

The module includes a test workflow that demonstrates the consensus processing task:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl -i inputs-test.json

# Using miniWDL
miniwdl run testrun.wdl -i inputs-test.json

# Using Sprocket
sprocket run testrun.wdl --entrypoint consensus_example -i inputs-test.json
```

Note: Testing requires annotated variant tables from the three callers as input.

## Input Requirements

The consensus processing task expects annotated variant tables in the format produced by Annovar's `table_annovar.pl` script. These tables should contain:

- Variant coordinates (chromosome, start, end, ref, alt)
- Annotation columns from the specified protocols
- Caller-specific information (quality scores, depth, etc.)

## Output Format

The consensus TSV file contains variants that have supporting evidence from multiple callers, with columns indicating:

- Variant information (position, alleles)
- Annotations from the input tables
- Caller support indicators
- Combined quality metrics

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Annotated variant tables from HaplotypeCaller, Mutect2, and bcftools
- Sufficient computational resources (8GB RAM recommended)

## Support and Documentation

For questions about this module:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines
