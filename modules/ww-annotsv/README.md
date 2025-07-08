# ww-annotsv
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for comprehensive structural variant annotation using AnnotSV.

## Overview

This module provides reusable WDL tasks for annotating structural variants using AnnotSV, a comprehensive annotation tool that integrates genomic, functional, and clinical information. AnnotSV provides detailed annotations including gene overlap, regulatory regions, pathogenicity predictions, and population frequency data.

The module is designed to be a foundational component within the WILDS ecosystem, suitable for use in larger structural variant analysis pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `annotsv_annotate`, `validate_outputs`
- **Workflow**: `annotsv_example` (demonstration workflow that executes all tasks)
- **Container**: `getwilds/annotsv:3.4.4`
- **Test Data**: Automatically downloads test VCF when no input files provided

## Tasks

### `annotsv_annotate`
Annotates structural variants with comprehensive genomic and clinical information.

**Inputs:**
- `raw_vcf` (File): Input VCF file containing structural variants to annotate
- `genome_build` (String): Reference genome build - "GRCh37" or "GRCh38" (default: "GRCh38")
- `tx_source` (String): Transcript annotation source - "ENSEMBL" or "RefSeq" (default: "RefSeq")
- `annotation_mode` (String): Annotation mode - "full" (comprehensive) or "split" (one line per SV) (default: "full")
- `include_ci` (Boolean): Include confidence intervals in breakpoint coordinates (default: true)
- `exclude_benign` (Boolean): Filter out likely benign variants from output (default: false)
- `sv_min_size` (Int): Minimum SV size in bp to consider for annotation (default: 50)
- `overlap_threshold` (Int): Minimum percentage overlap with genomic features (default: 70)
- `cpu_cores` (Int): Number of CPU cores to use (default: 4)
- `memory_gb` (Int): Memory allocation in GB (default: 8)

**Outputs:**
- `annotated_tsv` (File): Tab-delimited file with detailed annotations per SV

### `validate_outputs`
Validates AnnotSV outputs and generates comprehensive statistics.

**Inputs:**
- `annotated_tsv_files` (Array[File]): Array of annotated TSV files to validate

**Outputs:**
- `report` (File): Validation summary with structural variant annotation statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-annotsv/ww-annotsv.wdl" as annotsv_tasks

workflow my_sv_annotation_pipeline {
  input {
    Array[File] sv_vcfs
    String genome_build = "GRCh38"
    String annotation_mode = "full"
    String tx_source = "RefSeq"
  }
  
  scatter (vcf in sv_vcfs) {
    call annotsv_tasks.annotsv_annotate {
      input:
        raw_vcf = vcf,
        genome_build = genome_build,
        annotation_mode = annotation_mode,
        tx_source = tx_source,
        include_ci = true,
        sv_min_size = 50
    }
  }
  
  call annotsv_tasks.validate_outputs {
    input:
      annotated_tsv_files = annotsv_annotate.annotated_tsv
  }
  
  output {
    Array[File] annotated_variants = annotsv_annotate.annotated_tsv
    File validation_report = validate_outputs.report
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-manta**: Structural variant calling followed by annotation
- **ww-bcftools**: Variant manipulation before annotation
- **Custom pipelines**: Any workflow producing structural variant VCF files

## Testing the Module

The module includes a demonstration workflow with comprehensive testing:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-annotsv.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-annotsv.wdl -i inputs.json

# Using Sprocket
sprocket run ww-annotsv.wdl inputs.json
```

### Test Input Format

When providing your own VCF files:

```json
{
  "annotsv_example.vcfs": [
    "/path/to/structural_variants.vcf"
  ],
  "annotsv_example.genome_build": "GRCh38",
  "annotsv_example.sv_min_size": 50,
  "annotsv_example.annotation_mode": "full",
  "annotsv_example.include_ci": true,
  "annotsv_example.overlap_threshold": 70,
  "annotsv_example.cpus": 4,
  "annotsv_example.memory_gb": 8
}
```

**Note**: If no `vcfs` array is provided in the inputs, the workflow will automatically download test data using the `ww-testdata` module.

## Configuration Guidelines

### Genome Build Selection

Choose the appropriate genome build to match your variant calling:
- **GRCh38**: Recommended for new analyses (current reference)
- **GRCh37**: For compatibility with older datasets or specific requirements

### Transcript Source Selection

- **RefSeq**: More conservative annotation set, often preferred for clinical applications
- **ENSEMBL**: More comprehensive annotation set, includes more transcript variants

### Annotation Modes

- **full**: Provides comprehensive annotations with one line per annotation (recommended for detailed analysis)
- **split**: Provides one line per structural variant (useful for simpler downstream processing)

### Advanced Parameters

- `overlap_threshold`: Controls the minimum overlap percentage required for annotation (70% is conservative)
- `include_ci`: Include confidence intervals for more precise breakpoint annotation
- `exclude_benign`: Filter variants with high population frequency (set `exclude_benign` to true)
- `sv_min_size`: Minimum size threshold for structural variant consideration

## Output Format

AnnotSV generates tab-delimited files with comprehensive annotations including:

### Gene-level Annotations
- Gene symbols and identifiers
- Transcript impacts and consequences
- Protein domain information
- Expression data

### Regulatory Annotations
- Promoter and enhancer regions
- TAD (Topologically Associating Domain) boundaries
- ChromHMM functional annotations

### Clinical Annotations
- Pathogenicity predictions
- Disease associations (OMIM, ClinVar)
- Population frequency data (gnomAD-SV, 1000 Genomes)
- Dosage sensitivity scores

### Structural Information
- SV type and size classification
- Breakpoint precision and confidence intervals
- Repetitive element overlap

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Input VCF files containing structural variants with proper format (or use automatic test data)
- Sufficient computational resources (8GB RAM recommended)

## Support and Documentation

For questions about this module:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

For AnnotSV-specific questions:
- [AnnotSV Documentation](https://lbgi.fr/AnnotSV/)
- [AnnotSV GitHub Repository](https://github.com/lgmgeo/AnnotSV)
