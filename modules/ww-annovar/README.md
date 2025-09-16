# ww-annovar

A WILDS WDL module for variant annotation using Annovar. This module provides comprehensive annotation of genetic variants with customizable protocols and operations.

## Overview

Annovar is a widely-used tool for functionally annotating genetic variants detected from high-throughput sequencing data. This module provides a WDL interface to Annovar with built-in validation and flexible parameter configuration.

## Features

- **Flexible annotation protocols**: Support for RefGene, dbSNP, ClinVar, gnomAD, COSMIC, and many other databases
- **Multiple genome builds**: Compatible with both hg19 and hg38 reference genomes
- **Automatic test data**: Downloads gnomAD test data when no input VCFs provided
- **Output validation**: Built-in validation of annotated files with summary statistics
- **Parallel processing**: Scatter execution across multiple VCF files

## Task Reference

### annovar_annotate

Annotates variants using Annovar with customizable protocols and operations.

**Inputs**:
- `vcf_to_annotate` (File): Input VCF file to be annotated
- `ref_name` (String): Reference genome build name (hg19 or hg38)
- `annovar_protocols` (String): Comma-separated list of annotation protocols
- `annovar_operation` (String): Comma-separated list of operations corresponding to protocols
- `cpu_cores` (Int, optional): Number of CPU cores to allocate (default: 2)
- `memory_gb` (Int, optional): Memory in GB to allocate (default: 8)

**Outputs**:
- `annotated_vcf` (File): VCF file with Annovar annotations added
- `annotated_table` (File): Tab-delimited table with variant annotations

### validate_outputs

Validates Annovar outputs and generates summary statistics.

**Inputs**:
- `annotated_vcf_files` (Array[File]): Array of annotated VCF files to validate
- `annotated_table_files` (Array[File]): Array of annotated table files to validate

**Outputs**:
- `report` (File): Validation summary with annotation statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-annovar/ww-annovar.wdl" as annovar_tasks

workflow my_variant_pipeline {
  input {
    Array[File] variant_calls
    String genome_build = "hg38"
  }
  
  scatter (vcf in variant_calls) {
    call annovar_tasks.annovar_annotate {
      input:
        vcf_to_annotate = vcf,
        ref_name = genome_build,
        annovar_protocols = "refGene,gnomad211_exome,clinvar_20210123",
        annovar_operation = "g,f,f"
    }
  }
  
  call annovar_tasks.validate_outputs {
    input:
      annotated_vcf_files = annovar_annotate.annotated_vcf,
      annotated_table_files = annovar_annotate.annotated_table
  }
  
  output {
    Array[File] annotated_variants = annovar_annotate.annotated_vcf
    Array[File] annotation_tables = annovar_annotate.annotated_table
    File validation_report = validate_outputs.report
  }
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-gatk**: GATK variant calling followed by annotation
- **ww-bcftools**: Variant calling and manipulation before annotation
- **Custom pipelines**: Any workflow producing variant VCF files

## Testing the Module

The module includes a demonstration workflow that automatically downloads test data and runs without requiring input files:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-annovar.wdl

# Using miniWDL
miniwdl run ww-annovar.wdl

# Using Sprocket
sprocket run ww-annovar.wdl
```

## Annotation Configuration

### Common Protocol Combinations

**Basic gene annotation**:
```json
{
  "protocols": "refGene,knownGene",
  "operation": "g,f"
}
```

**Clinical variant analysis**:
```json
{
  "protocols": "refGene,gnomad211_exome,clinvar_20180603,cosmic70",
  "operation": "g,f,f,f"
}
```

**Comprehensive annotation**:
```json
{
  "protocols": "refGene,knownGene,cosmic70,esp6500siv2_all,clinvar_20180603,gnomad211_exome",
  "operation": "g,f,f,f,f,f"
}
```

### Protocol Reference

| Protocol | Type | Description |
|----------|------|-------------|
| `refGene` | Gene | RefSeq gene definitions |
| `knownGene` | Gene | UCSC knownGene definitions |
| `cosmic70` | Filter | COSMIC cancer variants |
| `clinvar_20180603` | Filter | ClinVar clinical annotations |
| `gnomad211_exome` | Filter | gnomAD population frequencies |
| `esp6500siv2_all` | Filter | ESP6500 population data |

### Operation Types

- `g`: Gene-based annotation
- `f`: Filter-based annotation (frequency, pathogenicity, etc.)

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Sufficient computational resources (recommend 2GB+ memory per task)

## Performance Considerations

- **Memory usage**: Typically requires 2-4GB RAM per annotation task
- **Database size**: Comprehensive protocols increase runtime and memory usage
- **Protocol selection**: Choose only needed protocols to optimize performance

## Output Description

- **Annotated VCF files**: Original VCF with added INFO field annotations
- **Annotation tables**: Tab-delimited files with detailed variant information and annotations
- **Validation report**: Summary of processing success and annotation statistics

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real variant data (chromosome 1 subset for efficiency)
- Comprehensive validation of all outputs and statistics

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to Annovar usage or configuration, please refer to the documentation present in the [Annovar website](https://annovar.openbioinformatics.org/). Please make sure to cite their work if you use Annovar in your analyses:

Wang K, Li M, Hakonarson H. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Research 2010; 38:e164.
