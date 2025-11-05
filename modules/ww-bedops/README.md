# ww-bedops Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for genomic file format conversions and operations using bedops.

## Overview

This module provides WDL tasks for using bedops, a high-performance toolkit for genomic interval operations and format conversions. Bedops is widely used in genomics for converting between different annotation formats and performing set operations on genomic coordinates.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-bedops.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `gtf_to_bed`

Convert GTF annotation files to BED12 format using bedops.

**Inputs:**
- `gtf_file` (File): GTF annotation file to convert
- `cpu_cores` (Int, default=1): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=2): Memory allocated for the task in GB

**Outputs:**
- `bed_file` (File): BED12-formatted annotation file

**Why BED12 format?**
BED12 format includes exon/intron structure, which is essential for:
- RNA-seq QC tools like RSeQC
- Visualization in genome browsers
- Gene body coverage analysis
- Splice junction analysis

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bedops/ww-bedops.wdl" as bedops_tasks

workflow my_analysis {
  input {
    File gtf_annotation
  }

  # Convert GTF to BED format
  call bedops_tasks.gtf_to_bed {
    input:
      gtf_file = gtf_annotation
  }

  output {
    File bed_annotation = gtf_to_bed.bed_file
  }
}
```

### Integration with Other WILDS Modules

**With ww-rseqc:**
```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bedops/ww-bedops.wdl" as bedops_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-rseqc/ww-rseqc.wdl" as rseqc_tasks

workflow rnaseq_qc {
  input {
    File bam_file
    File bam_index
    File gtf_file
    String sample_name
  }

  # Convert GTF to BED for RSeQC
  call bedops_tasks.gtf_to_bed {
    input:
      gtf_file = gtf_file
  }

  # Run RSeQC QC metrics
  call rseqc_tasks.run_rseqc {
    input:
      bam_file = bam_file,
      bam_index = bam_index,
      ref_bed = gtf_to_bed.bed_file,
      sample_name = sample_name
  }

  output {
    File qc_summary = run_rseqc.rseqc_summary
  }
}
```

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using Sprocket (recommended)
sprocket run testrun.wdl --entrypoint bedops_example

# Using miniWDL
miniwdl run testrun.wdl

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads test GTF data using `ww-testdata`
2. Converts GTF to BED format using `gtf_to_bed`
3. Validates the output
4. Demonstrates the module's functionality

## Docker Container

This module uses the **`getwilds/bedops:2.4.42`** container image from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library), which includes:
- bedops version 2.4.42
- All bedops utilities (`gtf2bed`, `gff2bed`, `vcf2bed`, `bedops`, `bedmap`, etc.)
- Optimized for genomic data processing

## Citation

If you use this module in your research, please cite bedops:

> **BEDOPS: high-performance genomic feature operations**
> Shane Neph, M. Scott Kuehn, Alex P. Reynolds, et al.
> Bioinformatics. 2012 Jul 15;28(14):1919-20.
> DOI: [10.1093/bioinformatics/bts277](https://doi.org/10.1093/bioinformatics/bts277)

## Bedops Capabilities

Bedops is a comprehensive toolkit with many utilities. This module currently includes:
- **`gtf_to_bed`**: Convert GTF to BED12 format

### Potential Future Tasks

We welcome contributions to add more bedops functionality:

**Format Conversions:**
- `gff_to_bed`: Convert GFF3 to BED format
- `vcf_to_bed`: Convert VCF to BED format
- `bam_to_bed`: Convert BAM to BED format

**Set Operations:**
- `bedops_intersect`: Find overlapping intervals
- `bedops_merge`: Merge adjacent/overlapping intervals
- `bedops_difference`: Find non-overlapping regions
- `bedops_complement`: Find complementary regions

**Mapping Operations:**
- `bedmap`: Map genomic features to regions
- `closest_features`: Find nearest genomic features

See the [bedops documentation](https://bedops.readthedocs.io/) for complete functionality.

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 1 core
- **Memory**: 2 GB
- **Runtime**: ~1-5 minutes depending on GTF size

### Resource Scaling
For larger annotation files:
- `cpu_cores`: Usually 1 is sufficient (bedops is single-threaded for most operations)
- `memory_gb`: Increase for very large annotation files (4-8 GB for human/mouse genomes)

### Input Requirements
- **GTF file**: Standard GTF format (GTF2.2 or GTF2.5)
  - Works with UCSC, Ensembl, GENCODE, and RefSeq annotations
  - Properly handles transcript_id, gene_id, and exon structures

## Contributing

We encourage contributions to expand this module! To add new bedops tasks:

1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Add your new task to `ww-bedops.wdl` following the existing pattern
3. Update this README with documentation for the new task
4. Add a test case to `testrun.wdl`
5. Test thoroughly with the demonstration workflow
6. Submit a pull request with detailed documentation

**Ideas for contributions:**
- GFF3 to BED conversion
- VCF to BED conversion
- Set operations (intersect, merge, complement)
- Region mapping operations

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[Bedops Documentation](https://bedops.readthedocs.io/)**: Official bedops documentation and tutorials
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
- **[ww-rseqc](../ww-rseqc/)**: RNA-seq QC module that uses bedops for format conversion
