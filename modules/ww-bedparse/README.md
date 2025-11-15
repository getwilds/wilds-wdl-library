# ww-bedparse

[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for genomic file format conversions using [bedparse](https://github.com/tleonardi/bedparse).

## Overview

This module provides a simple interface for converting GTF (Gene Transfer Format) files to BED12 format using the bedparse Python package. The BED12 format is commonly required by RNA-seq quality control tools like RSeQC.

## Key Features

- **Simple GTF to BED12 conversion**: One-command conversion that properly handles exon structures and CDS boundaries
- **Flexible filtering**: Optional filtering by GTF attributes
- **Extra fields**: Can include additional GTF fields in output

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library).

### Available Tasks

- **gtf2bed**: Converts GTF annotation files to BED12 format

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support

### Basic Usage

Import the module and call the `gtf2bed` task:

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bedparse/ww-bedparse.wdl" as ww_bedparse

workflow convert_gtf {
  call ww_bedparse.gtf2bed {
    input:
      gtf_file = "path/to/annotation.gtf"
  }

  output {
    File bed12 = gtf2bed.bed_file
  }
}
```

### Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `gtf_file` | GTF annotation file to convert | File | Yes | - |
| `extra_fields` | Comma-separated list of extra GTF fields to add after column 12 (e.g., "gene_id,gene_name") | String | No | - |
| `filter_key` | GTF extra field on which to apply filtering | String | No | - |
| `filter_type` | Comma-separated list of filterKey field values to retain | String | No | - |
| `cpu_cores` | Number of CPU cores | Int | No | 1 |
| `memory_gb` | Memory allocation in GB | Int | No | 2 |

### Advanced Usage with Filtering

```wdl
call ww_bedparse.gtf2bed {
  input:
    gtf_file = "annotation.gtf",
    extra_fields = "gene_id,gene_name",
    filter_key = "gene_type",
    filter_type = "protein_coding"
}
```

## Output Files

| Output | Description |
|--------|-------------|
| `bed_file` | BED12-formatted annotation file with proper exon blocks and CDS coordinates |

## BED12 Format

The output BED12 file contains 12 columns:

1. **chrom**: Chromosome name
2. **chromStart**: Start position (0-based)
3. **chromEnd**: End position
4. **name**: Gene/transcript name
5. **score**: Score (usually 0)
6. **strand**: Strand (+/-)
7. **thickStart**: CDS start position
8. **thickEnd**: CDS end position
9. **itemRgb**: RGB color (usually 0)
10. **blockCount**: Number of exons
11. **blockSizes**: Comma-separated exon sizes
12. **blockStarts**: Comma-separated exon starts (relative to chromStart)

## GTF Requirements

The GTF file must:
- Follow Ensembl GTF format
- Contain 'transcript' features in field 3
- Contain 'exon' features in field 3
- Optionally contain 'CDS', 'start_codon', or 'stop_codon' features for proper thickStart/thickEnd annotation

## Testing the Module

The module includes a test workflow that automatically runs with test data:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint bedparse_example
```

The test workflow:
1. Downloads a small GTF file (chromosome 1 subset)
2. Converts it to BED12 format
3. Outputs both the original GTF and converted BED file

## Use Cases

This module is commonly used with:
- **RSeQC**: RNA-seq quality control tools that require BED12 format
- **UCSC Genome Browser**: Visualization of gene structures
- **Custom analysis**: Any workflow requiring BED12 gene annotations

## Related WILDS Components

- **ww-rseqc module**: Uses BED12 files for RNA-seq QC
- **ww-star-deseq2 vignette**: Complete RNA-seq pipeline using this module for QC
- **ww-testdata module**: Provides test GTF files

## About bedparse

[bedparse](https://github.com/tleonardi/bedparse) is a Python module and CLI tool created by Tommaso Leonardi for performing operations on BED files. It provides reliable GTF to BED12 conversion that properly handles:
- Exon grouping by transcript
- CDS boundary annotation
- Coordinate system conversion (GTF 1-based to BED 0-based)

If you use this module in your research, please cite bedparse:

> **bedparse: feature extraction from BED files**
> Tommaso Leonardi
> Journal of Open Source Software, 2019, 4(34), 1228
> DOI: [10.21105/joss.01228](https://doi.org/10.21105/joss.01228)

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
