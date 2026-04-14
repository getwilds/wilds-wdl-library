# ww-gffread Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module wrapping the [`gffread`](https://github.com/gpertea/gffread) tool for GTF/GFF3 annotation manipulation. The primary purpose is normalizing bacterial NCBI GTFs so they work in eukaryotic-style RNA-seq pipelines.

## Overview

`gffread` is a small, fast utility from the Cufflinks/StringTie ecosystem for parsing, filtering, converting, and normalizing GTF and GFF3 annotation files. This module exposes the subset of `gffread` functionality most relevant to WILDS pipelines.

### Why this module exists: the bacterial GTF problem

NCBI RefSeq bacterial GTFs are structurally different from eukaryotic GTFs in a way that silently breaks standard RNA-seq pipelines. A typical NCBI bacterial GTF for *Pseudomonas aeruginosa* PAO1 contains:

- ~5500 `CDS` rows (one per protein-coding gene)
- ~5700 `gene` rows
- ~100 `exon` rows (only for tRNAs/rRNAs)

Tools like STAR (in GeneCounts mode) and RSeQC only consider `exon` features when building gene-level count tables and annotation summaries. When fed a bacterial GTF, these tools silently report per-gene counts for only the ~100 non-coding RNAs — dropping 98% of protein-coding genes from the results. Downstream DESeq2 then produces differential expression output for only ~100 genes instead of the ~5500 you expect.

The `normalize_gtf` task in this module fixes this by synthesizing `exon` features from `CDS` records using `gffread --force-exons`. The normalized GTF then works correctly with any downstream tool that expects eukaryotic-style exon annotations, with no user-visible changes or configuration required.

The task is safe to apply to any GTF: well-formed eukaryotic GTFs that already have proper `exon` annotations pass through effectively unchanged.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-gffread.wdl` - Contains task definitions
- **Test workflow**: `testrun.wdl` - Zero-config demonstration workflow
- **Documentation**: This README

## Available Tasks

### `normalize_gtf`

Normalizes a GTF file so downstream tools see exon features for every transcript. Synthesizes exon features from CDS records for GTFs that lack them (typical of NCBI bacterial GTFs). Eukaryotic GTFs with proper exon annotations pass through effectively unchanged.

**Inputs:**
- `input_gtf` (File): Input GTF file to normalize. Can be any GTF — bacterial or eukaryotic.
- `output_prefix` (String, default=`"normalized"`): Prefix used for output filenames
- `cpu_cores` (Int, default=1): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=2): Memory allocated for the task in GB

**Outputs:**
- `normalized_gtf` (File): GTF file with exon features synthesized from CDS records where needed

**A note on gffread + NCBI bacterial GTFs:** `gffread` 0.12.7 cannot parse NCBI bacterial GTFs directly — it reports `Error: no valid ID found for GFF record` because NCBI GTFs have `transcript_id ""` (empty string) on `gene`-type rows. The `normalize_gtf` task handles this automatically by stripping `gene`-type rows before invoking `gffread`. These rows are redundant with the corresponding `CDS` rows for coordinate and identifier information, so nothing is lost.

### `gff3_to_gtf`

Converts a GFF3 annotation file to GTF format. Useful when an upstream source (e.g. some Ensembl Bacteria releases) only publishes GFF3 but downstream tools expect GTF.

**Inputs:**
- `input_gff3` (File): Input GFF3 file to convert
- `output_prefix` (String, default=`"converted"`): Prefix used for output filenames
- `cpu_cores` (Int, default=1): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=2): Memory allocated for the task in GB

**Outputs:**
- `gtf_file` (File): GTF-format annotation converted from the input GFF3

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gffread/ww-gffread.wdl" as gffread_tasks

workflow my_rnaseq_pipeline {
  input {
    File reference_gtf
    # ... other inputs
  }

  # Ensure the GTF has proper exon features before handing it to STAR / RSeQC / etc.
  call gffread_tasks.normalize_gtf { input:
      input_gtf = reference_gtf
  }

  # Downstream tasks then use normalize_gtf.normalized_gtf instead of reference_gtf
  # ...
}
```

### Integration Example: fixing bacterial RNA-seq in a STAR-based pipeline

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gffread/ww-gffread.wdl" as gffread_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-star/ww-star.wdl" as star_tasks

workflow bacterial_rnaseq {
  input {
    File reference_fasta
    File reference_gtf  # e.g. NCBI PAO1 GTF
    # ...
  }

  call gffread_tasks.normalize_gtf { input:
      input_gtf = reference_gtf
  }

  call star_tasks.build_index { input:
      reference_fasta = reference_fasta,
      reference_gtf = normalize_gtf.normalized_gtf
  }

  # ... alignment, counting, DE analysis
}
```

## Testing the Module

### Automatic Demo Mode

The `testrun.wdl` workflow runs `normalize_gtf` against two contrasting inputs with no user configuration:

1. **Bacterial case**: the NCBI PAO1 GTF (via `ww-testdata.download_pao1_ref`). Exercises the CDS→exon synthesis path.
2. **Eukaryotic case**: the Ensembl human chromosome 15 GTF (via `ww-testdata.download_jcast_test_data`). Exercises the pass-through path.

Run locally with:

```bash
make run_sprocket NAME=ww-gffread
# or
make run_miniwdl NAME=ww-gffread
```

After the run completes, inspect the normalized GTF output files to verify:
- The bacterial case has ~5500 `exon` rows (up from ~100 in the input).
- The eukaryotic case has roughly the same number of `exon` rows as the input.

## Docker Container

This module uses the **`getwilds/gffread:0.12.7`** container image from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library), which includes gffread version 0.12.7.

## Citation

If you use `gffread` in your work, please cite:

> Pertea G and Pertea M. GFF Utilities: GffRead and GffCompare. *F1000Research* 2020, 9:304 (https://doi.org/10.12688/f1000research.23297.2)

## Parameters and Resource Requirements

### Default Resources

Both tasks default to 1 CPU and 2 GB memory, which is ample for typical annotation files (bacterial GTFs are ~5 MB, a human GTF is ~1 GB). Increase `memory_gb` only if processing extraordinarily large annotations.

## Support and Feedback

Report issues or request features at the [wilds-wdl-library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Related Resources

- [gffread GitHub](https://github.com/gpertea/gffread)
- [gffread documentation](https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)
- [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library)
