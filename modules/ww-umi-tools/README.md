# ww-umi-tools Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for [UMI-tools](https://github.com/CGATOxford/UMI-tools), a suite for handling Unique Molecular Identifiers (UMIs) in NGS data.

## Overview

Library prep methods like PRO-seq, single-cell RNA-seq, and many capture protocols use random UMIs to distinguish PCR duplicates from independent molecules that happen to map to the same position. After alignment, `umi_tools dedup` collapses reads that share both a mapping position and a UMI, leaving one representative read per unique molecule.

This module wraps `umi_tools dedup` so it can be plugged into WILDS pipelines downstream of any aligner module (`ww-bowtie2`, `ww-bwa`, `ww-star`, etc.) whose upstream FASTQ trimming step encoded UMIs into read names.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-umi-tools.wdl` — task definitions
- **Test workflow**: `testrun.wdl` — zero-config demonstration workflow
- **Documentation**: this README

## Available Tasks

### `dedup`

Deduplicates a coordinate-sorted, indexed BAM using `umi_tools dedup` and re-indexes the result.

**Inputs:**
- `input_bam` (File): Coordinate-sorted BAM with UMIs encoded in read names
- `input_bai` (File): Index for the input BAM
- `sample_name` (String): Sample name used for output file naming
- `paired` (Boolean, default=`true`): Pass `--paired` to dedup (required for paired-end data)
- `umi_separator` (String, default=`":"`): Character separating the read ID from the UMI in the read name (matches `fastp --umi` default)
- `method` (String, default=`"directional"`): Dedup clustering method — one of `unique`, `percentile`, `cluster`, `adjacency`, `directional`
- `extra_args` (String, default=`""`): Additional flags passed verbatim to `umi_tools dedup`
- `cpu_cores` (Int, default=2): CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB

**Outputs:**
- `deduped_bam` (File): Deduplicated BAM file
- `deduped_bai` (File): Index for the deduplicated BAM
- `log` (File): `umi_tools dedup` log file (input/output read counts and per-position stats)

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-umi-tools/ww-umi-tools.wdl" as umi_tools_tasks

workflow my_pipeline {
  input {
    File aligned_bam
    File aligned_bai
    String sample_name
  }

  call umi_tools_tasks.dedup {
    input:
      input_bam = aligned_bam,
      input_bai = aligned_bai,
      sample_name = sample_name,
      paired = true
  }

  output {
    File deduped_bam = dedup.deduped_bam
    File deduped_bai = dedup.deduped_bai
    File dedup_log = dedup.log
  }
}
```

### Advanced Usage Examples

**Single-end data with a non-default UMI separator:**

```wdl
call umi_tools_tasks.dedup {
  input:
    input_bam = aligned_bam,
    input_bai = aligned_bai,
    sample_name = "se_sample",
    paired = false,
    umi_separator = "_"
}
```

**Forward extra arguments (e.g., per-cell dedup for single-cell data):**

```wdl
call umi_tools_tasks.dedup {
  input:
    input_bam = aligned_bam,
    input_bai = aligned_bai,
    sample_name = "sc_sample",
    extra_args = "--per-cell --cell-tag=CB"
}
```

### Integration Examples

This module integrates with other WILDS components:
- **ww-fastp**: encodes UMIs into read names during adapter trimming (use `--umi` family of flags); `ww-umi-tools.dedup` reads them back out
- **ww-bowtie2 / ww-bwa / ww-star**: produce sorted BAMs that feed directly into `dedup`
- **ww-samtools**: useful for downstream filtering or stats on deduped output
- **ww-testdata**: provides `inject_synthetic_umis` for zero-config testing without real UMI data

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that runs end-to-end with no input files:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint umi_tools_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow:
1. Downloads a small paired-end BAM via `ww-testdata.download_bam_data`
2. Injects deterministic synthetic 6-mer UMIs into each read name via `ww-testdata.inject_synthetic_umis` (since the demo BAM has no real UMIs)
3. Runs `umi_tools dedup` on the UMI-tagged BAM and re-indexes the output

## Docker Container

This module uses the `getwilds/umitools:1.1.6` container image, which includes:
- `umi_tools` 1.1.6 (pip install on `python:3.12-bookworm`)
- `samtools` (used to index the deduplicated BAM)

## Citation

> **UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy**
> Smith T, Heger A, Sudbery I (2017)
> *Genome Research*, 27(3): 491–499
> DOI: [10.1101/gr.209601.116](https://doi.org/10.1101/gr.209601.116)

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 2 cores
- **Memory**: 8 GB
- **Runtime**: Scales with input BAM size and read density; small test BAMs complete in under a minute

### Resource Scaling
- `cpu_cores`: Mostly used by the post-dedup `samtools index` step; `umi_tools dedup` itself is single-threaded
- `memory_gb`: Increase for high-coverage BAMs — `umi_tools dedup` holds reads in memory per position bundle

## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for guidelines

## Related Resources

- **[UMI-tools documentation](https://umi-tools.readthedocs.io/)**
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
