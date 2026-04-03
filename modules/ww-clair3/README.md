# ww-clair3 Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for [Clair3](https://github.com/HKU-BAL/Clair3), a deep learning-based germline small variant caller. Clair3 uses a two-stage pileup and full-alignment approach for high-accuracy SNP and indel calling from Oxford Nanopore (ONT), PacBio HiFi, and Illumina sequencing data.

## Overview

This module wraps Clair3's variant calling functionality for use in WILDS WDL workflows. Clair3 combines a fast pileup model for initial candidate identification with a more thorough full-alignment model for difficult regions, achieving high accuracy across different sequencing platforms and coverage depths.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-clair3.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `run_clair3`

Runs Clair3 to call germline variants from aligned reads using deep learning pileup and full-alignment models.

**Inputs:**
- `sample_name` (String): Name identifier for the sample
- `input_bam` (File): Aligned reads in BAM format (must be sorted and indexed)
- `input_bam_index` (File): Index file for the input BAM
- `ref_fasta` (File): Reference genome FASTA file (must be indexed with .fai)
- `ref_fasta_index` (File): Index file for the reference FASTA
- `platform` (String, default="ont"): Sequencing platform: `ont`, `hifi`, or `ilmn`
- `model_path` (String, default="/opt/models/ont"): Path to Clair3 model directory
- `gvcf_enabled` (Boolean, default=false): Whether to produce a gVCF file for joint genotyping
- `bed_file` (File?, optional): BED file to restrict variant calling to specific regions
- `ctg_name` (String?, optional): Contig/chromosome name to restrict calling
- `include_all_ctgs` (Boolean, default=false): Call on all contigs (required for non-human species)
- `pileup_only` (Boolean, default=false): Use only the pileup model for faster variant calling (skips full-alignment stage)
- `gpu_enabled` (Boolean, default=false): Enable GPU acceleration for Clair3 inference
- `cpu_cores` (Int, default=8): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=32): Memory allocated for the task in GB

**Outputs:**
- `output_vcf` (File): VCF file containing merged variant calls
- `output_vcf_index` (File): Index file for the output VCF
- `output_gvcf` (Array[File]): gVCF file for joint genotyping (empty when gvcf_enabled is false)
- `output_gvcf_index` (Array[File]): gVCF index file (empty when gvcf_enabled is false)

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-clair3/ww-clair3.wdl" as clair3_tasks

workflow my_variant_calling {
  input {
    String sample_name
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
  }

  call clair3_tasks.run_clair3 { input:
    sample_name = sample_name,
    input_bam = input_bam,
    input_bam_index = input_bam_index,
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    platform = "ont"
  }

  output {
    File variants = run_clair3.output_vcf
    File variants_index = run_clair3.output_vcf_index
  }
}
```

### Advanced Usage Examples

**PacBio HiFi variant calling with gVCF output:**
```wdl
call clair3_tasks.run_clair3 { input:
  sample_name = "my_sample",
  input_bam = hifi_bam,
  input_bam_index = hifi_bam_index,
  ref_fasta = reference_fasta,
  ref_fasta_index = reference_fasta_index,
  platform = "hifi",
  model_path = "/opt/models/hifi",
  gvcf_enabled = true,
  cpu_cores = 16,
  memory_gb = 64
}
```

**GPU-accelerated variant calling:**
```wdl
call clair3_tasks.run_clair3 { input:
  sample_name = "gpu_sample",
  input_bam = ont_bam,
  input_bam_index = ont_bam_index,
  ref_fasta = reference_fasta,
  ref_fasta_index = reference_fasta_index,
  platform = "ont",
  gpu_enabled = true,
  cpu_cores = 8,
  memory_gb = 32
}
```

**Targeted variant calling with BED regions:**
```wdl
call clair3_tasks.run_clair3 { input:
  sample_name = "targeted_sample",
  input_bam = targeted_bam,
  input_bam_index = targeted_bam_index,
  ref_fasta = reference_fasta,
  ref_fasta_index = reference_fasta_index,
  platform = "ilmn",
  model_path = "/opt/models/ilmn",
  bed_file = target_regions_bed
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of test BAM and reference data
- **ww-gatk**: Can be used alongside GATK for comparison variant calling pipelines
- **ww-deepvariant**: Alternative variant caller for benchmarking

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint clair3_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads reference genome (chr1) and BAM test data using `ww-testdata`
2. Runs Clair3 variant calling on the test sample using the Illumina model (pileup-only mode for speed)
3. Outputs a VCF file with germline variant calls

## Docker Container

This module uses the `getwilds/clair3:2.0.0` container image, which includes:
- Clair3 v2.0.0 (PyTorch-based)
- Pre-trained models for ONT, PacBio HiFi, and Illumina at `/opt/models/`
- Samtools and other required dependencies

## Citation

> Clair3: Symphonizing pileup and full-alignment for deep learning-based long-read variant calling
> Zheng Z, Li S, Su J, Leung AW, Lam TW, Luo R.
> Nature Computational Science (2022)
> DOI: https://doi.org/10.1038/s43588-022-00387-x

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 8 cores
- **Memory**: 32 GB
- **Runtime**: Varies by data size and platform; typically 30-60 minutes for 30x WGS

### Resource Scaling
- `cpu_cores`: Clair3 parallelizes by splitting chromosomes into chunks (4 threads each)
- `memory_gb`: Increase for high-coverage samples or when calling on many contigs
- For targeted sequencing with BED files, lower resources may be sufficient

## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for detailed guidelines

## Related Resources

- **[Clair3 GitHub](https://github.com/HKU-BAL/Clair3)**: Official Clair3 repository with documentation
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
