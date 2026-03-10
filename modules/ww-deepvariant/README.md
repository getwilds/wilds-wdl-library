# ww-deepvariant Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for [DeepVariant](https://github.com/google/deepvariant), Google's deep-learning-based germline variant caller. DeepVariant uses a convolutional neural network to call variants from aligned sequencing reads, supporting WGS, WES, PacBio HiFi, and Oxford Nanopore technologies.

## Overview

This module wraps DeepVariant's `run_deepvariant` pipeline, which internally performs three steps:
1. **make_examples** - converts read pileups into image tensors
2. **call_variants** - classifies each tensor using a CNN model
3. **postprocess_variants** - converts model output to VCF/gVCF format

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-deepvariant.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `run_deepvariant`

Runs DeepVariant to call germline variants from aligned reads. Optionally produces a gVCF file for downstream joint genotyping.

**Inputs:**
- `sample_name` (String): Name identifier for the sample
- `input_bam` (File): Aligned reads in BAM format
- `input_bam_index` (File): Index file for the input BAM
- `ref_fasta` (File): Reference genome FASTA file
- `ref_fasta_index` (File): Index file for the reference FASTA
- `model_type` (String, default="WGS"): Sequencing model type (WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA, MASSEQ)
- `output_gvcf_enabled` (Boolean, default=false): Whether to also produce a gVCF file for downstream joint genotyping
- `regions` (String?, optional): Genomic regions to restrict variant calling
- `cpu_cores` (Int, default=8): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=32): Memory allocated for the task in GB

**Outputs:**
- `output_vcf` (File): VCF file containing variant calls
- `output_vcf_index` (File): Index file for the output VCF
- `output_gvcf` (Array[File]): gVCF file for downstream joint genotyping (empty when `output_gvcf_enabled` is false)
- `output_gvcf_index` (Array[File]): Index file for the gVCF (empty when `output_gvcf_enabled` is false)

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deepvariant/ww-deepvariant.wdl" as deepvariant_tasks

struct DeepVariantSample {
    String name
    File bam
    File bai
}

workflow my_variant_calling {
  input {
    Array[DeepVariantSample] samples
    File ref_fasta
    File ref_fasta_index
  }

  scatter (sample in samples) {
    call deepvariant_tasks.run_deepvariant {
      input:
        sample_name = sample.name,
        input_bam = sample.bam,
        input_bam_index = sample.bai,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }
  }

  output {
    Array[File] vcfs = run_deepvariant.output_vcf
  }
}
```

### Advanced Usage Examples

**WES mode with custom resources:**
```wdl
call deepvariant_tasks.run_deepvariant {
  input:
    sample_name = "exome_sample",
    input_bam = exome_bam,
    input_bam_index = exome_bai,
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    model_type = "WES",
    cpu_cores = 16,
    memory_gb = 64
}
```

**PacBio HiFi mode with gVCF output:**
```wdl
call deepvariant_tasks.run_deepvariant {
  input:
    sample_name = "pacbio_sample",
    input_bam = hifi_bam,
    input_bam_index = hifi_bai,
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    model_type = "PACBIO",
    output_gvcf_enabled = true
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of test BAM and reference data for demonstrations
- **ww-bwa / ww-bowtie2**: Upstream alignment modules to produce BAM inputs
- **ww-gatk**: Can be used alongside GATK for comparison or joint genotyping pipelines

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint deepvariant_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads a reference genome (chr1) using `ww-testdata`
2. Downloads test BAM files using `ww-testdata`
3. Calls variants on each sample using DeepVariant (VCF only)
4. Runs an additional call with `output_gvcf_enabled = true` to test gVCF output

## Docker Container

This module uses the `google/deepvariant:1.10.0` container image, which includes:
- DeepVariant variant caller (v1.10.0)
- Pre-trained CNN models for WGS, WES, PacBio, and ONT
- All necessary dependencies for variant calling

## Citation

> Poplin, R., Chang, P.C., Alexander, D. et al.
> A universal SNP and small-indel variant caller using deep neural networks.
> Nat Biotechnol 36, 983-987 (2018).
> DOI: https://doi.org/10.1038/nbt.4235

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 8 cores
- **Memory**: 32 GB
- **Runtime**: Varies by genome size and coverage depth

### Resource Scaling
- `cpu_cores`: DeepVariant benefits significantly from multi-core processing via `--num_shards`
- `memory_gb`: WGS typically needs 32+ GB; WES can use less (16 GB)
- For production WGS at 30x coverage, consider 16+ cores and 64 GB memory

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
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for detailed guidelines

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
- **[DeepVariant GitHub](https://github.com/google/deepvariant)**: Official DeepVariant repository
