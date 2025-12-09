# ww-varscan
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for somatic variant calling using VarScan2.

## Overview

This module provides a reusable WDL task for somatic variant detection in tumor-normal sample pairs using VarScan2. VarScan2 is a platform-independent variant caller that compares normal and tumor samples at the pileup level to identify somatic mutations. The module uses pre-generated mpileup files from Samtools as input and produces VCF files containing somatic SNV and indel calls.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `somatic`
- **Container**: `getwilds/varscan:2.4.6` (WILDS Docker image with VarScan2 installed)

## Tasks

### `somatic`
Performs somatic variant calling on tumor-normal pairs using pre-generated mpileup files. Uses default parameters except for producing output in VCF format, instead of the VarScan custom format.

**Inputs:**
- `sample_name` (String): Name of the sample (used in output file names)
- `normal_pileup` (File): Samtools mpileup file for the normal sample
- `tumor_pileup` (File): Samtools mpileup file for the tumor sample
- `cpu_cores` (Int): Number of CPU cores (default: 4)
- `memory_gb` (Int): Memory allocation in GB (default: 16)

**Important**: VarScan2 requires **mpileup files** as input, not BAM or FASTQ files. Mpileup files can be generated using SAMtools from aligned BAM files. See the `ww-samtools` [module](https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-samtools) and the [Samtools documentation](http://www.htslib.org/doc/samtools-mpileup.html) for details on creating mpileup files.

**Outputs:**
- `somatic_snvs_vcf` (File): VCF file containing somatic single nucleotide variant calls (`{sample_name}.snp.vcf`)
- `somatic_indels_vcf` (File): VCF file containing somatic insertion/deletion calls (`{sample_name}.indel.vcf`)

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-varscan/ww-varscan.wdl" as varscan_tasks

workflow my_somatic_analysis {
  input {
    String sample_name
    File normal_pileup
    File tumor_pileup
  }

  call varscan_tasks.somatic {
    input:
      sample_name = sample_name,
      normal_pileup = normal_pileup,
      tumor_pileup = tumor_pileup
  }

  output {
    File snvs = somatic.somatic_snvs_vcf
    File indels = somatic.somatic_indels_vcf
  }
}
```

### Advanced Usage Examples

**Custom resource allocation:**
```wdl
call varscan_tasks.somatic {
  input:
    sample_name = "my_sample",
    normal_pileup = normal_pileup,
    tumor_pileup = tumor_pileup,
    cpu_cores = 8,
    memory_gb = 32
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Download test BAM files to generate mpileup files
- **ww-samtools**: Generate mpileup files from BAM alignments before variant calling

**Complete pipeline with mpileup generation:**
```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-samtools/ww-samtools.wdl" as samtools_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-varscan/ww-varscan.wdl" as varscan_tasks

workflow somatic_variant_calling_pipeline {
  input {
    String sample_name
    File normal_bam
    File tumor_bam
    File ref_fasta
  }

  # Generate mpileup for normal sample
  call samtools_tasks.mpileup as normal_mpileup {
    input:
      bamfile = normal_bam,
      ref_fasta = ref_fasta,
      sample_name = sample_name + "_normal"
  }

  # Generate mpileup for tumor sample
  call samtools_tasks.mpileup as tumor_mpileup {
    input:
      bamfile = tumor_bam,
      ref_fasta = ref_fasta,
      sample_name = sample_name + "_tumor"
  }

  # Call somatic variants
  call varscan_tasks.somatic {
    input:
      sample_name = sample_name,
      normal_pileup = normal_mpileup.pileup,
      tumor_pileup = tumor_mpileup.pileup
  }

  output {
    File somatic_snvs = somatic.somatic_snvs_vcf
    File somatic_indels = somatic.somatic_indels_vcf
  }
}
```

**Multi-sample batch processing:**
```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-samtools/ww-samtools.wdl" as samtools_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-varscan/ww-varscan.wdl" as varscan_tasks

struct TumorNormalPair {
  String sample_name
  File normal_bam
  File tumor_bam
}

workflow batch_somatic_variant_calling {
  input {
    Array[TumorNormalPair] samples
    File ref_fasta
  }

  scatter (sample in samples) {
    # Generate mpileups
    call samtools_tasks.mpileup as normal_mpileup {
      input:
        bamfile = sample.normal_bam,
        ref_fasta = ref_fasta,
        sample_name = sample.sample_name + "_normal"
    }

    call samtools_tasks.mpileup as tumor_mpileup {
      input:
        bamfile = sample.tumor_bam,
        ref_fasta = ref_fasta,
        sample_name = sample.sample_name + "_tumor"
    }

    # Call variants
    call varscan_tasks.somatic {
      input:
        sample_name = sample.sample_name,
        normal_pileup = normal_mpileup.pileup,
        tumor_pileup = tumor_mpileup.pileup
    }
  }

  output {
    Array[File] all_snvs = somatic.somatic_snvs_vcf
    Array[File] all_indels = somatic.somatic_indels_vcf
  }
}
```

## Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 16-32 GB recommended (scales with dataset size and genome coverage)
- **CPUs**: 4-8 cores recommended
- **Runtime**: Depends on pileup file size
- **Storage requirements**: Ensure adequate disk space for input pileup files and output VCF files

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support for containerized execution
- Pre-generated mpileup files for both tumor and normal samples (can be generated using Samtools)
- Java runtime (included in Docker container)

## Module Development

This module follows WILDS WDL Library standards. For development and testing:

1. Ensure the Docker image `getwilds/varscan:2.4.6` is available
2. Test with small pileup files first
3. Inspect outputs

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to VarScan usage or configuration, please refer to the [VarScan documentation](https://dkoboldt.github.io/varscan/). Please make sure to cite their work if you use VarScan in your analyses:

Koboldt, D., Zhang, Q., Larson, D. et al. VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Res 22, 568–576 (2012). https://doi.org/10.1101/gr.129684.111

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
