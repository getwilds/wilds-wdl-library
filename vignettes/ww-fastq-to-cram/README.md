# ww-fastq-to-cram
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL vignette for converting paired FASTQ files to unmapped CRAM format with validation.

## Overview

This vignette demonstrates how to convert paired-end FASTQ files from multiple sources (flowcells, lanes, runs, etc.) into unmapped CRAM format using WILDS WDL modules. It's designed for scenarios where sequencing data from the same library was run across multiple sequencing groups and needs to be consolidated into a single CRAM file per sample with proper read group metadata.

The workflow handles:
- Converting FASTQ files to unmapped BAM with proper read group headers
- Merging multiple unmapped BAMs (one per FASTQ group) into a single CRAM per sample
- Validating the final CRAM files for formatting issues

## Workflow Structure

This vignette imports tasks from WILDS WDL modules:
- **ww-gatk**: `fastq_to_sam`, `validate_sam_file`
- **ww-samtools**: `merge_bams_to_cram`

### Input Structure

The workflow uses nested structs to organize input data:

```wdl
struct FastqGroup {
  String group_name
  Array[File] fastq_r1_locations
  Array[File] fastq_r2_locations
}

struct SampleData {
  String sample_name
  String? library_name
  String? sequencing_center
  Array[FastqGroup] fastq_groups
}
```

### Workflow Logic

1. **Per-group conversion**: For each FASTQ group in each sample, convert FASTQ pairs to unmapped BAM with read group: `{sample_name}_{group_name}`
2. **Merge to CRAM**: Merge all group BAMs for each sample into a single CRAM file
3. **Validate**: Run GATK ValidateSamFile on each final CRAM to ensure formatting is correct

## Usage

### Basic Example

Create an `inputs.json` file:

```json
{
  "fastq_to_cram.batch_info": [
    {
      "sample_name": "Sample_A",
      "library_name": "Lib_A",
      "sequencing_center": "SeqCenter",
      "fastq_groups": [
        {
          "group_name": "FC001",
          "fastq_r1_locations": ["/path/to/data/sample_A_FC001_R1.fastq.gz"],
          "fastq_r2_locations": ["/path/to/data/sample_A_FC001_R2.fastq.gz"]
        },
        {
          "group_name": "FC002",
          "fastq_r1_locations": ["/path/to/data/sample_A_FC002_R1.fastq.gz"],
          "fastq_r2_locations": ["/path/to/data/sample_A_FC002_R2.fastq.gz"]
        }
      ]
    }
  ],
  "fastq_to_cram.cpu_cores": 6,
  "fastq_to_cram.memory_gb": 12
}
```

### Running the Workflow

```bash
# Using Cromwell
java -jar cromwell.jar run ww-fastq-to-cram.wdl -i inputs.json

# Using miniWDL
miniwdl run ww-fastq-to-cram.wdl -i inputs.json

# Using Sprocket
sprocket run ww-fastq-to-cram.wdl -i inputs.json
```

## Inputs

### Required

- `batch_info` (Array[SampleData]): Array of samples with their associated FASTQ files and metadata

### Optional

- `cpu_cores` (Int): Number of CPU cores to use for processing (default: 6)
- `memory_gb` (Int): Memory allocation in GB (default: 12)

### SampleData Fields

- `sample_name` (String): Sample name for read group header
- `library_name` (String?): Library name for read group LB tag (optional, defaults to sample_name)
- `sequencing_center` (String?): Sequencing center for read group CN tag (optional, defaults to '.')
- `fastq_groups` (Array[FastqGroup]): Array of FASTQ file groups for this sample

### FastqGroup Fields

- `group_name` (String): Name identifier for this group (e.g., flowcell, lane, or run name)
- `fastq_r1_locations` (Array[File]): R1 FASTQ files for this group
- `fastq_r2_locations` (Array[File]): R2 FASTQ files for this group

## Outputs

- `unmapped_crams` (Array[File]): Unmapped CRAM files (one per sample)
- `unmapped_cram_indexes` (Array[File]): Index files for each unmapped CRAM
- `validation_reports` (Array[File]): Validation reports for each CRAM

## Use Cases

This vignette is ideal for:

1. **Storage optimization**: Reducing storage costs by converting FASTQs to more efficiently compressed unmapped CRAMs while retaining all sequencing information
2. **Multi-source sequencing**: When the same library was sequenced across multiple flowcells, lanes, or runs and needs to be consolidated with proper read group tracking
3. **Metadata consolidation**: Preserving sample, library, and sequencing center information in standardized SAM/BAM/CRAM format headers
4. **Pre-alignment processing**: Creating unmapped CRAMs with proper read group metadata before alignment (required input format for many GATK workflows)
5. **Quality control**: Validating file formats before proceeding with downstream analysis
6. **Pipeline compatibility**: Converting FASTQs to the SAM/BAM/CRAM format expected by alignment and variant calling pipelines

## Read Group Convention

Read groups are created with the convention: `{sample_name}_{group_name}`

Example:
- Sample: "NA12878"
- Group: "HISEQ001" (could be a flowcell, lane, or run identifier)
- Read Group Name: "NA12878_HISEQ001"

This ensures unique read group identifiers even when the same sample is sequenced across multiple groups.

## Performance Considerations

- **CPU allocation**: The workflow benefits from multiple cores for merging operations (recommend 4-8 cores)
- **Memory usage**: GATK FastqToSam typically requires 8GB RAM per task
- **Storage efficiency**: Unmapped CRAMs provide significant storage savings compared to FASTQs due to efficient compression algorithms, while also being smaller than unmapped BAMs. This can result in substantial cost savings for long-term storage of sequencing data, especially for large-scale projects with many samples.
- **Parallelization**: Tasks run in parallel across samples and FASTQ groups, providing good scalability for multi-sample batches

## Technical Details

### GATK FastqToSam Parameters

- Compression level: 5 (balance between speed and size)
- Platform: illumina (default, can be customized)
- Read group attributes: RG, SM, LB, CN populated from input metadata

### Samtools Merge Parameters

- Merging: Uses `samtools merge` with `cpu_cores - 1` threads
- Output format: CRAM with index file
- Threading: Uses `cpu_cores - 1` threads

## Example: Multiple Samples with Multiple FASTQ Groups

```json
{
  "fastq_to_cram.batch_info": [
    {
      "sample_name": "Patient_001",
      "library_name": "WGS_Lib_001",
      "sequencing_center": "WILDS_SeqCore",
      "fastq_groups": [
        {
          "group_name": "HISEQ_FC001",
          "fastq_r1_locations": ["/path/to/data/p001_fc001_L001_R1.fastq.gz", "/path/to/data/p001_fc001_L002_R1.fastq.gz"],
          "fastq_r2_locations": ["/path/to/data/p001_fc001_L001_R2.fastq.gz", "/path/to/data/p001_fc001_L002_R2.fastq.gz"]
        },
        {
          "group_name": "HISEQ_FC002",
          "fastq_r1_locations": ["/path/to/data/p001_fc002_L001_R1.fastq.gz"],
          "fastq_r2_locations": ["/path/to/data/p001_fc002_L001_R2.fastq.gz"]
        }
      ]
    },
    {
      "sample_name": "Patient_002",
      "library_name": "WGS_Lib_002",
      "sequencing_center": "WILDS_SeqCore",
      "fastq_groups": [
        {
          "group_name": "HISEQ_FC001",
          "fastq_r1_locations": ["/path/to/data/p002_fc001_L001_R1.fastq.gz"],
          "fastq_r2_locations": ["/path/to/data/p002_fc001_L001_R2.fastq.gz"]
        }
      ]
    }
  ]
}
```

This example will produce:
- `Patient_001.merged.cram` (combining 2 FASTQ groups)
- `Patient_002.merged.cram` (single FASTQ group)

## Related WILDS Components

- **ww-gatk module**: Provides `fastq_to_sam` and `validate_sam_file` tasks
- **ww-samtools module**: Provides `merge_bams_to_cram` task
- **ww-bwa-gatk vignette**: Demonstrates alignment workflow that could follow this preprocessing

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL vignette, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
