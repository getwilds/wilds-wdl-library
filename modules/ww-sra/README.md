# ww-sra
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for downloading genomic data using the SRA toolkit.

## Overview

This workflow enables automated downloading of sequencing data from the NCBI Sequence Read Archive (SRA) using the SRA toolkit. It handles both single-end and paired-end reads, automatically detecting the read type and processing accordingly.

The workflow uses `parallel-fastq-dump` for efficient, multi-threaded downloading of FASTQ files from SRA accessions.

## Features

- Automated downloading of FASTQ files from SRA accessions
- Automatic detection of paired-end vs. single-end reads
- Multi-threaded downloading for improved performance
- Standardized output naming for downstream processing
- Compatible with the WILDS workflow ecosystem

## Usage

### Requirements

- [Cromwell](https://cromwell.readthedocs.io/), [MiniWDL](https://github.com/chanzuckerberg/miniwdl), [Sprocket](https://sprocket.bio/), or another WDL-compatible workflow executor
- Docker/Apptainer (the workflow uses `getwilds/sra-tools:3.1.1` container)

### Basic Usage

1. Create an inputs JSON file with your SRA accessions:

```json
{
  "SRA_Download.sra_id_list": ["SRR32657057", "SRR33450906"]
}
```

2. Run the workflow using Cromwell:

```bash
# Cromwell
java -jar cromwell.jar run ww-sra.wdl --inputs ww-sra-inputs.json

# miniWDL
miniwdl run ww-sra.wdl -i ww-sra-inputs.json

# Sprocket
sprocket run ww-sra.wdl ww-sra-inputs.json
```

### Detailed Options

The workflow accepts the following inputs:

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `sra_id_list` | List of SRA accessions to download | Array[String] | Yes | - |

### Output Files

The workflow produces the following outputs:

| Output | Description | Type |
|--------|-------------|------|
| `output_R1` | R1 FASTQ files for each sample | Array[File] |
| `output_R2` | R2 FASTQ files for each sample (or empty files for single-end) | Array[File] |
| `is_paired_end` | Boolean indicators for paired-end status | Array[Boolean] |

## For Fred Hutch Users

For Fred Hutch users, we recommend using [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this workflow directly to the on-premise HPC cluster. To do this:

1. Clone or download this repository
2. Update `ww-sra-inputs.json` with your desired SRA accessions
3. Update `ww-sra-options.json` with your preferred output location (`final_workflow_outputs_dir`)
4. Submit the WDL file along with your custom JSONs to the Fred Hutch cluster via PROOF

### Example Options File

```json
{
    "workflow_failure_mode": "ContinueWhilePossible",
    "write_to_cache": true,
    "read_from_cache": true,
    "default_runtime_attributes": {
        "maxRetries": 1
    },
    "final_workflow_outputs_dir": "/your/output/path/",
    "use_relative_output_paths": true
}
```

## Advanced Usage

### Integrating with Other Workflows

This workflow is designed to be modular and can be easily integrated with other WILDS WDL workflows, such as alignment or analysis pipelines. The output files are named in a standardized format (`SRR_ID_1.fastq.gz` and `SRR_ID_2.fastq.gz`) for easy downstream processing.

### Future Improvements

Future versions of this workflow will include:
- Metadata retrieval using `efetch`
- Checksum validation for downloaded files
- Enhanced logging for better visibility

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on our [issue tracker](https://github.com/getwilds/ww-sra/issues).

## Contributing

If you would like to contribute to this WILDS WDL workflow, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.

