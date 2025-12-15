# ww-sourmash
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for genomic sequence comparison and metagenome analysis using [sourmash](https://sourmash.readthedocs.io/).

## Overview

This module provides reusable WDL tasks for generating sourmash signatures from genomic sequences, analyzing metagenome composition and performing similarity searches.

Designed to be a modular component in the WILDS ecosystem, this module is suitable for both standalone execution and integration into larger comparative genomics and metagenomics pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `sketch`, `gather`, `compare`
- **Test workflow**: `testrun.wdl` (demonstration workflow with automatic test data support)
- **Container**: `getwilds/sourmash:4.8.2`

## Tasks

### `sketch`

Generate sourmash MinHash signatures from BAM or FASTA files.

**Inputs:**
- `infile` (File): Input BAM or FASTA file
- `bam_as_input` (Boolean): Whether the input file is a BAM
- `k_value` (Int): k-mer size for sourmash sketch (default: 31)
- `scaled` (Int): Scaled value for sourmash sketch (default: 1000)
- `track_abundance` (Boolean): Whether to track k-mer abundance (default: true)
- `cpu_cores` (Int): Number of CPU cores to use (only used for BAM input) (default: 4)
- `memory_gb` (Int): Memory allocated for the task in GB (default: 8)
- `output_name` (String?): Optional custom output name (defaults to input file basename)

**Important:** This WDL assumes DNA input (not RNA or protein)

**Outputs:**
- `sig` (File): Sourmash sketch file`

### `gather`

Run sourmash gather to decompose metagenome samples and identify constituent genomes.

**Inputs:**
- `query_sig` (File): Query sourmash sketch file (typically from a metagenome sample)
- `reference_databases_sigs` (Array[File]): Array of reference database(s) and signature(s) to search against
- `threshold_bp` (Int): Minimum number of base pairs to report a match (default: 50000)
- `memory_gb` (Int): Memory allocated for the task in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 4)
- `output_name` (String?): Optional custom output name (defaults to query basename)

**Important:** Signatures and databases must use identical parameters (same k-mer size, same scaled factor)

**Outputs:**
- `result` (File): CSV file of sourmash gather results

### `compare`

Generate similarity matrix from signature files using sourmash compare.

**Inputs:**
- `sig_inputs` (Array[File]): Array of input signature files to compare
- `save_name` (String): Name to use for output files
- `k_value` (Int): Value of k used for sourmash sketch files
- `memory_gb` (Int): Memory allocated for the task in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 4)

**Outputs:**
- `npy` (File): Numpy binary matrix file of angular similarity matrix
- `csv` (File): CSV file of angular similarity matrix

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sourmash/ww-sourmash.wdl" as sourmash_tasks

struct MetagenomeSample {
    String name
    File bam_file
}

workflow my_metagenome_analysis {
  input {
    Array[MetagenomeSample] samples
    Array[File] reference_databases
    Int k_value = 31
    Int scaled = 1000
  }

  # Generate signatures for samples
  scatter (sample in samples) {
    call sourmash_tasks.sketch {
      input:
        infile = sample.bam_file,
        bam_as_input = true,
        k_value = k_value,
        scaled = scaled,
        track_abundance = true,
        output_name = sample.name
    }
  }

  # Run gather to identify genomes
  scatter (sig in sketch.sig) {
    call sourmash_tasks.gather {
      input:
        query_sig = sig,
        reference_databases_sigs = reference_databases
    }
  }

  # Compare samples to each other
  call sourmash_tasks.compare {
    input:
      sig_inputs = sketch.sig,
      save_name = "sample_comparison",
      k_value = k_value
  }

  output {
    Array[File] signatures = sketch.sig
    Array[File] gather_results = gather.result
    File similarity_matrix_npy = compare.npy
    File similarity_matrix_csv = compare.csv
  }
}
```

### Common Use Cases

**Metagenome analysis with bacterial reference database:**
```wdl
# Generate signatures from unmapped reads
call sourmash_tasks.sketch {
  input:
    infile = unmapped_reads_bam,
    bam_as_input = true,
    k_value = 31,
    scaled = 1000,
    track_abundance = true
}

# Identify bacterial genomes using GTDB database
call sourmash_tasks.gather {
  input:
    query_sig = sketch.sig,
    reference_databases_sigs = ["gtdb-reps-rs226-k31.dna.zip"]
}
```

**Reference genome signature generation:**
```wdl
call sourmash_tasks.sketch  {
  input:
    infile = ref_genome_fasta,
    bam_as_input = false,
    k_value = 31,
    scaled = 1000,
    track_abundance = false  # Not needed for references
}
```

### Integration Examples

This module pairs well with other WILDS modules:

- **ww-testdata**: Download human reference genomes to check for matches to human contamination

## Resource allocation

- **Memory**: Scales with scaled parameter (lower = more memory) and number of signatures
- **CPUs**: BAM sketching benefits from multiple cores
- **Storage requirements**: Signature files are typically much smaller than input sequences (scaled=1000 gives ~1000x compression)


## Docker Image

This module uses `getwilds/sourmash:4.8.2` from DockerHub.

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support for containerized execution
- Sufficient computational resources

## Features

- **Flexible input formats**: Supports both FASTA and BAM files
- **Parallel processing**: Can perform multi-sample parallelization via scatter-gather logic

## Module Development

This module is part of the WILDS WDL Library and follows WILDS development practices.

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to sourmash usage, please refer to the [sourmash documentation](https://sourmash.readthedocs.io/) and [GitHub repository](https://github.com/sourmash-bio/sourmash).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
