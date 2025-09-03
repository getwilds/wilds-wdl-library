# ww-sourmash
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for genomic sequence comparison and metagenome analysis using [sourmash](https://sourmash.readthedocs.io/).

## Overview

This module provides reusable WDL tasks for generating MinHash signatures from genomic sequences, performing similarity searches, and analyzing metagenome composition using **sourmash**. It includes comprehensive workflows for both individual sequence comparison and metagenome analysis.

Designed to be a modular component in the WILDS ecosystem, this workflow is suitable for both standalone execution and integration into larger comparative genomics and metagenomics pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `sourmash_sketch`, `sourmash_search`, `sourmash_gather`, `validate_outputs`
- **Workflow**: `sourmash_example` (demonstration workflow executing all tasks)
- **Container**: `getwilds/sourmash:4.8.2`

## Tasks

### `sourmash_sketch`

Generate MinHash signatures from FASTA/FASTQ sequence files.

**Inputs:**
- `input_file` (File): Input sequence file (FASTA or FASTQ format, optionally compressed)
- `output_name` (String): Name for the output signature file
- `ksize` (String): k-mer size for signature generation (recommended: 21, 31, 51)
- `molecule_type` (String): Molecule type: dna, protein, dayhoff, hp, or nucleotide
- `scaled` (Boolean): Use scaled MinHash (recommended for most applications)
- `scale_factor` (Int): Scale factor for scaled MinHash (lower = more precise, higher = faster)
- `memory_gb` (Int): Memory allocation in GB (default: 4)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `signature_file` (File): MinHash signature file in JSON format

### `sourmash_search`

Compare MinHash signatures to find sequence similarities using containment analysis.

**Inputs:**
- `query_sig` (File): Query MinHash signature file
- `database_sig` (File): Database MinHash signature file to search against
- `ksize` (String): k-mer size to use for search (must match signature k-mer size)
- `threshold` (Float): Minimum similarity threshold for reporting matches (0.0-1.0)
- `memory_gb` (Int): Memory allocation in GB (default: 4)
- `cpu_cores` (Int): Number of CPU cores to use (default: 1)

**Outputs:**
- `results_file` (File): CSV file containing similarity search results

### `sourmash_gather`

Analyze metagenome samples by identifying constituent genomes using sourmash gather.

**Inputs:**
- `query_sig` (File): Query MinHash signature file (typically from a metagenome sample)
- `database_sigs` (Array[File]): Array of database MinHash signature files representing reference genomes
- `ksize` (String): k-mer size to use for gather analysis (must match signature k-mer size)
- `threshold` (Float): Minimum abundance threshold for reporting genome matches (0.0-1.0)
- `memory_gb` (Int): Memory allocation in GB (default: 8)
- `cpu_cores` (Int): Number of CPU cores to use (default: 2)

**Outputs:**
- `results_file` (File): CSV file containing gather analysis results with genome matches and abundances

### `validate_outputs`

Validate sourmash outputs including signatures and analysis results.

**Inputs:**
- `query_signatures` (Array[File]): Query MinHash signature files to validate
- `database_signatures` (Array[File]): Database MinHash signature files to validate
- `search_results` (Array[File]): Search results CSV files to validate
- `gather_results` (Array[File]): Gather results CSV files to validate
- `cpu_cores` (Int): Number of CPU cores to use for validation (default: 1)
- `memory_gb` (Int): Memory allocation in GB for the task (default: 2)

**Outputs:**
- `report` (File): Validation report confirming output integrity and providing basic statistics

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sourmash/ww-sourmash.wdl" as sourmash_tasks

struct GenomeSample {
    String name
    File sequence_file
}

workflow my_comparative_genomics {
  input {
    Array[GenomeSample] query_samples
    Array[File] reference_genomes
    String ksize = "31"
  }

  # Generate signatures for query samples
  scatter (sample in query_samples) {
    call sourmash_tasks.sourmash_sketch as sketch_queries {
      input:
        input_file = sample.sequence_file,
        output_name = sample.name + ".sig",
        ksize = ksize,
        molecule_type = "dna",
        scaled = true,
        scale_factor = 1000
    }
  }

  # Generate signatures for reference genomes
  scatter (ref_genome in reference_genomes) {
    call sourmash_tasks.sourmash_sketch as sketch_refs {
      input:
        input_file = ref_genome,
        output_name = "ref_" + basename(ref_genome) + ".sig",
        ksize = ksize,
        molecule_type = "dna",
        scaled = true,
        scale_factor = 1000
    }
  }

  # Perform similarity searches
  scatter (query_sig in sketch_queries.signature_file) {
    scatter (ref_sig in sketch_refs.signature_file) {
      call sourmash_tasks.sourmash_search {
        input:
          query_sig = query_sig,
          database_sig = ref_sig,
          ksize = ksize,
          threshold = 0.1
      }
    }
  }

  output {
    Array[File] query_signatures = sketch_queries.signature_file
    Array[File] reference_signatures = sketch_refs.signature_file
    Array[Array[File]] similarity_results = sourmash_search.results_file
  }
}
```

### Common Use Cases

**Metagenome analysis:**
```wdl
call sourmash_tasks.sourmash_gather {
  input:
    query_sig = metagenome_signature,
    database_sigs = reference_genome_signatures,
    ksize = "31",
    threshold = 0.05
}
```

**Large-scale genome comparison:**
```wdl
scatter (genome_pair in all_genome_pairs) {
  call sourmash_tasks.sourmash_search {
    input:
      query_sig = genome_pair.query,
      database_sig = genome_pair.reference,
      threshold = 0.1
  }
}
```

**Custom k-mer analysis:**
```wdl
call sourmash_tasks.sourmash_sketch {
  input:
    input_file = sequences,
    ksize = "21",
    molecule_type = "protein",
    scaled = false,
    scale_factor = 500
}
```

### Integration Examples

This module pairs seamlessly with other WILDS modules:
- **ww-sra**: Download sequencing data for comparative analysis (built into demo workflow)
- **ww-testdata**: Automatic test data generation (built into demo workflow)
- **ww-bcftools**: Variant-based comparisons as complement to sequence similarity
- **Custom workflows**: Foundation for any comparative genomics or metagenomics pipeline

## Testing the Module

The module includes a demonstration workflow that can be tested independently:

```bash
# Using Cromwell
java -jar cromwell.jar run ww-sourmash.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-sourmash.wdl -i inputs.json

# Using Sprocket
sprocket run ww-sourmash.wdl inputs.json
```

### Automatic Demo Mode

When no samples or database files are provided, the workflow automatically:
1. Downloads test FASTQ data using `ww-testdata`
2. Downloads reference genome data using `ww-testdata`
3. Generates MinHash signatures for both queries and database
4. Performs similarity searches between all sample-database pairs
5. Conducts metagenome gather analysis
6. Validates all outputs with comprehensive reporting

### Test Input Format

**Minimal input (uses automatic demo data):**
```json
{}
```

**Basic configuration with your own data:**
```json
{
  "sourmash_example.samples": [
    {
      "name": "sample1",
      "fastq_file": "path/to/sample1.fastq.gz"
    },
    {
      "name": "sample2", 
      "fastq_file": "path/to/sample2.fastq.gz"
    }
  ],
  "sourmash_example.database_fastas": [
    "path/to/reference1.fasta",
    "path/to/reference2.fasta"
  ]
}
```

**Advanced configuration:**
```json
{
  "sourmash_example.samples": [
    {
      "name": "metagenome_sample",
      "fastq_file": "path/to/metagenome.fastq.gz"
    }
  ],
  "sourmash_example.database_fastas": [
    "path/to/genome1.fasta",
    "path/to/genome2.fasta",
    "path/to/genome3.fasta"
  ],
  "sourmash_example.ksize": "21",
  "sourmash_example.threshold": 0.1,
  "sourmash_example.molecule_type": "dna",
  "sourmash_example.scaled": true,
  "sourmash_example.scale_factor": 2000,
  "sourmash_example.gather_threshold": 0.01
}
```

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Sufficient storage for signature files and results
- Input sequences in FASTA or FASTQ format (compressed or uncompressed)

## Features

- **Flexible input formats**: Supports FASTA, FASTQ, compressed and uncompressed files
- **Scalable analysis**: Efficient MinHash sketching scales to large datasets
- **Metagenome support**: Specialized gather analysis for microbiome studies
- **Comprehensive validation**: Built-in output validation and reporting
- **Parallel processing**: Multi-sample and multi-database parallelization
- **Memory efficient**: Configurable scaling parameters for different dataset sizes

## Performance Considerations

- **Memory usage**: Scales with scale_factor (lower values = more memory)
- **Compute time**: Depends on sequence file size and number of comparisons
- **Storage requirements**: Signature files are typically much smaller than input sequences
- **Parallelization**: Both scatter-gather and internal parallelization supported

## Output Description

- **Signature files**: MinHash sketches in JSON format for fast similarity comparisons
- **Search results**: CSV files with similarity scores and containment metrics
- **Gather results**: Detailed metagenome composition with abundance estimates
- **Validation report**: Comprehensive validation with signature counts and result statistics

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real genomic datasets for integration testing
- Comprehensive validation of all outputs including signature integrity
- Tests with both single genome and metagenome scenarios

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to sourmash usage, please refer to the [sourmash documentation](https://sourmash.readthedocs.io/) and [GitHub repository](https://github.com/sourmash-bio/sourmash).

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.

