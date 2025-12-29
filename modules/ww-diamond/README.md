# ww-diamond
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for fast protein sequence alignment using [DIAMOND](https://github.com/bbuchfink/diamond).

## Overview

This module provides reusable WDL tasks for performing rapid protein sequence alignment with DIAMOND. DIAMOND is a high-performance alternative to BLAST that is up to 100 - 10,000 times faster while maintaining similar sensitivity. The module supports database creation from protein FASTA files and BLASTP alignment for protein-protein searches.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `make_database`, `diamond_blastp`
- **Test workflow**: `testrun.wdl` (demonstration workflow with automatic test data support)
- **Container**: `staphb/diamond:2.1.16` (StaPH-B Docker image with DIAMOND installed)

## Tasks

### `make_database`

Creates a DIAMOND database from a protein FASTA file.

**Inputs:**
- `fasta` (File): Input FASTA file containing protein sequences
- `memory_gb` (Int): Memory allocation in GB (default: 2)
- `cpu_cores` (Int): Number of CPU cores (default: 1)

**Outputs:**
- `diamond_db` (File): DIAMOND database file (.dmnd) created from input FASTA

### `diamond_blastp`

Aligns protein sequences to a DIAMOND database using the BLASTP algorithm.

**Inputs:**
- `diamond_db` (File): DIAMOND database file (.dmnd) from `make_database`
- `query` (File): Query FASTA file containing protein sequences to align
- `align_id` (String): Minimum identity percentage for alignment (default: "50")
- `top_pct` (String): Top percentage of hits to keep; 0 = keep all (default: "0")
- `query_cover` (String): Minimum query coverage percentage (default: "50")
- `subject_cover` (String): Minimum subject coverage percentage (default: "50")
- `blocksize` (Float): Billions of sequence letters to process at a time (default: "2.0")
- `memory_gb` (Int): Memory allocation in GB (default: 2)
- `cpu_cores` (Int): Number of CPU cores (default: 1)

**Outputs:**
- `aln` (File): Compressed alignment file in tabular format (outfmt 6), gzip-compressed

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-diamond/ww-diamond.wdl" as diamond_tasks

workflow my_protein_alignment {
  input {
    File reference_proteome
    File query_proteins
  }

  # Build DIAMOND database
  call diamond_tasks.make_database {
    input:
      fasta = reference_proteome,
      memory_gb = 4,
      cpu_cores = 2
  }

  # Align query proteins
  call diamond_tasks.diamond_blastp {
    input:
      diamond_db = make_database.diamond_db,
      query = query_proteins,
      memory_gb = 4,
      cpu_cores = 4
  }

  output {
    File database = make_database.diamond_db
    File alignments = diamond_blastp.aln
  }
}
```

### Advanced Usage Examples

**Custom alignment parameters:**
```wdl
call diamond_tasks.diamond_blastp {
  input:
    diamond_db = reference_db,
    query = my_proteins,
    align_id = "80",        # Higher identity threshold
    top_pct = "10",         # Keep only top 10% of hits
    query_cover = "70",     # Higher query coverage requirement
    subject_cover = "70",   # Higher subject coverage requirement
    blocksize = "0.4"       # Lower block size to save memory
    cpu_cores = 8,
    memory_gb = 8
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of reference proteomes and test data

## Testing the Module

The module includes a demonstration workflow that can be tested independently. The workflow in `testrun.wdl` automatically downloads test data and runs without requiring input files:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl
```

The test workflow (`diamond_example`) automatically:
1. Downloads reference proteome and query protein data using `ww-testdata`
2. Creates a DIAMOND database from the reference proteome
3. Performs BLASTP alignment of query proteins against the database
4. Validates all outputs
5. Generates a validation report

## Configuration Guidelines

### Resource Allocation

The module supports flexible resource configuration:

- **Memory**: 2-8 GB recommended (scales with database and query size)
- **CPUs**: 1-8 cores recommended; DIAMOND benefits significantly from multi-threading
- **Database creation**: Resource requirements depend on proteome size
  - Small bacterial proteomes: 2 GB sufficient
  - Eukaryotic proteomes: 4-8 GB recommended
  - Large metagenomic datasets: 8+ GB may be needed

### Advanced Considerations

- **Sensitivity modes**: This module uses the default BLASTP mode; DIAMOND also supports `blastx` (DNA-to-protein) and `blastn` modes not included in this module
- **Output format**: The task outputs tabular format (outfmt 6) with 12 standard BLAST columns
- **Filtering parameters**: The `align_id`, `query_cover`, and `subject_cover` parameters allow fine-tuning of alignment stringency
- **Performance**: DIAMOND is optimized for speed; for maximum sensitivity, consider adjusting the sensitivity mode in custom workflows

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support for containerized execution
- Sufficient computational resources (DIAMOND is lightweight compared to traditional BLAST)

## Performance Considerations

- **Speed**: DIAMOND provides extremely fast protein alignment, making it suitable for large-scale comparative genomics and metagenomics
- **Memory usage**: Generally modest; 2-8GB typically sufficient for most use cases
- **CPU scaling**: Performance improves linearly with additional cores
- **Storage requirements**: Minimal disk space needed; databases and alignments are compressed

## Output Description

### DIAMOND Database
- **{basename}.db.dmnd**: Binary DIAMOND database file that can be reused for multiple alignment operations

### Alignment Results

The alignment output file (`.aln.gz`) contains tab-separated values with the following 12 columns (BLAST outfmt 6):

1. **qseqid**: Query sequence ID
2. **sseqid**: Subject (reference) sequence ID
3. **pident**: Percentage of identical matches
4. **length**: Alignment length
5. **mismatch**: Number of mismatches
6. **gapopen**: Number of gap openings
7. **qstart**: Start position in query
8. **qend**: End position in query
9. **sstart**: Start position in subject
10. **send**: End position in subject
11. **evalue**: Expect value
12. **bitscore**: Bit score

The file is compressed with gzip to reduce storage requirements.

## Module Development

This module is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real protein sequence data (demonstration datasets for integration testing)
- Comprehensive validation of all outputs
- Integration testing with ww-testdata modules

For questions specific to this module or to contribute improvements, please see the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

For questions specific to DIAMOND usage or configuration, please refer to the [DIAMOND documentation](https://github.com/bbuchfink/diamond). Please make sure to cite their work if you use DIAMOND in your analyses:

Buchfink, B., Reuter, K., & Drost, H. G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods, 18(4), 366-368.

## Contributing

If you would like to contribute to this WILDS WDL module, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
