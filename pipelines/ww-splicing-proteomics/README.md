# ww-splicing-proteomics Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline for generating custom protein databases from RNA-seq alternative splicing data for mass spectrometry proteomics analysis.

## Overview

This pipeline implements the alternative splicing proteomics workflow described in the [JCAST documentation](https://ed-lau.github.io/jcast/install.html). It integrates the `ww-star`, `ww-rmats-turbo`, and `ww-jcast` modules to perform RNA-seq alignment, differential splicing detection, and protein sequence translation.

The pipeline enables discovery of tissue-specific or condition-specific protein isoforms by:
1. Aligning RNA-seq reads to detect splice junctions
2. Identifying differential alternative splicing events between sample groups
3. Translating splice variants into protein sequences for proteomics database searching

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and demonstrates:

- **Module Integration**: Combining alignment, splicing analysis, and proteomics modules
- **Data Flow**: Seamless passing of BAM files to rMATS and splicing events to JCAST
- **Best Practices**: Modular workflow design for alternative splicing proteomics

## Pipeline Steps

1. **STAR Index Building** (using `ww-star` module):
   - Builds STAR genome index from reference FASTA and GTF files
   - Optimizes parameters for splice-aware alignment

2. **Two-Pass Alignment** (using `ww-star` module):
   - Performs STAR two-pass alignment for each sample
   - Generates sorted BAM files optimized for splice junction detection

3. **Differential Splicing Analysis** (using `ww-rmats-turbo` module):
   - Detects alternative splicing events between two sample groups
   - Identifies skipped exons (SE), mutually exclusive exons (MXE), retained introns (RI), and alternative 5'/3' splice sites (A5SS/A3SS)
   - Performs statistical testing for differential splicing

4. **Protein Sequence Translation** (using `ww-jcast` module):
   - Translates alternative splicing events into protein sequences
   - Generates tiered FASTA files based on translation confidence
   - Creates custom protein databases for mass spectrometry analysis

## Module Dependencies

This pipeline imports and uses:
- **ww-star module**: For genome indexing and read alignment (`build_index`, `align_two_pass` tasks)
- **ww-rmats-turbo module**: For differential alternative splicing detection (`rmats` task)
- **ww-jcast module**: For splice junction to protein translation (`jcast` task)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for module imports
- Sufficient compute resources for STAR alignment and rMATS analysis

### Important Notes

- **JCAST requires Ensembl-format GTF files** with `transcript_type` attributes. Standard UCSC or NCBI RefSeq GTF files will not work with JCAST.
- The pipeline requires two sample groups for differential splicing analysis (e.g., treatment vs. control)
- BAM files must be coordinate-sorted (STAR produces this by default)

### Input Configuration

Create an inputs JSON file with your samples and reference genome:

```json
{
  "splicing_proteomics.samples": [
    {
      "name": "control_1",
      "r1": "/path/to/control1_R1.fastq.gz",
      "r2": "/path/to/control1_R2.fastq.gz",
      "group": "group1"
    },
    {
      "name": "control_2",
      "r1": "/path/to/control2_R1.fastq.gz",
      "r2": "/path/to/control2_R2.fastq.gz",
      "group": "group1"
    },
    {
      "name": "treatment_1",
      "r1": "/path/to/treatment1_R1.fastq.gz",
      "r2": "/path/to/treatment1_R2.fastq.gz",
      "group": "group2"
    },
    {
      "name": "treatment_2",
      "r1": "/path/to/treatment2_R1.fastq.gz",
      "r2": "/path/to/treatment2_R2.fastq.gz",
      "group": "group2"
    }
  ],
  "splicing_proteomics.reference_genome": {
    "name": "GRCh38",
    "fasta": "/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    "gtf": "/path/to/Homo_sapiens.GRCh38.gtf"
  },
  "splicing_proteomics.ensembl_gtf": "/path/to/Homo_sapiens.GRCh38.gtf",
  "splicing_proteomics.read_length": 150,
  "splicing_proteomics.star_cpu": 8,
  "splicing_proteomics.star_memory_gb": 64
}
```

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-splicing-proteomics.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-splicing-proteomics.wdl -i inputs.json

# Using Sprocket
sprocket run ww-splicing-proteomics.wdl inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this pipeline directly to the on-premise HPC cluster:

1. Ensure your input files are accessible by the cluster
2. Update the inputs JSON with sample information and group assignments
3. Submit through the PROOF interface

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `samples` | List of SampleInfo objects with paired FASTQ files and group assignments | Array[SampleInfo] | Yes | - |
| `reference_genome` | Reference genome information (name, fasta, gtf) | RefGenome | Yes | - |
| `ensembl_gtf` | Ensembl-format GTF file required for JCAST | File | Yes | - |
| `read_length` | RNA-seq read length for rMATS analysis | Int | Yes | - |
| `read_type` | Type of reads: 'paired' or 'single' | String | No | "paired" |
| `library_type` | Library type for rMATS | String | No | "fr-unstranded" |
| `output_prefix` | Prefix for output file names | String | No | "splicing_proteomics" |
| `sjdb_overhang` | STAR splice junction overhang parameter | Int | No | 100 |
| `genome_sa_index_nbases` | STAR genome index parameter | Int | No | 14 |
| `rmats_anchor_length` | Minimum nucleotides at junction ends | Int | No | 1 |
| `rmats_novel_splice_sites` | Enable unannotated splice site detection | Boolean | No | false |
| `rmats_cstat` | Cutoff splicing difference for statistical testing | Float | No | 0.0001 |
| `jcast_min_read_count` | Minimum junction read count for translation | Int | No | 1 |
| `jcast_qvalue_min` | Minimum FDR q-value threshold | Float | No | 0 |
| `jcast_qvalue_max` | Maximum FDR q-value threshold | Float | No | 1 |
| `jcast_splice_types` | Comma-separated splice types (MXE,RI,SE,A3SS,A5SS) | String | No | "" (all) |
| `star_cpu` | CPU cores for STAR alignment | Int | No | 8 |
| `star_memory_gb` | Memory in GB for STAR alignment | Int | No | 64 |
| `rmats_cpu` | CPU cores for rMATS | Int | No | 4 |
| `rmats_memory_gb` | Memory in GB for rMATS | Int | No | 16 |
| `jcast_cpu` | CPU cores for JCAST | Int | No | 2 |
| `jcast_memory_gb` | Memory in GB for JCAST | Int | No | 8 |

### SampleInfo Structure

```json
{
  "name": "sample_name",
  "r1": "/path/to/R1.fastq.gz",
  "r2": "/path/to/R2.fastq.gz",
  "group": "group1"
}
```

Note: Samples must be assigned to either "group1" or "group2" for differential splicing analysis.

### RefGenome Structure

```json
{
  "name": "genome_name",
  "fasta": "/path/to/genome.fasta",
  "gtf": "/path/to/annotation.gtf"
}
```

## Output Files

| Output | Description | Source Module |
|--------|-------------|---------------|
| `star_bam` | Aligned BAM files for each sample | ww-star |
| `star_bai` | BAM index files for each sample | ww-star |
| `star_gene_counts` | Gene-level read counts for each sample | ww-star |
| `star_log_final` | STAR final alignment logs for each sample | ww-star |
| `star_sj` | Splice junction files for each sample | ww-star |
| `rmats_output` | Tarball containing rMATS differential splicing results | ww-rmats-turbo |
| `jcast_protein_fasta` | Combined FASTA file with translated splice variant proteins | ww-jcast |
| `jcast_output` | Tarball containing all JCAST output files | ww-jcast |

### rMATS Output Contents

The `rmats_output` tarball contains:
- `SE.MATS.JC.txt` - Skipped exon events (junction counts)
- `MXE.MATS.JC.txt` - Mutually exclusive exon events
- `RI.MATS.JC.txt` - Retained intron events
- `A5SS.MATS.JC.txt` - Alternative 5' splice site events
- `A3SS.MATS.JC.txt` - Alternative 3' splice site events
- `*.MATS.JCEC.txt` - Junction + exon count versions of each event type
- `summary.txt` - Summary statistics

### JCAST Output Contents

The `jcast_output` tarball contains tiered FASTA files:
- **T1**: High-confidence splice variants
- **T2-T4**: Lower confidence tiers
- **canonical**: Canonical protein sequences
- **orphan**: Untranslatable variants

## Resource Considerations

### Compute Requirements
- **Memory**: 64GB recommended for human genome STAR alignment
- **CPUs**: 8+ cores recommended for efficient processing
- **Storage**: Sufficient space for reference genome, STAR index, BAM files, and analysis outputs
- **Network**: Stable internet connection for module imports

### Optimization Tips
- Use call caching to save intermediate files if the pipeline crashes partway
- Adjust CPU and memory parameters based on available resources
- For testing, use a chromosome subset and reduce resource allocation
- Ensure adequate replicates per group (2+ recommended) for robust statistical analysis

## Testing the Pipeline

The pipeline includes a test workflow that uses real RNA-seq data:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint splicing_proteomics_example
```

The test workflow automatically:
1. Downloads Ensembl chr15 reference data (GTF and FASTA from JCAST test data)
2. Downloads real RNA-seq data from SRA (4 samples from airway study)
3. Runs the complete splicing proteomics pipeline
4. Generates protein FASTA files from detected splice variants

**Test Dataset Details:**
- Uses 4 samples from SRA (SRR1039508, SRR1039509, SRR1039512, SRR1039513)
- Real paired-end RNA-seq from human airway smooth muscle cells ( a subset of 500,000 reads per sample)
- 2 untreated samples (group1) + 2 dexamethasone-treated samples (group2)
- Aligns to Ensembl chr15 for efficiency

## Use Cases

This pipeline is designed for:
- **Proteogenomics**: Creating sample-specific protein databases for mass spectrometry
- **Alternative splicing analysis**: Identifying condition-specific splice variants
- **Biomarker discovery**: Finding differentially spliced proteins between conditions
- **Cancer proteomics**: Detecting tumor-specific splice isoforms

## Related WILDS Components

- **ww-star module**: STAR alignment functionality
- **ww-rmats-turbo module**: rMATS differential splicing analysis
- **ww-jcast module**: JCAST protein sequence translation
- **ww-star-deseq2 pipeline**: RNA-seq differential expression analysis

## References

- STAR: Dobin et al. (2013) Bioinformatics
- rMATS: Shen et al. (2014) PNAS
- JCAST: Lau et al. (2019) J Proteome Res

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
