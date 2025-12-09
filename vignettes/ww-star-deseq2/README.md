# ww-star-deseq2 Vignette
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL vignette demonstrating a complete RNA-seq analysis pipeline from alignment to differential expression using modular WILDS components.

## Overview

This vignette showcases how to combine WILDS WDL modules to create a comprehensive RNA-seq analysis workflow. It integrates the `ww-star`, `ww-rseqc`, `ww-bedparse`, and `ww-deseq2` modules to perform alignment, quality control, and differential expression analysis.

This vignette serves as both a functional workflow and a demonstration of modular WDL design patterns within the WILDS ecosystem.

## Vignette Structure

This vignette is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and demonstrates:

- **Module Integration**: Combining multiple modules for a complete analysis pipeline
- **Data Flow**: Seamless passing of outputs between alignment, QC, and analysis modules
- **Best Practices**: Modular workflow design patterns for RNA-seq analysis

## Pipeline Steps

1. **STAR Index Building** (using `ww-star` module):
   - Builds STAR genome index from reference FASTA and GTF files
   - Optimizes parameters for the reference genome

2. **GTF to BED Conversion** (using `ww-bedparse` module):
   - Converts GTF annotation to BED12 format for RSeQC compatibility
   - Preserves gene structure information including exons and CDS

3. **Two-Pass Alignment** (using `ww-star` module):
   - Performs STAR two-pass alignment for each sample
   - Generates BAM files, BAI indices, gene counts, and alignment metrics

4. **Quality Control** (using `ww-rseqc` module):
   - Runs comprehensive RNA-seq quality control metrics
   - Analyzes read distribution, gene body coverage, and other QC metrics

5. **Count Matrix Assembly** (using `ww-deseq2` module):
   - Combines gene counts from all samples into a single matrix
   - Creates sample metadata file with condition information

6. **Differential Expression Analysis** (using `ww-deseq2` module):
   - Performs DESeq2 statistical analysis
   - Generates normalized counts, statistical results, and visualizations

## Module Dependencies

This vignette imports and uses:
- **ww-star module**: For genome indexing and read alignment (`build_index`, `align_two_pass` tasks)
- **ww-bedparse module**: For GTF to BED12 conversion (`gtf2bed` task)
- **ww-rseqc module**: For RNA-seq quality control (`run_rseqc` task)
- **ww-deseq2 module**: For count matrix assembly and differential expression analysis (`combine_count_matrices`, `run_deseq2` tasks)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for module imports
- Sufficient compute resources for STAR alignment and DESeq2 analysis

### Input Configuration

Create an inputs JSON file with your samples and reference genome:

```json
{
  "star_deseq2.samples": [
    {
      "name": "sample1",
      "r1": "/path/to/sample1_R1.fastq.gz",
      "r2": "/path/to/sample1_R2.fastq.gz",
      "condition": "control"
    },
    {
      "name": "sample2",
      "r1": "/path/to/sample2_R1.fastq.gz",
      "r2": "/path/to/sample2_R2.fastq.gz",
      "condition": "treatment"
    }
  ],
  "star_deseq2.reference_genome": {
    "name": "hg38",
    "fasta": "/path/to/genome.fasta",
    "gtf": "/path/to/annotation.gtf"
  },
  "star_deseq2.reference_level": "control",
  "star_deseq2.contrast": "condition,treatment,control",
  "star_deseq2.star_cpu": 8,
  "star_deseq2.star_memory_gb": 64
}
```

### Running the Vignette

```bash
# Using Cromwell
java -jar cromwell.jar run ww-star-deseq2.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-star-deseq2.wdl -i inputs.json

# Using Sprocket
sprocket run ww-star-deseq2.wdl inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this workflow directly to the on-premise HPC cluster:

1. Ensure your input files are accessible by the cluster
2. Update the inputs JSON with sample information and conditions
3. Submit through the PROOF interface

### Running on Cirro

This vignette includes [Cirro](https://cirro.bio/) platform configuration files in the [.cirro/](.cirro/) directory for cloud execution. To run this workflow on Cirro:

1. Ensure you have access to a Cirro instance
2. Follow the [Cirro pipeline documentation](https://docs.cirro.bio/pipelines/adding-pipelines/) to add this workflow to your Cirro instance
3. The configuration files in `.cirro/` define:
   - Input form fields for the Cirro web interface
   - Compute resource requirements
   - Output file handling and organization
   - Data preprocessing steps

For detailed information on configuring and using Cirro pipelines, see the [official Cirro documentation](https://docs.cirro.bio/).

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `samples` | List of SampleInfo objects with paired FASTQ files and condition labels | Array[SampleInfo] | Yes | - |
| `reference_genome` | Reference genome information (name, fasta, gtf) | RefGenome | Yes | - |
| `reference_level` | Reference level for DESeq2 contrast (e.g., 'control') | String | No | "" |
| `contrast` | DESeq2 contrast string (e.g., 'condition,treatment,control') | String | No | "" |
| `star_cpu` | Number of CPU cores for STAR alignment tasks | Int | No | 8 |
| `star_memory_gb` | Memory allocation in GB for STAR alignment tasks | Int | No | 64 |

### SampleInfo Structure

```json
{
  "name": "sample_name",
  "r1": "/path/to/R1.fastq.gz",
  "r2": "/path/to/R2.fastq.gz",
  "condition": "control_or_treatment"
}
```

### RefGenome Structure

```json
{
  "name": "genome_name",
  "fasta": "/path/to/genome.fasta",
  "gtf": "/path/to/annotation.gtf"
}
```

## Output Files

The vignette produces comprehensive outputs from all modules:

| Output | Description | Source Module |
|--------|-------------|---------------|
| `star_bam` | Aligned BAM files for each sample | ww-star |
| `star_bai` | BAM index files for each sample | ww-star |
| `star_gene_counts` | Gene-level read counts for each sample | ww-star |
| `star_log_final` | STAR final alignment logs for each sample | ww-star |
| `star_log_progress` | STAR progress logs for each sample | ww-star |
| `star_log` | STAR main logs for each sample | ww-star |
| `star_sj` | Splice junction files for each sample | ww-star |
| `rseqc_qc_summary` | Quality control summary reports for each sample | ww-rseqc |
| `combined_counts_matrix` | Combined gene count matrix for all samples | ww-deseq2 |
| `sample_metadata` | Sample metadata file with condition information | ww-deseq2 |
| `deseq2_all_results` | Complete DESeq2 differential expression results | ww-deseq2 |
| `deseq2_significant_results` | Filtered results with only significant genes | ww-deseq2 |
| `deseq2_normalized_counts` | DESeq2 normalized count matrix | ww-deseq2 |
| `deseq2_pca_plot` | Principal Component Analysis plot | ww-deseq2 |
| `deseq2_volcano_plot` | Volcano plot of differential expression | ww-deseq2 |
| `deseq2_heatmap` | Heatmap of differentially expressed genes | ww-deseq2 |

## Resource Considerations

### Compute Requirements
- **Memory**: 64GB recommended for human genome alignment (6GB minimum for chr1 testing)
- **CPUs**: 8+ cores recommended for efficient processing (2 cores minimum for testing)
- **Storage**: Sufficient space for reference genome, STAR index, BAM files, and analysis outputs
- **Network**: Stable internet connection for module imports and SRA downloads

### Optimization Tips
- Use call caching to save intermediate files if the pipeline crashes partway
- Adjust `star_cpu` and `star_memory_gb` parameters based on available resources
- For testing on limited hardware (e.g., GitHub Actions runners), use chr1 subset and reduce CPU/memory allocation
- Ensure adequate replicates per condition (3+ recommended) for robust DESeq2 analysis
- The workflow supports resource customization via `star_cpu` and `star_memory_gb` input parameters

## Testing the Vignette

The vignette includes a test workflow that uses real RNA-seq data from the DESeq2 airway dataset:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint star_deseq2_example
```

The test workflow automatically:
1. Downloads small reference genome data (chromosome 1 subset)
2. Downloads real RNA-seq data from SRA (4 samples from DESeq2 airway study: 2 untreated + 2 dexamethasone-treated)
3. Builds STAR genome index
4. Performs two-pass alignment for all samples
5. Converts GTF to BED12 format for RSeQC compatibility
6. Runs RSeQC quality control on aligned BAM files
7. Combines gene count matrices from all samples
8. Performs DESeq2 differential expression analysis
9. Generates PCA plots, volcano plots, and heatmaps
10. Validates all outputs

**Test Dataset Details:**
- Uses 4 samples from SRA (SRR1039508, SRR1039509, SRR1039512, SRR1039513)
- Real paired-end RNA-seq from human airway smooth muscle cells
- 2 untreated samples + 2 dexamethasone-treated samples
- Aligns to chromosome 1 subset for efficiency
- Resource requirements: 2 CPUs / 6GB RAM (fits GitHub Actions runners)

The vignette is automatically tested as part of the WILDS WDL Library CI/CD pipeline using:
- Multiple WDL executors (Cromwell, miniWDL, Sprocket)
- Real RNA-seq data that validates RSeQC and DESeq2 functionality
- Comprehensive output validation including plots and statistical results
- Performance benchmarking

## Integration Patterns

This vignette demonstrates several key WDL patterns:
- **Module Composition**: Combining multiple analysis modules effectively
- **Data Passing**: Seamless transfer of outputs between modules
- **Struct Usage**: Efficient organization of complex sample and reference data
- **Resource Management**: Coordinated resource allocation across modules
- **Scatter-Gather**: Parallel processing of multiple samples with aggregation

## Extending the Vignette

This vignette can be extended by:
- Adding upstream quality control (e.g., FastQC, Trim Galore)
- Including additional QC tools (e.g., Qualimap, featureCounts)
- Adding gene set enrichment analysis (GSEA)
- Integrating isoform-level quantification (e.g., Salmon, Kallisto)
- Including visualization and reporting modules (e.g., MultiQC)
- Adding batch effect correction steps

## Related WILDS Components

- **ww-star module**: STAR alignment functionality
- **ww-rseqc module**: RNA-seq quality control
- **ww-bedparse module**: Genomic file format conversion
- **ww-deseq2 module**: Differential expression analysis
- **ww-sra-star vignette**: SRA download to alignment pipeline
- **Other vignettes**: Additional integration examples

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL vignette, please see our [WILDS Contributor Guide](https://getwilds.org/guide/) and the [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
