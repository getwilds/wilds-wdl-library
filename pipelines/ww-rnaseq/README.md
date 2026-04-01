# ww-rnaseq Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive WILDS WDL pipeline for RNA-seq analysis, covering the core workflow from raw FASTQs through differential expression with integrated quality control at every stage. Inspired by nf-core/rnaseq but focused on the essential happy path.

## Overview

This pipeline extends the simpler `ww-star-deseq2` pipeline with upstream read QC, adapter trimming, and aggregated MultiQC reporting. It integrates 7 WILDS WDL modules to perform a production-ready RNA-seq analysis workflow:

1. **Pre-trim QC** — Assess raw read quality
2. **Adapter Trimming** — Remove adapters and low-quality bases
3. **Post-trim QC** — Verify trimming improved read quality
4. **Genome Alignment** — Two-pass STAR alignment
5. **Alignment QC** — Comprehensive post-alignment quality metrics
6. **Differential Expression** — Count aggregation and DESeq2 analysis
7. **Aggregated Reporting** — MultiQC summary of all QC metrics

## Pipeline Structure

**Complexity Level: Advanced** (8 steps, 7 distinct modules)

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and demonstrates:

- **Module Integration**: Combining 7 modules for a complete end-to-end RNA-seq workflow
- **Data Flow**: Seamless passing of outputs from QC through alignment to differential expression
- **Best Practices**: Pre/post-trim QC, adapter trimming, two-pass alignment, and aggregated reporting

## Pipeline Steps

1. **FastQC — Pre-trim QC** (using `ww-fastqc` module):
   - Generates quality reports for raw FASTQ files
   - Provides baseline quality metrics before trimming

2. **Trim Galore — Adapter Trimming** (using `ww-trimgalore` module):
   - Removes adapter sequences and low-quality bases
   - Filters reads below minimum length threshold
   - Produces trimming reports with adapter detection statistics

3. **FastQC — Post-trim QC** (using `ww-fastqc` module):
   - Generates quality reports for trimmed reads
   - Confirms quality improvement after trimming

4. **STAR Index Building** (using `ww-star` module):
   - Builds STAR genome index from reference FASTA and GTF files
   - Runs once and is reused across all samples

5. **STAR Two-Pass Alignment** (using `ww-star` module):
   - Performs two-pass alignment on trimmed reads for each sample
   - Generates BAM files, BAI indices, gene counts, and alignment metrics

6. **GTF to BED Conversion** (using `ww-bedparse` module):
   - Converts GTF annotation to BED12 format for RSeQC compatibility

7. **RSeQC — Alignment QC** (using `ww-rseqc` module):
   - Analyzes read distribution, gene body coverage, and strand specificity
   - Produces comprehensive alignment quality metrics

8. **DESeq2 — Count Assembly and Differential Expression** (using `ww-deseq2` module):
   - Combines gene counts from all samples into a single matrix
   - Performs statistical analysis and generates visualizations (PCA, volcano, heatmap)

9. **MultiQC — Aggregated Reporting** (using `ww-multiqc` module):
   - Collects reports from FastQC, Trim Galore, STAR, and RSeQC
   - Generates a single interactive HTML report summarizing all QC metrics

## Module Dependencies

This pipeline imports and uses:
- **ww-fastqc module**: Pre-trim and post-trim read quality assessment (`run_fastqc` task)
- **ww-trimgalore module**: Adapter and quality trimming (`trimgalore_paired` task)
- **ww-star module**: Genome indexing and read alignment (`build_index`, `align_two_pass` tasks)
- **ww-bedparse module**: GTF to BED12 conversion (`gtf2bed` task)
- **ww-rseqc module**: RNA-seq alignment quality control (`run_rseqc` task)
- **ww-deseq2 module**: Count matrix assembly and differential expression (`combine_count_matrices`, `run_deseq2` tasks)
- **ww-multiqc module**: Aggregated QC reporting (`run_multiqc` task)

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
  "rnaseq.samples": [
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
  "rnaseq.reference_genome": {
    "name": "hg38",
    "fasta": "/path/to/genome.fasta",
    "gtf": "/path/to/annotation.gtf"
  },
  "rnaseq.reference_level": "control",
  "rnaseq.contrast": "condition,treatment,control",
  "rnaseq.star_cpu": 8,
  "rnaseq.star_memory_gb": 64
}
```

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-rnaseq.wdl --inputs inputs.json --options options.json

# Using miniWDL
miniwdl run ww-rnaseq.wdl -i inputs.json

# Using Sprocket
sprocket run ww-rnaseq.wdl inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/datademos/proof-how-to/) to submit this pipeline directly to the on-premise HPC cluster:

1. Ensure your input files are accessible by the cluster
2. Update the inputs JSON with sample information and conditions
3. Submit through the PROOF interface

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `samples` | List of SampleInfo objects with paired FASTQ files and condition labels | Array[SampleInfo] | Yes | - |
| `reference_genome` | Reference genome information (name, fasta, gtf) | RefGenome | Yes | - |
| `reference_level` | Reference level for DESeq2 contrast (e.g., 'control') | String | No | "" |
| `contrast` | DESeq2 contrast string (e.g., 'condition,treatment,control') | String | No | "" |
| `trim_quality` | Quality threshold for Trim Galore trimming | Int | No | 20 |
| `trim_length` | Minimum read length to retain after trimming | Int | No | 20 |
| `star_cpu` | Number of CPU cores for STAR alignment tasks | Int | No | 8 |
| `star_memory_gb` | Memory allocation in GB for STAR alignment tasks | Int | No | 64 |
| `genome_sa_index_nbases` | STAR SA pre-indexing string length (scales with genome size) | Int | No | 14 |

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

The pipeline produces comprehensive outputs from all modules:

| Output | Description | Source Module |
|--------|-------------|---------------|
| `pretrim_fastqc_html` | FastQC HTML reports for raw reads per sample | ww-fastqc |
| `pretrim_fastqc_zip` | FastQC ZIP archives for raw reads per sample | ww-fastqc |
| `trimgalore_r1_trimmed` | Trimmed R1 FASTQ files per sample | ww-trimgalore |
| `trimgalore_r2_trimmed` | Trimmed R2 FASTQ files per sample | ww-trimgalore |
| `trimgalore_r1_report` | Trim Galore R1 trimming reports per sample | ww-trimgalore |
| `trimgalore_r2_report` | Trim Galore R2 trimming reports per sample | ww-trimgalore |
| `posttrim_fastqc_html` | FastQC HTML reports for trimmed reads per sample | ww-fastqc |
| `posttrim_fastqc_zip` | FastQC ZIP archives for trimmed reads per sample | ww-fastqc |
| `star_bam` | Aligned BAM files per sample | ww-star |
| `star_bai` | BAM index files per sample | ww-star |
| `star_gene_counts` | Gene-level read counts per sample | ww-star |
| `star_log_final` | STAR final alignment logs per sample | ww-star |
| `star_log_progress` | STAR progress logs per sample | ww-star |
| `star_log` | STAR main logs per sample | ww-star |
| `star_sj` | Splice junction files per sample | ww-star |
| `rseqc_qc_summary` | RSeQC quality control summary per sample | ww-rseqc |
| `combined_counts_matrix` | Combined gene count matrix for all samples | ww-deseq2 |
| `sample_metadata` | Sample metadata with condition information | ww-deseq2 |
| `deseq2_all_results` | Complete DESeq2 differential expression results | ww-deseq2 |
| `deseq2_significant_results` | Filtered results with only significant genes | ww-deseq2 |
| `deseq2_normalized_counts` | DESeq2 normalized count matrix | ww-deseq2 |
| `deseq2_pca_plot` | PCA plot of samples | ww-deseq2 |
| `deseq2_volcano_plot` | Volcano plot of differential expression | ww-deseq2 |
| `deseq2_heatmap` | Heatmap of differentially expressed genes | ww-deseq2 |
| `multiqc_report` | Aggregated interactive HTML QC report | ww-multiqc |
| `multiqc_data` | MultiQC parsed metrics data directory | ww-multiqc |

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
- Trim Galore defaults (quality=20, length=20) work well for most Illumina data

## Testing the Pipeline

The pipeline includes a test workflow that uses real RNA-seq data from the DESeq2 airway dataset:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint rnaseq_example
```

The test workflow automatically:
1. Downloads small reference genome data (chromosome 1 subset, 50Mbp)
2. Downloads real RNA-seq data from SRA (4 samples from DESeq2 airway study: 2 untreated + 2 dexamethasone-treated)
3. Runs pre-trim FastQC on raw reads
4. Trims adapters with Trim Galore
5. Runs post-trim FastQC on trimmed reads
6. Builds STAR genome index
7. Performs two-pass alignment for all samples
8. Converts GTF to BED12 format for RSeQC
9. Runs RSeQC quality control on aligned BAMs
10. Combines gene count matrices from all samples
11. Performs DESeq2 differential expression analysis
12. Aggregates all QC reports with MultiQC

**Test Dataset Details:**
- Uses 4 samples from SRA (SRR1039508, SRR1039509, SRR1039512, SRR1039513)
- Real paired-end RNA-seq from human airway smooth muscle cells
- 2 untreated samples + 2 dexamethasone-treated samples
- Aligns to chromosome 1 subset for efficiency
- Resource requirements: 2 CPUs / 6GB RAM (fits GitHub Actions runners)

## Integration Patterns

This pipeline demonstrates several key WDL patterns:
- **Module Composition**: Combining 7 analysis modules in a coherent workflow
- **Data Passing**: Seamless transfer of outputs between QC, trimming, alignment, and analysis
- **Struct Usage**: Efficient organization of complex sample and reference data
- **Resource Management**: Coordinated resource allocation across modules
- **Scatter-Gather**: Parallel processing of multiple samples with aggregation
- **QC Aggregation**: Collecting reports from multiple tools into a single MultiQC summary

## Extending the Pipeline

This pipeline can be extended by:
- Adding UMI handling for deduplication
- Including pseudoalignment quantification (e.g., Salmon, Kallisto)
- Adding gene set enrichment analysis (GSEA)
- Including contamination screening
- Adding batch effect correction steps
- Integrating isoform-level analysis

## Related WILDS Components

- **ww-fastqc module**: Read quality assessment
- **ww-trimgalore module**: Adapter and quality trimming
- **ww-star module**: STAR alignment functionality
- **ww-bedparse module**: Genomic file format conversion
- **ww-rseqc module**: RNA-seq quality control
- **ww-deseq2 module**: Differential expression analysis
- **ww-multiqc module**: Aggregated QC reporting
- **ww-star-deseq2 pipeline**: Simpler RNA-seq pipeline (alignment + DE only)
- **ww-sra-star pipeline**: SRA download to alignment pipeline

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
