# ww-rnaseq Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive WILDS WDL pipeline for RNA-seq analysis, covering the core workflow from raw FASTQs through differential expression with integrated quality control at every stage. Inspired by nf-core/rnaseq but focused on the core workflow steps.

## Overview

This pipeline extends the simpler `ww-star-deseq2` pipeline with upstream read QC, adapter trimming, and aggregated MultiQC reporting. It integrates 8 WILDS WDL modules to perform a production-ready RNA-seq analysis workflow:

1. **GTF Normalization** — Ensure annotation has proper exon features (critical for bacterial genomes)
2. **Pre-trim QC** — Assess raw read quality
3. **Adapter Trimming** — Remove adapters and low-quality bases
4. **Post-trim QC** — Verify trimming improved read quality
5. **Genome Alignment** — Two-pass STAR alignment
6. **Alignment QC** — Comprehensive post-alignment quality metrics
7. **Differential Expression** — Count aggregation and DESeq2 analysis
8. **Aggregated Reporting** — MultiQC summary of all QC metrics

## Pipeline Structure

**Complexity Level: Advanced** (12 steps, 8 distinct modules + 1 pipeline-local task)

1. **GTF Normalization** (using `ww-gffread` module):
   - Ensures the reference GTF has proper `exon` features for every transcript
   - Critical for bacterial NCBI GTFs which only have CDS rows — without this step, STAR and RSeQC silently ignore 98% of protein-coding genes
   - Eukaryotic GTFs with proper exon annotations pass through unchanged

2. **FastQC — Pre-trim QC** (using `ww-fastqc` module):
   - Generates quality reports for raw FASTQ files
   - Provides baseline quality metrics before trimming

3. **Trim Galore — Adapter Trimming** (using `ww-trimgalore` module):
   - Removes adapter sequences and low-quality bases
   - Filters reads below minimum length threshold
   - Produces trimming reports with adapter detection statistics

4. **FastQC — Post-trim QC** (using `ww-fastqc` module):
   - Generates quality reports for trimmed reads
   - Confirms quality improvement after trimming

5. **STAR Index Building** (using `ww-star` module):
   - Builds STAR genome index from reference FASTA and normalized GTF
   - Runs once and is reused across all samples

6. **STAR Two-Pass Alignment** (using `ww-star` module):
   - Performs two-pass alignment on trimmed reads for each sample
   - Generates BAM files, BAI indices, gene counts, and alignment metrics

7. **GTF to BED Conversion** (using `ww-bedparse` module):
   - Converts normalized GTF to BED12 format for RSeQC compatibility

8. **RSeQC — Alignment QC** (using `ww-rseqc` module):
   - Analyzes read distribution, gene body coverage, and strand specificity
   - Produces comprehensive alignment quality metrics

9. **DESeq2 — Count Assembly and Differential Expression** (using `ww-deseq2` module):
   - Combines gene counts from all samples into a single matrix
   - Performs statistical analysis and generates visualizations (PCA, volcano, heatmap)

10. **DESeq2 — Compiled Results** (using `ww-deseq2` module):
    - Merges differential expression results with normalized counts
    - Annotates genes with descriptions from the original (non-normalized) GTF (gene_name, gene_biotype, product, locus_tag, db_xref where available)
    - Walks every feature row and aggregates per gene_id, so it works on GTFs without `gene` rows (e.g. UCSC ncbiRefSeq) and on NCBI GTFs that use the `gene "name"` attribute instead of `gene_name`
    - Produces a single comprehensive CSV for downstream use

11. **MultiQC — Aggregated Reporting** (using `ww-multiqc` module):
    - Collects reports from FastQC, Trim Galore, STAR, and RSeQC
    - Generates a single interactive HTML report summarizing all QC metrics

12. **Organize Outputs** (optional, pipeline-local task):
    - Reorganizes all outputs into a clean numbered directory structure with sample-named subdirectories
    - Packages results as a single `.tar.gz` for easy delivery
    - Disabled by default (`organize_results = false`) to avoid data duplication
    - BAM files and trimmed FASTQs can be optionally included via `include_bams` and `include_trimmed_fastqs`

## Module Dependencies

This pipeline imports and uses:
- **ww-gffread module**: GTF normalization to ensure exon features exist (`normalize_gtf` task)
- **ww-fastqc module**: Pre-trim and post-trim read quality assessment (`run_fastqc` task)
- **ww-trimgalore module**: Adapter and quality trimming (`trimgalore_paired` task)
- **ww-star module**: Genome indexing and read alignment (`build_index`, `align_two_pass` tasks)
- **ww-bedparse module**: GTF to BED12 conversion (`gtf2bed` task)
- **ww-rseqc module**: RNA-seq alignment quality control (`run_rseqc` task)
- **ww-deseq2 module**: Count matrix assembly, differential expression, and compiled results (`combine_count_matrices`, `run_deseq2`, `compile_deseq2_results` tasks)
- **ww-multiqc module**: Aggregated QC reporting (`run_multiqc` task)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Internet access for module imports
- Sufficient compute resources for STAR alignment and DESeq2 analysis

### Input Configuration

Sample information can be provided either as a **TSV file** (recommended) or as **separate arrays** in the inputs JSON.

#### Option 1: TSV File (Recommended)

Create a tab-separated file (`samples.tsv`) with a header row and one row per sample:

```tsv
name	r1	r2	condition
sample1	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz	control
sample2	/path/to/sample2_R1.fastq.gz	/path/to/sample2_R2.fastq.gz	control
sample3	/path/to/sample3_R1.fastq.gz	/path/to/sample3_R2.fastq.gz	treatment
sample4	/path/to/sample4_R1.fastq.gz	/path/to/sample4_R2.fastq.gz	treatment
```

Then reference it in your inputs JSON:

```json
{
  "rnaseq.samples_tsv": "/path/to/samples.tsv",
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

An example `samples.tsv` template is included in this directory.

#### Option 2: Separate Arrays

Alternatively, provide sample information as parallel arrays directly in the inputs JSON. This approach is primarily used by the `testrun.wdl` test workflow, where constructing a TSV file within WDL itself is not straightforward.

```json
{
  "rnaseq.sample_names": ["sample1", "sample2", "sample3", "sample4"],
  "rnaseq.r1_fastqs": ["/path/to/sample1_R1.fastq.gz", "/path/to/sample2_R1.fastq.gz", "/path/to/sample3_R1.fastq.gz", "/path/to/sample4_R1.fastq.gz"],
  "rnaseq.r2_fastqs": ["/path/to/sample1_R2.fastq.gz", "/path/to/sample2_R2.fastq.gz", "/path/to/sample3_R2.fastq.gz", "/path/to/sample4_R2.fastq.gz"],
  "rnaseq.conditions": ["control", "control", "treatment", "treatment"],
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
| `samples_tsv` | TSV file with columns: name, r1, r2, condition (header required) | File | No* | - |
| `sample_names` | Array of sample names | Array[String] | No* | - |
| `r1_fastqs` | Array of R1 FASTQ file paths | Array[File] | No* | - |
| `r2_fastqs` | Array of R2 FASTQ file paths | Array[File] | No* | - |
| `conditions` | Array of condition labels per sample | Array[String] | No* | - |
| `reference_genome` | Reference genome information (name, fasta, gtf) | RefGenome | Yes | - |
| `reference_level` | Reference level for DESeq2 contrast (e.g., 'control') | String | No | "" |
| `contrast` | DESeq2 contrast string (e.g., 'condition,treatment,control') | String | No | "" |
| `trim_quality` | Quality threshold for Trim Galore trimming | Int | No | 20 |
| `trim_length` | Minimum read length to retain after trimming | Int | No | 20 |
| `star_cpu` | Number of CPU cores for STAR alignment tasks | Int | No | 8 |
| `star_memory_gb` | Memory allocation in GB for STAR alignment tasks | Int | No | 64 |
| `genome_sa_index_nbases` | STAR SA pre-indexing string length (scales with genome size) | Int | No | 14 |
| `min_counts` | Minimum counts a gene must have to pass DESeq2 pre-filtering | Int | No | 10 |
| `min_samples` | A gene must meet the `min_counts` threshold in this many samples to be kept. 0 = Gene is kept if its count across all samples meets the `min_count` threshold (default: 0)| Int | No | 0 |
| `shrinkage_method` | LFC shrinkage method: `apeglm`, `ashr`, `normal`, or empty for none | String | No | "" |
| `organize_results` | Reorganize outputs into a clean directory tarball | Boolean | No | false |
| `include_bams` | Include BAM/BAI files in the organized tarball | Boolean | No | false |
| `include_trimmed_fastqs` | Include trimmed FASTQ files in the organized tarball | Boolean | No | false |

*Provide either `samples_tsv` OR all four separate arrays (`sample_names`, `r1_fastqs`, `r2_fastqs`, `conditions`).

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
| `deseq2_ma_plot` | MA plot of log fold change vs. mean expression | ww-deseq2 |
| `deseq2_ma_plot_shrunk` | MA plot with shrunken LFC (empty if shrinkage not applied) | ww-deseq2 |
| `deseq2_results_shrunk` | DESeq2 results with shrunken LFC (empty if shrinkage not applied) | ww-deseq2 |
| `deseq2_compiled_results` | Combined CSV with DESeq2 stats, normalized counts, and gene annotations (gene_name, gene_biotype, product, locus_tag, db_xref where available) | ww-deseq2 |
| `multiqc_report` | Aggregated interactive HTML QC report | ww-multiqc |
| `multiqc_data` | MultiQC parsed metrics data directory | ww-multiqc |
| `organized_results` | (Optional) Tarball with all outputs organized by step and sample | pipeline-local |

## Resource Considerations

### Compute Requirements
- **Memory**: 64GB recommended for human genome alignment (6GB minimum for chr1 testing)
- **CPUs**: 8+ cores recommended for efficient processing (2 cores minimum for testing)
- **Storage**: Sufficient space for reference genome, STAR index, BAM files, and analysis outputs
- **Network**: Stable internet connection for module imports and SRA downloads

### Optimization Tips
- Use call caching to save intermediate files if the pipeline crashes partway
- Adjust `star_cpu` and `star_memory_gb` parameters based on available resources

## Testing the Pipeline

The pipeline includes a test workflow that uses real RNA-seq data referenced in the DESeq2 airway dataset:

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
3. Normalizes the reference GTF (ensures exon features exist for all transcripts)
4. Runs pre-trim FastQC on raw reads
5. Trims adapters with Trim Galore
6. Runs post-trim FastQC on trimmed reads
7. Builds STAR genome index using normalized GTF
8. Performs two-pass alignment for all samples
9. Converts normalized GTF to BED12 format for RSeQC
10. Runs RSeQC quality control on aligned BAMs
11. Combines gene count matrices from all samples
12. Performs DESeq2 differential expression analysis
13. Aggregates all QC reports with MultiQC

**Test Dataset Details:**
- Uses 4 samples from SRA (SRR1039508, SRR1039509, SRR1039512, SRR1039513)
- Real paired-end RNA-seq from human airway smooth muscle cells
- 2 untreated samples + 2 dexamethasone-treated samples
- Aligns to chromosome 1 subset for efficiency
- Resource requirements: 2 CPUs / 6GB RAM (fits GitHub Actions runners)

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

Special thanks to [Elle Glenny](https://github.com/eglenny) for extensive testing of this pipeline, identifying bugs, and suggesting feature improvements that have shaped the current design. Thank you for your contributions!

If you would like to contribute to this WILDS WDL pipeline, please see our [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.

