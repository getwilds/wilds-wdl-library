# ww-proseq Pipeline

[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline for paired-end PRO-seq data with UMI deduplication and spike-in normalization. Adapted from the canonical [JAJ256/PROseq_alignment.sh](https://github.com/JAJ256/PROseq_alignment.sh) shell pipeline (Judd Lab / qPRO-seq protocol).

## Overview

PRO-seq (Precision Run-On sequencing) maps the 3' end of nascent RNA to identify the position of engaged RNA Pol II at single-base resolution. This pipeline wraps the standard qPRO-seq processing flow:

1. UMI-aware adapter trimming
2. rRNA depletion via alignment to an rDNA reference
3. Parallel alignment to a spike-in-merged reference and the experimental genome
4. UMI-based PCR-duplicate removal on both alignments
5. Strand-specific 3' single-base bigWig generation (one bigWig per strand, scoring the 3' end of each read)
6. Aggregated QC reporting

Spike-in alignments are kept (filtered to spike-in contigs only) so downstream analyses can derive cross-sample normalization factors.

## Pipeline Structure

**Complexity Level: Advanced** (6 distinct modules + multi-step per-sample processing)

1. **fastp UMI extraction + adapter trimming** — `ww-fastp.fastp_paired` with `umi_loc = "per_read"` (default). Extracted UMIs are appended to read names with `:` separator.
2. **rRNA depletion** — `ww-bowtie2.bowtie2_align` against the rDNA reference with `capture_unaligned = true`. The aligned BAM is discarded; the unaligned reads feed both downstream alignment steps.
3. **Spike-in alignment + chromosome filter** — `ww-bowtie2.bowtie2_align` against the merged experimental+spike-in reference with `--very-sensitive-local --no-mixed --no-discordant`, MAPQ ≥ 10, proper-pair filter. Then `ww-samtools.filter_bam_by_chrom_prefix` keeps only spike-in contigs.
4. **Experimental alignment** — `ww-bowtie2.bowtie2_align` against the experimental genome with `--sensitive-local`, MAPQ ≥ 10, proper-pair filter.
5. **UMI deduplication** — `ww-umi-tools.dedup` (`--paired`, `--umi-separator=":"`) on both BAMs.
6. **Strand-specific bigWigs** — `ww-deeptools.bam_coverage` twice per sample with `--Offset 1 --binSize 1 --normalizeUsing None` and `--samFlagInclude 82` (forward) / `98` (reverse).
7. **MultiQC** — aggregates fastp metrics across samples.

## Module Dependencies

| Module | Purpose |
|--------|---------|
| `ww-fastp` | UMI-aware adapter trimming |
| `ww-bowtie2` | All three alignment steps (rDNA, spike-in, experimental) |
| `ww-samtools` | Spike-in chromosome filter via `filter_bam_by_chrom_prefix` |
| `ww-umi-tools` | PCR duplicate removal |
| `ww-deeptools` | Strand-specific 3' bigWig generation |
| `ww-multiqc` | Aggregated QC reporting |

## Usage

### Requirements

- A WDL executor (Cromwell, miniWDL, or Sprocket)
- Docker or Apptainer/Singularity for container execution
- Network access to pull container images and (for the testrun) reference data

### Input Configuration

Three things to configure:

1. **Samples** — array of `ProseqSample` structs (`name`, `r1`, `r2`)
2. **References** — a `ProseqReferences` struct containing:
   - `experimental_fasta` — the experimental organism's genome (e.g., dm6 for Drosophila S2 cells)
   - `spikein_merged_fasta` — experimental + spike-in genomes concatenated, with a distinguishing prefix on spike-in contig names
   - `rdna_fasta` — the rDNA reference for rRNA depletion (typically the host's 45S precursor)
   - `spikein_chrom_prefix` — the prefix used on spike-in contig names in the merged FASTA
3. **Pipeline parameters** — UMI config, MAPQ threshold, resource allocation

See [`inputs.json`](inputs.json) for a template.

### Building the merged spike-in reference

The `ww-testdata.merge_fastas_with_prefix` task can build the merged FASTA programmatically:

```wdl
call ww_testdata.merge_fastas_with_prefix as build_merged { input:
  first_fasta = experimental_fasta,
  second_fasta = spikein_fasta,
  second_prefix = "hg38",
  output_name = "exp_spikein_merged"
}
```

Or build it offline once with:
```bash
cat experimental.fa <(sed 's/^>/>hg38/' spikein.fa) > merged.fa
samtools faidx merged.fa
```

### Running the Pipeline

```bash
# Using miniWDL
miniwdl run ww-proseq.wdl --input inputs.json

# Using Sprocket
sprocket run ww-proseq.wdl --inputs inputs.json --entrypoint proseq

# Using Cromwell
java -jar cromwell.jar run ww-proseq.wdl --inputs inputs.json
```

## Input Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `samples` | `Array[ProseqSample]` | required | Per-sample name + R1 + R2 |
| `references` | `ProseqReferences` | required | Three FASTAs + spike-in chromosome prefix |
| `umi_loc` | `String` | `"per_read"` | fastp UMI location. `per_read` for qPRO-seq (UMIs on both ends), `read1` / `read2` for one-sided UMIs, `""` to disable |
| `umi_len` | `Int` | `6` | UMI length in basepairs |
| `mapq_threshold` | `Int` | `10` | Minimum MAPQ for retained alignments |
| `align_cpu` | `Int` | `4` | CPU cores per bowtie2 alignment task |
| `align_memory_gb` | `Int` | `8` | Memory per bowtie2 alignment task |

### `ProseqSample` Structure

```wdl
struct ProseqSample {
    String name
    File r1
    File r2
}
```

### `ProseqReferences` Structure

```wdl
struct ProseqReferences {
    File experimental_fasta
    File spikein_merged_fasta
    File rdna_fasta
    String spikein_chrom_prefix
}
```

## Output Files

Per sample (arrays scattered over `samples`):
- `fastp_r1_trimmed`, `fastp_r2_trimmed` — UMI-extracted, adapter-trimmed FASTQs
- `fastp_html`, `fastp_json` — fastp reports
- `experimental_bam_dedup`, `experimental_bai_dedup` — deduplicated experimental BAM
- `experimental_dedup_log` — `umi_tools dedup` log (input/output read counts, PCR-dup stats)
- `spikein_bam_dedup`, `spikein_bai_dedup` — deduplicated spike-in BAM (for cross-sample normalization)
- `spikein_dedup_log` — `umi_tools dedup` log for the spike-in BAM
- `bigwig_forward`, `bigwig_reverse` — strand-specific 3' single-base bigWigs

Aggregated:
- `multiqc_report` — MultiQC HTML
- `multiqc_data` — MultiQC parsed-metrics directory

## Resource Considerations

### Compute Requirements

| Step | CPU | Memory | Notes |
|------|-----|--------|-------|
| `bowtie2_build` × 3 | 4 cores | 16 GB | Runs once per pipeline run for each reference FASTA |
| `bowtie2_align` × 3 per sample | 4-8 cores | 8-16 GB | Memory scales with reference size, especially the spike-in-merged reference |
| `umi_tools dedup` × 2 per sample | 2 cores | 8 GB | Single-threaded; memory scales with read density |
| `bam_coverage` × 2 per sample | 4 cores | 8 GB | Runs on two strand-specific bigWigs per sample |

### Optimization Tips

- The spike-in alignment is run against the **merged** reference (experimental + spike-in concatenated), then filtered. This matches the original PROseq_alignment.sh script and ensures spike-in normalization factors are computed on reads that are unambiguously spike-in.
- For very deep PRO-seq runs (>100M reads/sample), bump `align_memory_gb` to 32 GB.
- The `bam_coverage` bigWig step is the only place that does any biological math — `--Offset 1` is what makes this PRO-seq specific (3' read end = active site of RNA Polymerase II).

## Testing the Pipeline

The `testrun.wdl` is zero-config:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint proseq_example
```

It:
1. Pulls 50K read pairs from SRA accession **SRR11607571** (Judd et al. 2020 qPRO-seq, Drosophila S2 cells — same lab/protocol as the source script)
2. Downloads small dm6 (chr2L:1-1Mb) and hg38 (chr1:1-1Mb) reference fragments via `ww-testdata.download_ref_data`
3. Builds the merged spike-in reference via `ww-testdata.merge_fastas_with_prefix`
4. Pulls the human 45S rRNA precursor via `ww-testdata.download_rrna_reference`
5. Runs the full pipeline end-to-end on the synthetic-but-realistic dataset

Total testrun duration: ~10-15 minutes (dominated by the SRA download).

## Citation

If you use this pipeline, please cite the qPRO-seq protocol it implements:

> **qPRO-seq paper:** Judd J, Wojenski LA, Wainman LM, et al. A rapid, sensitive, scalable method for Precision Run-On sequencing (PRO-seq). *bioRxiv* 2020. DOI: [10.1101/2020.05.18.102277](https://doi.org/10.1101/2020.05.18.102277)

> **qPRO-seq wet-lab protocol:** [protocols.io/view/a-rapid-sensitive-scalable-method-for-precision-ru-3byl42p7ovo5](https://www.protocols.io/view/a-rapid-sensitive-scalable-method-for-precision-ru-3byl42p7ovo5/v1)

> **Reference shell pipeline:** Judd J, JAJ256/PROseq_alignment.sh. https://github.com/JAJ256/PROseq_alignment.sh

The qPRO-seq protocol is itself a dual-UMI extension of the original PRO-seq method:

> **Original PRO-seq method:** Mahat DB, Kwak H, Booth GT, et al. Base-pair-resolution genome-wide mapping of active RNA polymerases using precision nuclear run-on (PRO-seq). *Nat Protoc.* 2016;11(8):1455-1476. DOI: [10.1038/nprot.2016.086](https://doi.org/10.1038/nprot.2016.086)

Plus the underlying tool citations: fastp, bowtie2, samtools, UMI-tools, deepTools, MultiQC.

## Support

- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the library's [Contributor Guide](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md)

## License

Distributed under the MIT License. See `LICENSE` for details.
