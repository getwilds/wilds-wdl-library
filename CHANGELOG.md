# Changelog

All notable changes to the WILDS WDL Library will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-01-08

We're excited to announce the first official release of the WILDS WDL Library — a comprehensive, centralized collection of reusable WDL modules and production-ready pipelines for bioinformatics research.

### Highlights

- **32 Modules** covering alignment, variant calling, structural variants, annotation, RNA-seq, QC, and more
- **8 Pipelines** ranging from basic (2-3 modules) to advanced research-grade workflows
- **Multi-executor support** — validated on Cromwell, miniWDL, and Sprocket
- **Zero-configuration testing** — every module includes a testrun that works out of the box
- **PROOF integration** — ready for Fred Hutch HPC infrastructure

### Modules

| Category | Modules |
|----------|---------|
| Alignment | `ww-bwa`, `ww-star` |
| Variant Calling | `ww-gatk`, `ww-bcftools`, `ww-varscan`, `ww-strelka` |
| Structural Variants | `ww-delly`, `ww-manta`, `ww-smoove` |
| Annotation | `ww-annovar`, `ww-annotsv`, `ww-bedparse`, `ww-bedtools` |
| Data Download | `ww-sra`, `ww-ena`, `ww-gdc` |
| RNA-seq | `ww-salmon`, `ww-rseqc`, `ww-deseq2` |
| QC & Processing | `ww-samtools`, `ww-fastqc` |
| Copy Number | `ww-ichorcna`, `ww-cnvkit` |
| Specialized | `ww-cellranger`, `ww-shapemapper`, `ww-sourmash`, `ww-spades`, `ww-megahit`, `ww-diamond` |
| Utilities | `ww-aws-sso`, `ww-testdata`, `ww-template`, `ww-tritonnp` |

### Pipelines

| Pipeline | Description | Complexity |
|----------|-------------|------------|
| `ww-bwa-gatk` | DNA alignment and variant calling | Basic |
| `ww-ena-star` | ENA download + RNA-seq alignment | Basic |
| `ww-fastq-to-cram` | FASTQ to CRAM conversion | Basic |
| `ww-sra-salmon` | SRA download + transcript quantification | Basic |
| `ww-sra-star` | SRA download + RNA-seq alignment | Basic |
| `ww-saturation` | Sequencing saturation analysis | Intermediate |
| `ww-star-deseq2` | Complete RNA-seq differential expression | Intermediate |
| `ww-leukemia` | Comprehensive genomic analysis | Advanced |

### Documentation

- [Documentation Site](https://getwilds.org/wilds-wdl-library/)
- [Fred Hutch SciWiki Article](https://sciwiki.fredhutch.org/datascience/wilds_wdl/)

[0.1.0]: https://github.com/getwilds/wilds-wdl-library/releases/tag/v0.1.0
