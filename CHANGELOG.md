# Changelog

All notable changes to the WILDS WDL Library will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2026-07-13

Third release of the WILDS WDL Library, expanding into single-cell genomics, structural biology, and improved HPC/infrastructure tooling.

### Highlights

- **53 Modules** (up from 44): 9 new modules added
- **15 Pipelines** (up from 12): 3 new pipelines added
- **HPC Test Run Infrastructure**: modules and pipelines can now be validated directly on Fred Hutch HPC via automated test runs

### Added

- New modules: `ww-cellbender`, `ww-clair3`, `ww-esmfold`, `ww-gffread`, `ww-multiqc`, `ww-seurat`, `ww-trimgalore`, `ww-umi-tools`, `ww-viennarna`
- New pipelines: `ww-proseq`, `ww-rnaseq`, `ww-sra-cellranger`
- HPC-specific test run workflow (`testrun_hpc.wdl`) for modules and pipelines
- Automated HPC test run GitHub Action with PR comment reporting
- Zenodo DOI citation support
- NGC key file support in `ww-sra` and related pipelines
- EDAM ontology metadata for Dockstore entries
- Repo stats GitHub Action
- WILDS WDL agent scaffolding for AI-assisted module development
- Cirro HTTP download support in `ww-ena-star`
- AI disclosure in README and CONTRIBUTING

### Changed

- Docker images moved from `runtime` block to task `input` block for all modules
- CellRanger `run_count` tasks now support graceful failure mode
- `ww-sra-cellranger` pipeline updated to accept Run Selector `.txt` file as input and include CellBender post-processing
- `ww-rnaseq` pipeline updated with collaborator feedback and GTF processing fixes
- Updated `sprocket` version and cleaned up doc post-processing
- Added concurrency limit to Cromwell CI/CD configuration
- Modules registered as tools in Dockstore with updated author info

## [0.2.0] - 2026-03-31

Second release of the WILDS WDL Library, featuring significant growth in both content and contributor community.

### Highlights

- **44 Modules** (up from 33) — 11 new modules added
- **12 Pipelines** (up from 8) — 4 new pipelines added
- **12 Contributors** — grown from a core team to a collaborative project driven by the [WILDS WDL Development Program](https://sciwiki.fredhutch.org/datascience/wilds_workflow_dev/)
- **Cirro Integration** — multiple pipelines now include cloud execution configurations for [Cirro](https://cirro.bio/)

### Added

- New modules: `ww-bowtie`, `ww-bowtie2`, `ww-colabfold`, `ww-deeptools`, `ww-deepvariant`, `ww-fastp`, `ww-glimpse2`, `ww-jcast`, `ww-rmats-turbo`, `ww-sjl`, `ww-starling`
- New pipelines: `ww-imputation`, `ww-jetlag`, `ww-splicing-proteomics`, `ww-starling-batch`
- Cirro platform configurations for `ww-bwa-gatk`, `ww-fastq-to-cram`, `ww-imputation`, `ww-saturation`, `ww-splicing-proteomics`, `ww-sra-salmon`, `ww-sra-star`
- Dockstore integration via `.dockstore.yml`
- CITATION.cff for proper citation support
- Automated Cirro configuration validation in CI/CD
- Documentation site auto-generation via Sprocket and GitHub Pages

### Changed

- Single-end alignment support in `ww-star`, `ww-salmon`, and related pipelines
- Updated execution engine versions in CI/CD
- Multi-sample phasing support in `ww-glimpse2` and `ww-imputation`
- Organization references updated from DaSL to OCDO

### Fixed

- Memory allocation fixes for `ww-gatk` tasks (mutect2, haplotype_caller)
- Import URL consistency across modules and pipelines
- Author ORCID corrections in Dockstore configuration
- Various bug fixes in `ww-glimpse2`, `ww-imputation`, and `ww-cellranger`

## [0.1.0] - 2025-01-09

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
| Assembly | `ww-spades`, `ww-megahit` |
| Specialized | `ww-cellranger`, `ww-shapemapper`, `ww-sourmash`, `ww-diamond`, `ww-tritonnp` |
| Utilities | `ww-aws-sso`, `ww-testdata` |

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

[0.3.0]: https://github.com/getwilds/wilds-wdl-library/releases/tag/v0.3.0
[0.2.0]: https://github.com/getwilds/wilds-wdl-library/releases/tag/v0.2.0
[0.1.0]: https://github.com/getwilds/wilds-wdl-library/releases/tag/v0.1.0
