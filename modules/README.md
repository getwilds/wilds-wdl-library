# WILDS WDL Modules

WILDS WDL modules are the building blocks of the WILDS WDL ecosystem. Each module is a tool-specific collection of reusable WDL tasks packaged with a demonstration workflow that serves both as an example and as a comprehensive test suite.

## Module Philosophy

The WILDS WDL module system is designed around several key principles:

### **Modularity and Reusability**
Each module encapsulates a specific bioinformatics tool or process, making it easy to:
- Import individual tasks into larger workflows
- Combine multiple modules to create complex pipelines
- Maintain and update tools independently
- Share common functionality across projects

### **Zero-Configuration Testing**
Every module includes a "unit test" workflow that:
- Requires no input parameters or configuration
- Automatically downloads realistic test data
- Executes each task at least once with hardcoded inputs
- Validates outputs to ensure correctness
- Serves as living documentation of proper usage
- Enables continuous integration testing

### **Standardization**
All modules follow consistent patterns for:
- File organization and naming conventions
- Input/output specifications
- Documentation structure
- Container image management
- Validation approaches

## Module Structure

Each module directory contains the following standard components:

```
modules/module-name/
├── module-name.wdl         # Main WDL file with task definitions
├── testrun.wdl             # Zero-configuration test workflow
└── README.md               # Module-specific documentation
```

### **Required Components**

- **WDL File**: Contains all task definitions for the module
- **Test Workflow**: `testrun.wdl` - Zero-configuration demonstration workflow that exercises all tasks
- **README**: Module-specific documentation with usage examples

## Available Modules

| Module | Tool | Container | Description |
|--------|------|-----------|-------------|
| `ww-annotsv` | Structural Variant Annotator | `getwilds/annotsv:3.4.4` | Annotate structural variants with AnnotSV |
| `ww-annovar` | Variant Annotator | `getwilds/annovar:GRCh38` | Annotate genetic variants with ANNOVAR |
| `ww-aws-sso` | AWS Operations | `getwilds/awscli:2.27.49` | AWS S3 operations with SSO and temporary credential support |
| `ww-bcftools` | Utilities for Variant Calls | `getwilds/bcftools:1.19` | Call and analyze variants with BCFtools |
| `ww-bedparse` | Bedparse Format Converter | `getwilds/bedparse:0.2.3` | Convert GTF annotation file to BED12 format using bedparse |
| `ww-bedtools` | Utilities for Genomic Intervals | `getwilds/bedtools:2.31.1` | Work with genomic intervals |
| `ww-bowtie` | Bowtie Short-Read Aligner | `getwilds/bowtie:1.3.1` | Build Bowtie index and align short reads to a reference genome |
| `ww-bowtie2` | Bowtie 2 Sequence Aligner | `getwilds/bowtie2:2.5.4` | Build Bowtie 2 index and align sequence reads to a reference genome |
| `ww-bwa` | BWA Aligner | `getwilds/bwa:0.7.17` | Alignment with the Burrows-Wheeler Aligner |
| `ww-cellranger` | Cell Ranger Gene Expression | `getwilds/cellranger:10.0.0` | Run cellranger count on gene expression reads from one GEM well |
| `ww-cnvkit` | CNVkit Copy Number Analysis | `getwilds/cnvkit:0.9.10` | Create CNVkit reference and run copy number analysis on tumor samples |
| `ww-colabfold` | Protein Structure Predictor | `getwilds/colabfold:1.5.5` | Predict protein structures with ColabFold (AlphaFold2 + MMseqs2) |
| `ww-consensus` | Consensus Variant Caller | `getwilds/consensus:0.1.1` | Multi-caller consensus variant integration |
| `ww-deeptools` | deepTools Coverage and Visualization | `getwilds/deeptools:3.5.6` | Analyze and visualize high-throughput sequencing data including coverage tracks and heatmaps |
| `ww-deepvariant` | DeepVariant Variant Caller | `google/deepvariant:1.10.0` | Call germline variants from aligned reads using deep learning |
| `ww-delly` | Structural Variant Caller | `getwilds/delly:1.2.9` | Call structural variants with Delly |
| `ww-deseq2` | DESeq2 Differential Expression | `getwilds/deseq2:1.40.2` | Differential expression analysis using DESeq2 |
| `ww-diamond` | DIAMOND Protein Aligner | `getwilds/diamond:2.1.16` | Create DIAMOND database and align protein sequences using BLASTP |
| `ww-ena` | ENA Data Downloader | `getwilds/ena-tools:2.1.1` | Download sequencing data from the European Nucleotide Archive (ENA) |
| `ww-fastp` | Fastp Quality Trimmer | `getwilds/fastp:1.1.0` | Quality filtering, adapter trimming, and QC reporting for FASTQ files |
| `ww-fastqc` | FastQC Quality Control | `getwilds/fastqc:0.12.1` | Run FastQC quality control analysis on FASTQ files |
| `ww-gatk` | GATK Variant Calling | `getwilds/gatk:4.6.1.0` | Variant calling and processing with GATK |
| `ww-gdc` | GDC Data Downloader | `getwilds/gdc-client:2.3.0` | Download files from GDC using manifest files or file UUIDs |
| `ww-glimpse2` | GLIMPSE2 Imputation | `getwilds/glimpse2:2.0.1` | Low-coverage whole genome sequencing imputation including reference panel preparation and phasing |
| `ww-ichorcna` | Tumor Fraction Estimator | `getwilds/ichorcna:0.2.0` | Estimate tumor fraction with ichorCNA |
| `ww-jcast` | JCAST Alternative Splicing | `getwilds/jcast:0.3.5` | Translate alternative splicing events from rMATS output into protein sequences |
| `ww-manta` | Structural Variant Caller | `getwilds/manta:1.6.0` | Call structural variants with Manta |
| `ww-megahit` | MEGAHIT Metagenome Assembler | `getwilds/megahit:1.2.9` | De novo metagenome assembly using MEGAHIT |
| `ww-rmats-turbo` | rMATS-turbo Splicing Analysis | `getwilds/rmats-turbo:latest` | Detect and quantify differential alternative splicing events from RNA-seq data |
| `ww-rseqc` | RSeQC Quality Control | `getwilds/rseqc:5.0.4` | Run comprehensive quality control metrics on aligned RNA-seq data |
| `ww-salmon` | Salmon Transcript Quantifier | `getwilds/salmon:1.10.3` | Build Salmon index and quantify transcript expression from RNA-seq reads |
| `ww-samtools` | Utilities for SAM/BAM/CRAM Files | `getwilds/samtools:1.11` | Work with Sequence Alignment/Map (SAM) format files |
| `ww-shapemapper` | ShapeMapper RNA Structure | `getwilds/shapemapper:2.3` | Analyze RNA structure probing data and generate reactivity profiles |
| `ww-sjl` | Solar Jetlag Tile Processor | `getwilds/r-utils:0.1.0` | Calculate sunrise/sunset times for geographic tiles as part of the SJL model pipeline |
| `ww-smoove` | Structural Variant Caller | `brentp/smoove:latest` | Call structural variants with Smoove |
| `ww-sourmash` | Sourmash K-mer Sketch | `getwilds/sourmash:4.8.2` | Generate sourmash sketches from BAM or FASTA files for sequence comparison |
| `ww-spades` | SPAdes Genome Assembler | `getwilds/spades:4.2.0` | De novo metagenomic assembly using metaSPAdes |
| `ww-sra` | SRA Toolkit | `getwilds/sra-tools:3.1.1` | Download sequencing data from NCBI SRA |
| `ww-star` | STAR Aligner | `getwilds/star:2.7.6a` | RNA-seq alignment with two-pass methodology |
| `ww-starling` | STARLING Ensemble Generator | `getwilds/starling:2.0.0a3` | Generate structural ensembles for intrinsically disordered proteins |
| `ww-strelka` | Strelka Variant Caller | `getwilds/strelka:2.9.10` | Run Strelka germline and somatic variant calling |
| `ww-testdata` | Test Data Downloader | `getwilds/awscli:2.27.49` | Download reference genomes and test datasets |
| `ww-tritonnp` | TritonNP Nucleosome Positioning | `python:bullseye` | Nucleosome positioning analysis from cfDNA using FFT-based fragment size analysis |
| `ww-varscan` | VarScan Variant Caller | `getwilds/varscan:2.4.6` | Run VarScan somatic and germline variant calling |

## Using Modules

### **Importing Tasks**

Modules are designed to be imported into other workflows:

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-star/ww-star.wdl" as star_tasks

workflow my_pipeline {
  call star_tasks.build_star_index { 
    input: reference_fasta = my_fasta 
  }
}
```

We recommend using GitHub URLs when importing WILDS WDL modules for a few reasons:

- **No local cloning required**: Use modules directly without downloading the repository
- **Version control**: Pin to specific commits or tags for reproducibility
- **Easy updates**: Switch between versions by changing the URL
- **Modular usage**: Import only the modules you need

### **Running Test Workflows**

Each module can be executed independently for testing or demonstration with zero configuration:

```bash
# Navigate to module directory
cd modules/ww-star

# Run test workflow with your preferred executor (no inputs needed)
miniwdl run testrun.wdl
java -jar cromwell.jar run testrun.wdl
sprocket run testrun.wdl
```

The test workflows automatically download test data and use hardcoded settings optimized for testing.

## Testing and Validation

### **Automated Testing**
All modules are automatically tested through GitHub Actions:
- **Multi-executor testing**: Cromwell, miniWDL, and Sprocket
- **Real data validation**: Uses actual bioinformatics datasets
- **Output verification**: Comprehensive validation of all outputs

## Integration with Pipelines

Modules serve as the foundation for WILDS pipelines, which demonstrate how to combine multiple modules into complete analysis workflows. See the `pipelines/` directory for examples of module integration.

## Contributing

### **Adding New Modules**
1. Fork the repository
2. Create a new module directory following the standard structure
3. Implement tasks and demonstration workflow
4. Add comprehensive tests and validation
5. Submit a pull request with detailed documentation

### **Improving Existing Modules**
- Bug fixes and improvements are always welcome
- New features should maintain backward compatibility
- All changes must pass the automated test suite

## Support and Feedback

For questions about modules or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org
- See the [contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for detailed guidelines
