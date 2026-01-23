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
| `ww-annovar` | Variant Annotator | `getwilds/annovar:GRCh38` | Annotate genetic variants with ANNOVAR |
| `ww-annotsv` | Structural Variant Annotator | `getwilds/annotsv:3.4.4` | Annotate structural variants with AnnotSV |
| `ww-aws-sso` | AWS Operations | `getwilds/awscli:2.27.49` | AWS S3 operations with SSO and temporary credential support |
| `ww-bcftools` | Utilities for Variant Calls | `getwilds/bcftools:1.19` | Call and analyze variants with BCFtools |
| `ww-bedtools` | Utilities for Genomic Intervals | `getwilds/bedtools:2.31.1` | Work with genomic intervals |
| `ww-bwa` | BWA Aligner | `getwilds/bwa:0.7.17` | Alignment with the Burrows-Wheeler Aligner |
| `ww-consensus` | Consensus Variant Caller | `getwilds/consensus:0.1.1` | Multi-caller consensus variant integration |
| `ww-delly` | Structural Variant Caller | `getwilds/delly:1.2.9` | Call structural variants with Delly |
| `ww-gatk` | GATK Variant Calling | `getwilds/gatk:4.6.1.0` | Variant calling and processing with GATK |
| `ww-ichorcna` | Tumor Fraction Estimator | `getwilds/ichorcna:0.2.0` | Estimate tumor fraction with ichorCNA |
| `ww-manta` | Structural Variant Caller | `getwilds/manta:1.6.0` | Call structural variants with Manta |
| `ww-samtools` | Utilities for SAM/BAM/CRAM Files | `getwilds/samtools:1.11` | Work with Sequence Alignment/Map (SAM) format files |
| `ww-smoove` | Structural Variant Caller | `brentp/smoove:latest` | Call structural variants with Smoove |
| `ww-sra` | SRA Toolkit | `getwilds/sra-tools:3.1.1` | Download sequencing data from NCBI SRA |
| `ww-star` | STAR Aligner | `getwilds/star:2.7.6a` | RNA-seq alignment with two-pass methodology |
| `ww-testdata` | Test Data Downloader | `getwilds/awscli:2.27.49` | Download reference genomes and test datasets |

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
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines
