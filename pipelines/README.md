# WILDS WDL Pipelines

WILDS WDL pipelines are complete analysis workflows that combine multiple modules to create functional bioinformatics pipelines. These range from compact educational examples to comprehensive production-ready analyses, all sharing a common structure and design philosophy.

## Pipeline Philosophy

The WILDS pipeline system is designed around several key principles:

### **Module Integration**
Each pipeline demonstrates how to:
- Import and combine multiple WILDS WDL modules effectively
- Pass data seamlessly between different tools
- Coordinate resource allocation across modules

### **Real-World Patterns**
Pipelines focus on workflows that are:
- Commonly used in bioinformatics research
- Suitable for production use with real datasets
- Extensible for specific research needs

### **Educational Value**
Every pipeline serves as:
- A working example of best practices
- A template for similar workflows
- A testing framework for module compatibility

## Complexity Levels

Pipelines vary in scope and complexity. Each pipeline and its complexity level is documented below to help you choose the right starting point:

| Level | Modules | Typical Runtime | Description |
|-------|---------|-----------------|-------------|
| **Basic** | 2-3 | < 30 minutes | Simple integrations ideal for learning |
| **Intermediate** | 4-6 | 1-4 hours | Multi-step analyses for common use cases |
| **Advanced** | 10+ | > 4 hours | Comprehensive production pipelines |

## Pipeline Structure

Each pipeline directory contains the following required components:

```
pipelines/pipeline-name/
├── pipeline-name.wdl    # Main workflow importing multiple modules
├── testrun.wdl          # Zero-configuration test workflow
├── inputs.json          # Example inputs for the complete pipeline
└── README.md            # Pipeline-specific documentation
```


### **Optional Components**

- **Options JSON**: Workflow execution configuration for different WDL engines
- **Additional Input Files**: Alternative inputs for different use cases
- **Platform Configurations**: Platform-specific execution configs in subdirectories (e.g., `.cirro/` for Cirro platform integration)

## Available Pipelines

| Pipeline | Complexity | Modules Used | Description |
|----------|------------|--------------|-------------|
| `ww-bwa-gatk` | Basic | `ww-bwa`, `ww-gatk` | DNA alignment and variant calling |
| `ww-ena-star` | Basic | `ww-ena`, `ww-star` | ENA download and RNA-seq alignment |
| `ww-fastq-to-cram` | Basic | `ww-bwa`, `ww-samtools` | FASTQ to CRAM conversion |
| `ww-saturation` | Intermediate | Multiple | Sequencing saturation analysis |
| `ww-sra-salmon` | Basic | `ww-sra`, `ww-salmon` | SRA download and transcript quantification |
| `ww-sra-star` | Basic | `ww-sra`, `ww-star` | SRA download and RNA-seq alignment |
| `ww-star-deseq2` | Intermediate | `ww-star`, `ww-deseq2` | RNA-seq alignment and differential expression |
| `ww-leukemia` | Advanced | Multiple | Complete leukemia analysis pipeline |

## Using Pipelines

### **Running Directly (No Clone Required)**

You can download and run any pipeline without cloning the entire repository:

```bash
# Download a pipeline and its example inputs
# Option 1: Use curl from the command line
curl -O https://raw.githubusercontent.com/getwilds/wilds-wdl-library/main/pipelines/ww-sra-star/ww-sra-star.wdl
curl -O https://raw.githubusercontent.com/getwilds/wilds-wdl-library/main/pipelines/ww-sra-star/inputs.json
# Option 2: Download directly from GitHub by navigating to the file and clicking the download button

# Modify inputs.json as necessary for your data, then run via the command line or PROOF's point-and-click interface
sprocket run ww-sra-star.wdl inputs.json
```

This works because all pipelines import modules using GitHub URLs, so your WDL executor fetches dependencies automatically.

### **Testing and Demonstration**

Pipelines also include zero-configuration test workflows for quick demonstrations:

```bash
# Download and run a test workflow (no inputs needed)
curl -O https://raw.githubusercontent.com/getwilds/wilds-wdl-library/main/pipelines/ww-sra-star/testrun.wdl
sprocket run testrun.wdl
```

If you have the repository cloned:

```bash
# Navigate to pipeline directory
cd pipelines/ww-sra-star

# Run test workflow with your preferred executor (no inputs needed):

# Sprocket
sprocket run testrun.wdl

# Cromwell
java -jar cromwell.jar run testrun.wdl

# miniWDL
miniwdl run testrun.wdl

### **As Complete Workflows**

Pipelines can be executed as standalone workflows with your own data:

```bash
# Navigate to pipeline directory
cd pipelines/ww-sra-star

# Run with your preferred executor using the example inputs.json as a starting point:

# Sprocket
sprocket run ww-sra-star.wdl inputs.json

# Cromwell
java -jar cromwell.jar run ww-sra-star.wdl --inputs inputs.json

# miniWDL
miniwdl run ww-sra-star.wdl -i inputs.json

### **As Templates**

Pipelines serve as starting points for custom workflows:
- Copy and modify for specific research needs
- Add additional modules or analysis steps
- Adjust parameters for different datasets
- Extend with quality control or visualization components

## Testing and Validation

### **Automated Testing**
All pipelines are automatically tested through GitHub Actions:
- **Multi-executor testing**: Cromwell, miniWDL, and Sprocket
- **End-to-end validation**: Complete pipeline testing with real data
- **Integration verification**: Ensures modules work together correctly

### **Test Data**
- Small but realistic datasets for efficient testing
- Publicly available reference data
- Automatically downloaded during continuous integration

## Integration Patterns

Pipelines demonstrate key WDL patterns:
- **Module Composition**: How to effectively combine different tools
- **Data Flow Management**: Passing complex outputs between modules
- **Resource Coordination**: Balancing compute requirements across steps

## Contributing

### **Adding New Pipelines**
1. Fork the repository
2. Create a new pipeline directory following the standard structure
3. Import existing modules (prefer existing modules over creating new tasks)
4. Document the complexity level in your README
5. Add comprehensive documentation and test inputs
6. Submit a pull request demonstrating the integration pattern

### **Improving Existing Pipelines**
- Performance optimizations and resource tuning
- Enhanced documentation and usage examples
- Additional input configurations for different use cases
- Integration of new modules as they become available

## Support and Feedback

For questions about pipelines or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines
