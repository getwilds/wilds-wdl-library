# WILDS WDL Vignettes

WILDS WDL vignettes are compact workflows that demonstrate how to combine multiple modules to create complete bioinformatics pipelines. These workflows showcase common computational patterns and serve as both functional pipelines and educational examples of modular workflow design.

## Vignette Philosophy

The WILDS WDL vignette system is designed around several key principles:

### **Module Integration**
Each vignette demonstrates how to:
- Import and combine multiple WILDS WDL modules effectively
- Pass data seamlessly between different tools
- Coordinate resource allocation across modules
- Handle complex data structures and dependencies

### **Real-World Patterns**
Vignettes focus on workflows that are:
- Commonly used in bioinformatics research
- Representative of standard analysis patterns
- Suitable for production use with real datasets
- Extensible for specific research needs

### **Educational Value**
Every vignette serves as:
- A working example of best practices
- Documentation of integration patterns
- A template for similar workflows
- A testing framework for module compatibility

## Vignette Structure

Each vignette directory contains the following standard components:

```
vignettes/vignette-name/
├── vignette-name.wdl    # Main workflow importing multiple modules
├── inputs.json          # Test inputs for the complete pipeline
└── README.md            # Vignette-specific documentation
```

### **Required Components**

- **WDL File**: Main workflow that imports and combines modules
- **Inputs Json**: Realistic test data for the complete pipeline
- **README**: Detailed documentation with usage examples and integration patterns

### **Optional Components**

- **Options Json**: Workflow execution configuration for different executors
- **Additional Input Files**: Alternative configurations for different use cases

## Available Vignettes

| Vignette | Modules Used | Description |
|----------|--------------|-------------|
| `ww-sra-star` | `ww-sra`, `ww-star` | Complete RNA-seq pipeline from SRA download to alignment |

## Using Vignettes

### **As Complete Workflows**

Vignettes can be executed as standalone pipelines:

```bash
# Navigate to vignette directory
cd vignettes/ww-sra-star

# Run with your preferred executor
miniwdl run ww-sra-star.wdl -i inputs.json
java -jar cromwell.jar run ww-sra-star.wdl --inputs inputs.json --options options.json
sprocket run ww-sra-star.wdl inputs.json
```

### **As Templates**

Vignettes serve as starting points for custom workflows:
- Copy and modify for specific research needs
- Add additional modules or analysis steps
- Adjust parameters for different datasets
- Extend with quality control or visualization components

## Testing and Validation

### **Automated Testing**
All vignettes are automatically tested through GitHub Actions:
- **Multi-executor testing**: Cromwell, miniWDL, and Sprocket
- **End-to-end validation**: Complete pipeline testing with real data
- **Integration verification**: Ensures modules work together correctly

### **Test Data**
- Small but realistic datasets for efficient testing
- Publicly available reference data
- Automated download and setup during CI

## Integration Patterns

Vignettes demonstrate key WDL patterns:
- **Module Composition**: How to effectively combine different tools
- **Data Flow Management**: Passing complex outputs between modules
- **Resource Coordination**: Balancing compute requirements across steps

## Contributing

### **Adding New Vignettes**
1. Fork the repository
2. Create a new vignette directory following the standard structure
3. Import existing modules (no new tasks should be created)
4. Add comprehensive documentation and test inputs
5. Submit a pull request demonstrating the integration pattern

### **Improving Existing Vignettes**
- Performance optimizations and resource tuning
- Enhanced documentation and usage examples
- Additional input configurations for different use cases
- Integration of new modules as they become available

## Support and Feedback

For questions about vignettes or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines
