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

### **Comprehensive Testing**
Every module includes a "unit test" workflow that:
- Executes each task at least once with realistic inputs
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
├── module-name.wdl         # Main WDL file with tasks and demo workflow
├── inputs.json             # Test inputs for the demo workflow
└── README.md               # Module-specific documentation
```

### **Required Components**

- **WDL File**: Contains all tasks plus a demonstration workflow
- **Inputs JSON**: Provides realistic test data for automated testing
- **README**: Module-specific documentation with usage examples

### **Optional Components**

- **Options JSON**: Workflow execution configuration
- **Additional Input Files**: Alternative input configurations for testing

## Available Modules

| Module | Tool/Process | Container | Description |
|--------|--------------|-----------|-------------|
| `ww-sra` | SRA Toolkit | `getwilds/sra-tools:3.1.1` | Download sequencing data from NCBI SRA |
| `ww-star` | STAR Aligner | `getwilds/star:2.7.6a` | RNA-seq alignment with two-pass methodology |

## Using Modules

### **Importing Tasks**

Modules are designed to be imported into other workflows:

```wdl
import "path/to/modules/ww-star/ww-star.wdl" as star_tasks

workflow my_pipeline {
  call star_tasks.build_star_index { 
    input: reference_fasta = my_fasta 
  }
}
```

### **Running Demo Workflows**

Each module can be executed independently for testing or demonstration:

```bash
# Navigate to module directory
cd modules/ww-star

# Run with your preferred executor
miniwdl run ww-star.wdl -i inputs.json
java -jar cromwell.jar run ww-star.wdl --inputs inputs.json --options options.json
sprocket run ww-star.wdl inputs.json
```

## Testing and Validation

### **Automated Testing**
All modules are automatically tested through GitHub Actions:
- **Multi-executor testing**: Cromwell, miniWDL, and Sprocket
- **Real data validation**: Uses actual bioinformatics datasets
- **Output verification**: Comprehensive validation of all outputs

## Integration with Vignettes

Modules serve as the foundation for WILDS vignettes, which demonstrate how to combine multiple modules into complete analysis pipelines. See the `vignettes/` directory for examples of module integration.

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
