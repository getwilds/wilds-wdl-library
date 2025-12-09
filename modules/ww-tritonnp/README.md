# ww-tritonnp Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module for running [TritonNP](https://github.com/caalo/TritonNP) for nucleosome positioning analysis.

## Overview

This module provides reusable WDL tasks for running **TritonNP**, which generates phasing features using FFT-based fragment size analysis. This module is designed to be a modular component in the WILDS ecosystem, suitable for integration into larger bioinformatics pipelines.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Tasks**: `triton_main`, `combine_fms`
- **Test workflow**: `testrun.wdl` (demonstration workflow executing all tasks)
- **Container**: `python:bullseye`

## Tasks

### `triton_main`

Runs TritonNP on a single sample to generate phasing feature matrices.

**Inputs:**
- `sample_name` (String): Sample name
- `bam_path` (File): BAM file
- `bam_index_path` (File): BAM index file
- `bias_path` (File): GC corrected file from Griffin
- `annotation` (File): BED file of genomic region to process on
- `reference_genome` (File): Reference genome file
- `reference_genome_index` (File): Reference genome file index
- `results_dir` (String): Output directory name
- `map_quality` (Int): Mapping quality threshold as a positive integer
- `size_range` (String): Size range as a space-delimited string, such as '15 500'
- `cpus` (Int): Number of CPUs to use
- `plot_list` (File): File containing names of genes to plot

**Outputs:**
- `fm_file` (File): Phasing feature matrix output file

### `combine_fms`

Combines phasing feature matrices from multiple samples.

**Inputs:**
- `fm_files` (Array[File]): Array of output files from TritonNP
- `results_dir` (String): Output directory name

**Outputs:**
- `final` (File): Aggregated output file from TritonNP

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-tritonnp/ww-tritonnp.wdl" as tritonnp

workflow my_cfdna_analysis {
  input {
    Array[String] sample_names
    Array[File] bam_files
    Array[File] bam_indices
    Array[File] bias_files
    File annotation_bed
    File reference_fasta
    File reference_fai
    File plot_list
  }

  scatter (i in range(length(sample_names))) {
    call tritonnp.triton_main {
      input:
        sample_name = sample_names[i],
        bam_path = bam_files[i],
        bam_index_path = bam_indices[i],
        bias_path = bias_files[i],
        annotation = annotation_bed,
        reference_genome = reference_fasta,
        reference_genome_index = reference_fai,
        results_dir = "tritonnp_output",
        map_quality = 20,
        size_range = "15 500",
        cpus = 4,
        plot_list = plot_list
    }
  }

  call tritonnp.combine_fms {
    input:
      fm_files = triton_main.fm_file,
      results_dir = "combined_results"
  }

  output {
    Array[File] individual_features = triton_main.fm_file
    File combined_features = combine_fms.final
  }
}
```

### Advanced Usage Examples

**Custom fragment size range and mapping quality:**
```wdl
call tritonnp.triton_main {
  input:
    sample_name = "patient_001",
    bam_path = aligned_bam,
    bam_index_path = aligned_bai,
    bias_path = gc_bias_file,
    annotation = promoter_regions_bed,
    reference_genome = hg38_fasta,
    reference_genome_index = hg38_fai,
    results_dir = "output",
    map_quality = 30,
    size_range = "100 220",  # Focus on nucleosome-sized fragments
    cpus = 8,
    plot_list = genes_of_interest
}
```

### Integration Examples

This module pairs well with other WILDS modules:
- **ww-bwa**: For aligning sequencing reads to generate input BAM files
- **ww-samtools**: For BAM file processing and quality control

## Testing the Module

The module includes a test workflow ([testrun.wdl](testrun.wdl)) that automatically downloads test data and runs without requiring input files:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl --entrypoint tritonnp_example

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Downloads test data using `ww-testdata`
2. Processes the test sample with TritonNP
3. Combines results into a final output file

## Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support

## Performance Considerations

### Default Resources
- **CPU**: User-specified
- **Memory**: 4 GB per task
- **Runtime**: Varies based on input size and CPU allocation

## Citation

If you use this module in your research, please cite:

> TritonNP
> https://github.com/caalo/TritonNP

## Additional Resources

- **[TritonNP GitHub Repository](https://github.com/caalo/TritonNP)**: Source code and documentation

## Contributing

To improve this module or report issues:
1. Fork the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library)
2. Make your changes following WILDS conventions
3. Test thoroughly with the demonstration workflow
4. Submit a pull request with detailed documentation

## Support and Feedback

For questions about this module or to report issues:
- Open an issue in the [WILDS WDL Library repository](https://github.com/getwilds/wilds-wdl-library/issues)
- Contact the Fred Hutch Data Science Lab at wilds@fredhutch.org
- See the [WILDS Contributor Guide](https://getwilds.org/guide/) for detailed guidelines

## Related Resources

- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation

## License

Distributed under the MIT License. See `LICENSE` for details.
