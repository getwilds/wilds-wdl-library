# ww-pairtree Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module wrapping [Pairtree](https://github.com/morrislab/pairtree), a tool for reconstructing cancer evolutionary histories from multi-sample bulk DNA sequencing data.

## Overview

Pairtree infers phylogenetic relationships between cancer cell subpopulations (subclones) by analyzing mutation frequencies across multiple tissue samples from the same tumor. It first computes pairwise probability distributions over possible evolutionary relationships between mutations, then uses MCMC to sample clone trees whose subclonal frequencies best explain the observed variant allele frequencies. Pairtree scales to substantially more samples and subclones than earlier clone-tree reconstruction methods.

This module wraps the three core steps of a Pairtree analysis: clustering mutations into subclones, sampling clone trees, and visualizing the results.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-pairtree.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `cluster_variants`

Clusters somatic mutations into subclones based on variant allele frequencies across samples (wraps `clustervars`).

**Inputs:**
- `ssm_file` (File): SSM file with variant/total read counts per mutation per sample
- `params_file` (File): Params JSON file specifying sample names (clusters/garbage may be empty)
- `model` (String, default="linfreq"): Clustering model to use (`linfreq` or `pairwise`)
- `output_name` (String, default="clustered"): Prefix for the output params.json file
- `cpu_cores` (Int, default=2): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=4): Memory allocated for the task in GB
- `docker_image` (String, default=`getwilds/pairtree:1.0.1`): Docker image to use for this task

**Outputs:**
- `clustered_params_file` (File): Params JSON file with samples, mutation clusters, and garbage mutations

### `run_pairtree`

Samples clone trees consistent with observed mutation frequencies via MCMC (wraps `pairtree`).

**Inputs:**
- `ssm_file` (File): SSM file with variant/total read counts per mutation per sample
- `params_file` (File): Params JSON file with samples and mutation clusters (e.g. from `cluster_variants`)
- `output_name` (String, default="pairtree_results"): Prefix for the output results.npz file
- `trees_per_chain` (Int, default=3000): Number of MCMC tree samples to draw per chain
- `parallel_chains` (Int, default=2): Number of parallel MCMC chains/processes to use
- `phi_fitter` (String, default="projection"): Method used to fit subclonal frequencies to observed data
- `cpu_cores` (Int, default=2): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=8): Memory allocated for the task in GB
- `docker_image` (String, default=`getwilds/pairtree:1.0.1`): Docker image to use for this task

**Outputs:**
- `results_file` (File): NPZ archive containing sampled tree structures, subclonal frequencies, and log-likelihoods

### `plot_tree`

Generates an interactive HTML report visualizing the sampled clone trees (wraps `plottree`).

**Inputs:**
- `ssm_file` (File): SSM file with variant/total read counts per mutation per sample
- `params_file` (File): Params JSON file with samples and mutation clusters
- `results_file` (File): NPZ results file produced by `run_pairtree`
- `run_id` (String, default="pairtree_run"): Identifier used to label the visualization
- `cpu_cores` (Int, default=1): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=4): Memory allocated for the task in GB
- `docker_image` (String, default=`getwilds/pairtree:1.0.1`): Docker image to use for this task

**Outputs:**
- `tree_html` (File): Interactive HTML visualization of the reconstructed clone trees
- `tree_json` (File): Tree structure and frequency data exported as JSON

## Usage as a Module

### Importing into Your Workflow

```wdl
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-pairtree/ww-pairtree.wdl" as pairtree_tasks

workflow my_analysis_pipeline {
  input {
    File ssm_file
    File params_file
  }

  call pairtree_tasks.cluster_variants {
    input:
      ssm_file = ssm_file,
      params_file = params_file
  }

  call pairtree_tasks.run_pairtree {
    input:
      ssm_file = ssm_file,
      params_file = cluster_variants.clustered_params_file
  }

  call pairtree_tasks.plot_tree {
    input:
      ssm_file = ssm_file,
      params_file = cluster_variants.clustered_params_file,
      results_file = run_pairtree.results_file
  }

  output {
    File clustered_params_file = cluster_variants.clustered_params_file
    File results_file = run_pairtree.results_file
    File tree_html = plot_tree.tree_html
  }
}
```

### Advanced Usage Examples

**Custom resource allocation and MCMC settings for larger cohorts:**
```wdl
call pairtree_tasks.run_pairtree {
  input:
    ssm_file = ssm_file,
    params_file = cluster_variants.clustered_params_file,
    trees_per_chain = 10000,
    parallel_chains = 8,
    cpu_cores = 8,
    memory_gb = 32
}
```

**Using the pairwise clustering model instead of the default linfreq model:**
```wdl
call pairtree_tasks.cluster_variants {
  input:
    ssm_file = ssm_file,
    params_file = params_file,
    model = "pairwise"
}
```

### Integration Examples

This module integrates seamlessly with other WILDS components:
- **ww-testdata**: Automatic provisioning of synthetic SSM/params test data via `create_pairtree_data`
- **Variant calling modules** (e.g. ww-gatk, ww-clair3): Somatic variant calls can be converted into SSM format upstream of this module

## Testing the Module

The module includes a test workflow (`testrun.wdl`) that can be run independently:

```bash
# Using miniWDL
miniwdl run testrun.wdl

# Using Sprocket
sprocket run testrun.wdl

# Using Cromwell
java -jar cromwell.jar run testrun.wdl
```

### Automatic Demo Mode

The test workflow automatically:
1. Generates a small synthetic SSM file and params.json using `ww-testdata`
2. Clusters the synthetic mutations into subclones
3. Samples clone trees via MCMC (reduced `trees_per_chain` for fast CI runs)
4. Generates an interactive HTML visualization of the results

## Docker Container

This module uses the `getwilds/pairtree:1.0.1` container image, which includes:
- Pairtree and its Python dependencies (numpy, scipy)
- All necessary system dependencies for clustering, tree sampling, and visualization

## Citation

> Wintersinger, J.A., Dobson, S.M., Stein, L.D., Dick, J.E., Morris, Q.D. (2022). Reconstructing complex cancer evolutionary histories from multiple bulk DNA samples using Pairtree. *Blood Cancer Discovery*, 3(3), 208-219.
> DOI: [10.1158/2643-3230.BCD-21-0092](https://doi.org/10.1158/2643-3230.BCD-21-0092)

> Kulman, E., Kumar, A., Morris, Q.D. (2022). Reconstructing cancer phylogenies using Pairtree, a clone tree reconstruction algorithm. *STAR Protocols*, 3(4), 101733.
> PMID: [36129821](https://pubmed.ncbi.nlm.nih.gov/36129821/)

## Parameters and Resource Requirements

### Default Resources
- **CPU**: 1-2 cores per task
- **Memory**: 4-8 GB per task
- **Runtime**: Less than 1 minute per task for demo data; scales with number of mutations, samples, and `trees_per_chain` for real datasets

### Resource Scaling
- `trees_per_chain` and `parallel_chains`: Increase for larger cohorts or when a more thoroughly sampled posterior is needed
- `cpu_cores` / `memory_gb`: Scale with the number of mutations and samples in the SSM file; large cohorts (100+ samples, thousands of mutations) may require substantially more memory

## Related Resources

- **[Pairtree GitHub Repository](https://github.com/morrislab/pairtree)**: Source code and detailed documentation
- **[WILDS Docker Library](https://github.com/getwilds/wilds-docker-library)**: Container images used by WDL workflows
- **[WILDS Documentation](https://getwilds.org/)**: Comprehensive guides and best practices
- **[WDL Specification](https://openwdl.org/)**: Official WDL language documentation
