# ww-pairtree Module

[![Project Status: Prototype – Useable, some support, open to feedback, unstable API.](https://getwilds.org/badges/badges/prototype.svg)](https://getwilds.org/badges/#prototype)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL module wrapping [Pairtree](https://github.com/morrislab/pairtree), a tool for reconstructing cancer evolutionary histories from multi-sample bulk DNA sequencing data.

## Overview

Pairtree infers phylogenetic relationships between cancer cell subpopulations (subclones) by analyzing mutation frequencies across multiple tissue samples from the same tumor. It first computes pairwise probability distributions over possible evolutionary relationships between mutations, then uses MCMC to sample clone trees whose subclonal frequencies best explain the observed variant allele frequencies. Pairtree scales to substantially more samples and subclones than earlier clone-tree reconstruction methods.

This module wraps the core steps of a Pairtree analysis: converting per-sample VCFs into Pairtree's SSM format, clustering mutations into subclones, sampling clone trees, and visualizing the results.

## Module Structure

This module is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and follows the standard WILDS module structure:

- **Main WDL file**: `ww-pairtree.wdl` - Contains task definitions for the module
- **Test workflow**: `testrun.wdl` - Demonstration workflow for testing and examples
- **Documentation**: This README with usage examples and parameter descriptions

## Available Tasks

### `vcf_to_ssm`

Converts single-sample SNV VCFs into a multi-sample SSM file and a matching skeleton params.json. Assumes diploid heterozygous SNVs with no copy-number correction (a fixed `var_read_prob` per input), suitable for clonal cell-line lineage analyses where purity is ~100% and ploidy is uniform across samples.

**Inputs:**
- `vcfs` (Array[File]): Single-sample VCFs containing SNV calls with AD (allelic depth) FORMAT fields, one per sample
- `sample_names` (Array[String]): Sample names, in the same order as `vcfs`, used to label the params.json samples array
- `var_read_prob` (Float, default=0.5): Expected fraction of reads supporting the variant allele given true heterozygosity
- `output_name` (String, default="converted"): Prefix for the output SSM and params.json files
- `cpu_cores` (Int, default=2): Number of CPU cores allocated for the task
- `memory_gb` (Int, default=4): Memory allocated for the task in GB
- `docker_image` (String, default=`getwilds/bcftools:1.19`): Docker image to use for this task

**Outputs:**
- `ssm_file` (File): SSM file with variant/total read counts per mutation per sample
- `params_file` (File): Skeleton params JSON with sample names and empty clusters/garbage arrays

**Note:** Mutations called in some samples but not others (private mutations) are recorded as 0 variant reads out of 0 total reads for the samples lacking that call, marking the site as uncovered/uninformative there rather than fabricating a reference genotype. This task does not adjust for differing ploidy or copy number between samples; if your samples have known karyotype differences, compute `var_read_prob` per-sample outside this task and build the SSM directly.

### `cluster_variants`

Clusters somatic mutations into subclones based on variant allele frequencies across samples (wraps `clustervars`).

**Inputs:**
- `ssm_file` (File): SSM file with variant/total read counts per mutation per sample
- `params_file` (File): Params JSON file specifying sample names (clusters/garbage may be empty)
- `model` (String, default="linfreq"): Clustering model to use (`linfreq` or `pairwise`)
- `parallel_chains` (Int, default=0): Number of Gibbs sampling chains to run in parallel via multiprocessing. Defaults to 0 (serial, no multiprocessing) since `clustervars`' multiprocessing.Manager() fails with an `AF_UNIX path too long` error under executors with deeply nested working directories, such as Cromwell
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
- `parallel_chains` (Int, default=0): Number of MCMC chains/processes to run in parallel via multiprocessing. Defaults to 0 (serial, no multiprocessing) since `pairtree`'s multiprocessing.Manager() fails with an `AF_UNIX path too long` error under executors with deeply nested working directories, such as Cromwell
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
    Array[File] sample_vcfs
    Array[String] sample_names
  }

  call pairtree_tasks.vcf_to_ssm {
    input:
      vcfs = sample_vcfs,
      sample_names = sample_names
  }

  call pairtree_tasks.cluster_variants {
    input:
      ssm_file = vcf_to_ssm.ssm_file,
      params_file = vcf_to_ssm.params_file
  }

  call pairtree_tasks.run_pairtree {
    input:
      ssm_file = vcf_to_ssm.ssm_file,
      params_file = cluster_variants.clustered_params_file
  }

  call pairtree_tasks.plot_tree {
    input:
      ssm_file = vcf_to_ssm.ssm_file,
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
**Note:** Only raise `parallel_chains` above 0 on executors that don't nest working directories deeply (e.g. Sprocket, miniWDL). Under Cromwell, `pairtree`'s multiprocessing.Manager() fails with `AF_UNIX path too long` once the container's working directory path exceeds the OS's Unix-socket path limit.

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
- **ww-testdata**: Automatic provisioning of synthetic daughter-cell-line VCF test data via `create_pairtree_vcfs`
- **Variant calling modules** (e.g. ww-gatk, ww-clair3): Somatic or germline SNV calls can be fed directly into `vcf_to_ssm`

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
1. Generates two synthetic daughter-cell-line VCFs using `ww-testdata`
2. Converts them into a multi-sample SSM file and skeleton params.json
3. Clusters the synthetic mutations into subclones
4. Samples clone trees via MCMC (reduced `trees_per_chain` for fast CI runs)
5. Generates an interactive HTML visualization of the results

## Docker Container

Most tasks in this module use the `getwilds/pairtree:1.0.1` container image, which includes:
- Pairtree and its Python dependencies (numpy, scipy, scikit-learn)
- All necessary system dependencies for clustering, tree sampling, and visualization

The `vcf_to_ssm` task instead uses `getwilds/bcftools:1.19`, since the conversion relies on `bcftools merge`/`query` rather than any Pairtree functionality.

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
