# ww-pairtree-lineage Pipeline
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A WILDS WDL pipeline for inferring clonal lineage relationships among related samples from single-nucleotide variant (SNV) VCFs, using [Pairtree](https://github.com/morrislab/pairtree).

## Overview

This pipeline answers questions like: *did two daughter cell lines arise from the same common ancestor clone, and if so, what does their shared evolutionary history look like?* It takes per-sample SNV VCFs, converts them into Pairtree's SSM format, clusters mutations into subclones, samples clone trees consistent with the observed mutation frequencies via MCMC, and renders an interactive visualization of the result.

Shared mutations that appear at consistent frequency across related samples indicate inheritance from a common ancestor; the reconstructed tree topology shows whether samples are siblings descending from a shared internal node, whether one is a direct descendant of another, or whether they show no meaningful shared ancestry. This makes the pipeline well suited to clonal cell-line lineage tracing, in addition to its original use case of reconstructing cancer subclone phylogenies from multi-region or multi-timepoint tumor sampling.

**Complexity level**: Basic (1 module, 4 sequential tasks)

## Pipeline Structure

This pipeline is part of the [WILDS WDL Library](https://github.com/getwilds/wilds-wdl-library) and contains:

- **Main WDL file**: `ww-pairtree-lineage.wdl` - Workflow definition
- **Test workflow**: `testrun.wdl` - Demonstration workflow using test data
- **Example inputs**: `inputs.json` - Starting point for input configuration

## Module Dependencies

This pipeline imports and uses:
- **ww-pairtree module**: For VCF-to-SSM conversion, mutation clustering, clone tree sampling, and visualization (`vcf_to_ssm`, `cluster_variants`, `run_pairtree`, `plot_tree` tasks)

## Usage

### Requirements

- WDL-compatible workflow executor (Cromwell, miniWDL, Sprocket, etc.)
- Docker/Apptainer support
- Per-sample SNV VCFs with AD (allelic depth) FORMAT fields

### Assumptions

This pipeline assumes diploid heterozygous SNVs with no copy-number correction (`var_read_prob` fixed at 0.5 for every mutation in every sample), which is the standard assumption for clonal cell lines at ~100% purity and uniform ploidy. If your samples have known karyotype differences, build the SSM file directly with per-sample `var_read_prob` values rather than using this pipeline's `vcf_to_ssm` step.

Mutations called in some samples but not others are treated as uncovered/uninformative in the samples lacking that call (0 variant reads out of 0 total reads), not as reference genotype calls.

### Input Configuration

```json
{
  "pairtree_lineage.sample_vcfs": ["/path/to/daughter_line_A.vcf", "/path/to/daughter_line_B.vcf"],
  "pairtree_lineage.sample_names": ["daughter_line_A", "daughter_line_B"],
  "pairtree_lineage.var_read_prob": 0.5,
  "pairtree_lineage.clustering_model": "linfreq",
  "pairtree_lineage.trees_per_chain": 3000,
  "pairtree_lineage.cpu_cores": 2,
  "pairtree_lineage.memory_gb": 8
}
```

`sample_vcfs` and `sample_names` must be given in the same order. Add a known ancestral/founder line's VCF as an additional sample when available -- it anchors the tree and makes "common ancestor" claims more interpretable than inferring ancestry from the daughter lines alone.

### Running the Pipeline

```bash
# Using Cromwell
java -jar cromwell.jar run ww-pairtree-lineage.wdl --inputs inputs.json

# Using miniWDL
miniwdl run ww-pairtree-lineage.wdl -i inputs.json

# Using Sprocket
sprocket run ww-pairtree-lineage.wdl @inputs.json
```

### For Fred Hutch Users

Fred Hutch users can use [PROOF](https://sciwiki.fredhutch.org/datademos/proof-how-to/) to submit this pipeline directly to the on-premise HPC cluster.

## Input Parameters

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `sample_vcfs` | Array of single-sample SNV VCFs, one per sample, each with AD FORMAT fields | Array[File] | Yes | - |
| `sample_names` | Sample names, in the same order as `sample_vcfs` | Array[String] | Yes | - |
| `var_read_prob` | Expected fraction of reads supporting the variant allele given true heterozygosity | Float | No | `0.5` |
| `clustering_model` | Clustering model to use (`linfreq` or `pairwise`) | String | No | `"linfreq"` |
| `trees_per_chain` | Number of MCMC tree samples to draw per chain | Int | No | `3000` |
| `cpu_cores` | Number of CPU cores allocated per task | Int | No | `2` |
| `memory_gb` | Memory allocated per task in GB | Int | No | `8` |

## Output Files

| Output | Description |
|--------|-------------|
| `ssm_file` | SSM file with variant/total read counts per mutation per sample |
| `clustered_params_file` | Params JSON file with samples, mutation clusters, and garbage mutations |
| `results_file` | NPZ archive containing sampled tree structures, subclonal frequencies, and log-likelihoods |
| `tree_html` | Interactive HTML visualization of the reconstructed clone trees |
| `tree_json` | Tree structure and frequency data exported as JSON |

## Interpreting Results

Inspect `tree_html` (or `tree_json`'s `phi`/`struct` fields) for the pattern of shared vs. private mutations across samples:
- **Shared common ancestor**: mutations present in both samples at consistent clonal frequency, with each sample also carrying its own private mutations -- shown as a branch point with the samples as siblings (or descendants) below a shared internal node
- **One descended from the other**: one sample's mutation set is a strict subset of the other's -- shown as a direct parent-child relationship in the tree
- **No shared ancestry**: little to no meaningfully shared derived mutation signal -- samples hang independently off the root

## Testing the Pipeline

The pipeline includes a test workflow that can be run independently:

```bash
# Using Cromwell
java -jar cromwell.jar run testrun.wdl

# Using miniWDL
miniwdl run testrun.wdl --entrypoint pairtree_lineage_example

# Using Sprocket
sprocket run testrun.wdl
```

The test workflow generates two synthetic daughter-cell-line VCFs (sharing 2 SNVs from a common ancestor, plus 1 private SNV each) via `ww-testdata`, then runs the full pipeline against them with reduced `trees_per_chain` for fast CI runs.

## Citation

> Wintersinger, J.A., Dobson, S.M., Stein, L.D., Dick, J.E., Morris, Q.D. (2022). Reconstructing complex cancer evolutionary histories from multiple bulk DNA samples using Pairtree. *Blood Cancer Discovery*, 3(3), 208-219.
> DOI: [10.1158/2643-3230.BCD-21-0092](https://doi.org/10.1158/2643-3230.BCD-21-0092)

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Office of the Chief Data Officer (OCDO) at wilds@fredhutch.org, or open an issue on the [WILDS WDL Library issue tracker](https://github.com/getwilds/wilds-wdl-library/issues).

## Contributing

If you would like to contribute to this WILDS WDL pipeline, please see our [WILDS WDL Library contributing guidelines](https://github.com/getwilds/wilds-wdl-library/blob/main/.github/CONTRIBUTING.md) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
