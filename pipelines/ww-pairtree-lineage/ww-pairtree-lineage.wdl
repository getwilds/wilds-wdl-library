## WILDS WDL pipeline for clonal lineage inference from single-nucleotide variants.
## Converts per-sample SNV VCFs (e.g. from a set of related cell lines) into Pairtree's
## SSM format, clusters mutations into subclones, samples clone trees via MCMC, and
## renders an interactive visualization -- answering questions like whether two daughter
## lines arose from the same common ancestor clone.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-pairtree/ww-pairtree.wdl" as ww_pairtree

workflow pairtree_lineage {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "WDL pipeline inferring clonal lineage relationships among related samples (e.g. daughter cell lines) from single-nucleotide variant VCFs, using Pairtree to cluster mutations into subclones and reconstruct their evolutionary tree"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-pairtree-lineage/ww-pairtree-lineage.wdl"
    outputs: {
      ssm_file: "SSM file with variant/total read counts per mutation per sample",
      clustered_params_file: "Params JSON file with samples, mutation clusters, and garbage mutations",
      results_file: "NPZ archive containing sampled tree structures, subclonal frequencies, and log-likelihoods",
      tree_html: "Interactive HTML visualization of the reconstructed clone trees",
      tree_json: "Tree structure and frequency data exported as JSON"
    }
  }

  parameter_meta {
    sample_vcfs: "Array of single-sample SNV VCFs, one per sample (e.g. daughter cell line), each with AD (allelic depth) FORMAT fields"
    sample_names: "Sample names, in the same order as sample_vcfs"
    var_read_prob: "Expected fraction of reads supporting the variant allele given true heterozygosity (0.5 for diploid heterozygous SNVs with no copy-number correction, the standard assumption for clonal cell lines)"
    clustering_model: "Clustering model to use in cluster_variants (linfreq or pairwise)"
    trees_per_chain: "Number of MCMC tree samples to draw per chain in run_pairtree"
    cpu_cores: "Number of CPU cores allocated per task"
    memory_gb: "Memory allocated per task in GB"
  }

  input {
    Array[File] sample_vcfs
    Array[String] sample_names
    Float var_read_prob = 0.5
    String clustering_model = "linfreq"
    Int trees_per_chain = 3000
    Int cpu_cores = 2
    Int memory_gb = 8
  }

  call ww_pairtree.vcf_to_ssm { input:
      vcfs = sample_vcfs,
      sample_names = sample_names,
      var_read_prob = var_read_prob,
      cpu_cores = cpu_cores,
      memory_gb = memory_gb
  }

  call ww_pairtree.cluster_variants { input:
      ssm_file = vcf_to_ssm.ssm_file,
      params_file = vcf_to_ssm.params_file,
      model = clustering_model,
      cpu_cores = cpu_cores,
      memory_gb = memory_gb
  }

  call ww_pairtree.run_pairtree { input:
      ssm_file = vcf_to_ssm.ssm_file,
      params_file = cluster_variants.clustered_params_file,
      trees_per_chain = trees_per_chain,
      cpu_cores = cpu_cores,
      memory_gb = memory_gb
  }

  call ww_pairtree.plot_tree { input:
      ssm_file = vcf_to_ssm.ssm_file,
      params_file = cluster_variants.clustered_params_file,
      results_file = run_pairtree.results_file,
      cpu_cores = cpu_cores,
      memory_gb = memory_gb
  }

  output {
    File ssm_file = vcf_to_ssm.ssm_file
    File clustered_params_file = cluster_variants.clustered_params_file
    File results_file = run_pairtree.results_file
    File tree_html = plot_tree.tree_html
    File tree_json = plot_tree.tree_json
  }
}
