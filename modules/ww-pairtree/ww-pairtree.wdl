## WILDS WDL module for Pairtree.
## Reconstructs cancer evolutionary trees from multi-sample bulk DNA sequencing data
## by clustering mutations into subclones and sampling clone trees consistent with
## observed variant allele frequencies.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

#### TASK DEFINITIONS ####
# Define tasks for each functionality of the tool represented by this module

task vcf_to_ssm {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Converts single-sample SNV VCFs into a multi-sample SSM file and a matching skeleton params.json consumable by cluster_variants. Assumes diploid heterozygous SNVs with no copy-number correction (var_read_prob fixed per input), suitable for clonal cell-line lineage analyses."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-pairtree/ww-pairtree.wdl"
    outputs: {
        ssm_file: "SSM file with variant/total read counts per mutation per sample",
        params_file: "Skeleton params JSON with sample names and empty clusters/garbage arrays"
    }
    topic: "oncology,genomics"
    species: "human"
    operation: "format_conversion"
    input_sample_required: "vcfs:variant_calling:vcf"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "ssm_file:variant_calling:tabular,params_file:variant_calling:json"
    output_reference: "none"
  }

  parameter_meta {
    vcfs: "Array of single-sample VCFs containing SNV calls with AD (allelic depth) FORMAT fields, one per sample"
    sample_names: "Sample names, in the same order as vcfs, used to label the params.json samples array"
    var_read_prob: "Expected fraction of reads supporting the variant allele given true heterozygosity (0.5 for diploid heterozygous SNVs with no copy-number correction)"
    output_name: "Prefix for the output SSM and params.json files"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
    docker_image: "Docker image to use for this task"
  }

  input {
    Array[File] vcfs
    Array[String] sample_names
    Float var_read_prob = 0.5
    String output_name = "converted"
    Int cpu_cores = 2
    Int memory_gb = 4
    String docker_image = "getwilds/bcftools:1.19"
  }

  command <<<
    set -eo pipefail

    # Compress and index each single-sample VCF so bcftools merge can join them
    vcf_array=(~{sep=' ' vcfs})
    compressed_vcfs=()
    for vcf in "${vcf_array[@]}"; do
      base="$(basename "$vcf")"
      bcftools view -O z -o "${base}.gz" "$vcf"
      bcftools index -t "${base}.gz"
      compressed_vcfs+=("${base}.gz")
    done

    # Join all samples on shared loci; sites absent from a given sample are
    # recorded as missing ("./.") genotypes rather than dropped
    bcftools merge -m none "${compressed_vcfs[@]}" -O v -o merged.vcf

    # Build the SSM header
    printf "id\tname\tvar_reads\ttotal_reads\tvar_read_prob\n" > "~{output_name}.ssm"

    # Extract per-sample allelic depths and assemble SSM rows. A missing AD
    # (private mutation not called in that sample) is recorded as 0 variant
    # reads out of 0 total reads -- an uninformative, uncovered site rather
    # than a fabricated reference call.
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' merged.vcf | \
      awk -v OFS='\t' -v vrp="~{var_read_prob}" '
      BEGIN { vrp = vrp + 0 }
      {
        var_str = ""
        tot_str = ""
        vrp_str = ""
        for (i = 5; i <= NF; i++) {
          split($i, ad, ",")
          if ($i == "." || ad[1] == "" || ad[2] == "") {
            v = 0; t = 0
          } else {
            v = ad[2] + 0
            t = ad[1] + ad[2]
          }
          sep = (i == 5) ? "" : ","
          var_str = var_str sep v
          tot_str = tot_str sep t
          vrp_str = vrp_str sep vrp
        }
        printf "s%d\t%s_%s_%s_%s\t%s\t%s\t%s\n", NR-1, $1, $2, $3, $4, var_str, tot_str, vrp_str
      }' >> "~{output_name}.ssm"

    # Build a skeleton params.json: sample names only, ready for cluster_variants
    sample_json=$(printf '"%s",' ~{sep=' ' sample_names} | sed 's/,$//')
    printf '{"samples": [%s], "clusters": [], "garbage": []}\n' "$sample_json" > "~{output_name}.params.json"
  >>>

  output {
    File ssm_file = "~{output_name}.ssm"
    File params_file = "~{output_name}.params.json"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task cluster_variants {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Clusters somatic mutations into subclones based on variant allele frequencies across samples, producing a params.json consumable by run_pairtree"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-pairtree/ww-pairtree.wdl"
    outputs: {
        clustered_params_file: "Params JSON file with samples, mutation clusters, and garbage mutations"
    }
    topic: "oncology,genomics"
    species: "human"
    operation: "clustering"
    input_sample_required: "ssm_file:variant_calling:tabular"
    input_sample_optional: "none"
    input_reference_required: "params_file:variant_calling:json"
    input_reference_optional: "none"
    output_sample: "clustered_params_file:variant_calling:json"
    output_reference: "none"
  }

  parameter_meta {
    ssm_file: "SSM file with variant/total read counts per mutation per sample"
    params_file: "Params JSON file specifying sample names (clusters/garbage may be empty)"
    model: "Clustering model to use (linfreq or pairwise)"
    parallel_chains: "Number of Gibbs sampling chains to run in parallel via multiprocessing. 0 runs chains serially with no multiprocessing, which avoids clustervars' multiprocessing.Manager() Unix-socket path-length failures under executors with deeply nested working directories (e.g. Cromwell); increase for large real datasets on executors without that issue."
    output_name: "Prefix for the output params.json file"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
    docker_image: "Docker image to use for this task"
  }

  input {
    File ssm_file
    File params_file
    String model = "linfreq"
    Int parallel_chains = 0
    String output_name = "clustered"
    Int cpu_cores = 2
    Int memory_gb = 4
    String docker_image = "getwilds/pairtree:1.0.1"
  }

  command <<<
    set -eo pipefail

    clustervars \
      --model ~{model} \
      --parallel ~{parallel_chains} \
      "~{ssm_file}" \
      "~{params_file}" \
      "~{output_name}.params.json"
  >>>

  output {
    File clustered_params_file = "~{output_name}.params.json"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task run_pairtree {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Samples clone trees consistent with observed mutation frequencies via MCMC, producing a posterior distribution over cancer evolutionary histories"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-pairtree/ww-pairtree.wdl"
    outputs: {
        results_file: "NPZ archive containing sampled tree structures, subclonal frequencies, and log-likelihoods"
    }
    topic: "oncology,genomics"
    species: "human"
    operation: "phylogenetic_analysis"
    input_sample_required: "ssm_file:variant_calling:tabular"
    input_sample_optional: "none"
    input_reference_required: "params_file:variant_calling:json"
    input_reference_optional: "none"
    output_sample: "results_file:phylogenetic_tree:npz"
    output_reference: "none"
  }

  parameter_meta {
    ssm_file: "SSM file with variant/total read counts per mutation per sample"
    params_file: "Params JSON file with samples and mutation clusters (e.g. from cluster_variants)"
    output_name: "Prefix for the output results.npz file"
    trees_per_chain: "Number of MCMC tree samples to draw per chain"
    parallel_chains: "Number of MCMC chains/processes to run in parallel via multiprocessing. 0 runs chains serially with no multiprocessing, which avoids pairtree's multiprocessing.Manager() Unix-socket path-length failures under executors with deeply nested working directories (e.g. Cromwell); increase for large real datasets on executors without that issue."
    phi_fitter: "Method used to fit subclonal frequencies to observed data"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
    docker_image: "Docker image to use for this task"
  }

  input {
    File ssm_file
    File params_file
    String output_name = "pairtree_results"
    Int trees_per_chain = 3000
    Int parallel_chains = 0
    String phi_fitter = "projection"
    Int cpu_cores = 2
    Int memory_gb = 8
    String docker_image = "getwilds/pairtree:1.0.1"
  }

  command <<<
    set -eo pipefail

    pairtree \
      --params "~{params_file}" \
      --trees-per-chain ~{trees_per_chain} \
      --parallel ~{parallel_chains} \
      --phi-fitter ~{phi_fitter} \
      "~{ssm_file}" \
      "~{output_name}.results.npz"
  >>>

  output {
    File results_file = "~{output_name}.results.npz"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task plot_tree {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Generates an interactive HTML report visualizing the sampled clone trees, subclonal frequencies, and pairwise mutation relationships"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-pairtree/ww-pairtree.wdl"
    outputs: {
        tree_html: "Interactive HTML visualization of the reconstructed clone trees",
        tree_json: "Tree structure and frequency data exported as JSON"
    }
    topic: "oncology,genomics"
    species: "human"
    operation: "visualization"
    input_sample_required: "ssm_file:variant_calling:tabular,results_file:phylogenetic_tree:npz"
    input_sample_optional: "none"
    input_reference_required: "params_file:variant_calling:json"
    input_reference_optional: "none"
    output_sample: "tree_html:phylogenetic_tree:html,tree_json:phylogenetic_tree:json"
    output_reference: "none"
  }

  parameter_meta {
    ssm_file: "SSM file with variant/total read counts per mutation per sample"
    params_file: "Params JSON file with samples and mutation clusters"
    results_file: "NPZ results file produced by run_pairtree"
    run_id: "Identifier used to label the visualization"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
    docker_image: "Docker image to use for this task"
  }

  input {
    File ssm_file
    File params_file
    File results_file
    String run_id = "pairtree_run"
    Int cpu_cores = 1
    Int memory_gb = 4
    String docker_image = "getwilds/pairtree:1.0.1"
  }

  command <<<
    set -eo pipefail

    plottree \
      --runid "~{run_id}" \
      --tree-json "~{run_id}.tree.json" \
      "~{ssm_file}" \
      "~{params_file}" \
      "~{results_file}" \
      "~{run_id}.results.html"
  >>>

  output {
    File tree_html = "~{run_id}.results.html"
    File tree_json = "~{run_id}.tree.json"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
