## WILDS WDL module for consensus variant calling.
## Combines results from multiple variant callers (HaplotypeCaller, Mutect2, bcftools)
## to generate high-confidence consensus variant calls.
##
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task consensus_processing {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Generate consensus variant calls by combining results from multiple variant callers (HaplotypeCaller, Mutect2, bcftools)"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-consensus/ww-consensus.wdl"
    outputs: {
        consensus_tsv: "Tab-separated file containing consensus variant calls with evidence from all callers"
    }
  }

  parameter_meta {
    gatk_vars: "Annotated variant table from GATK HaplotypeCaller"
    sam_vars: "Annotated variant table from samtools/bcftools"
    mutect_vars: "Annotated variant table from GATK Mutect2"
    base_file_name: "Base name for output files"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File gatk_vars
    File sam_vars
    File mutect_vars
    String base_file_name
    Int cpu_cores = 1
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Pull consensus script from GitHub
    wget -q "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/move-consensus/modules/ww-consensus/consensus-trio.R" \
      -O consensus-trio.R

    Rscript consensus-trio.R \
      "~{gatk_vars}" "~{sam_vars}" "~{mutect_vars}" "~{base_file_name}"
  >>>

  output {
    File consensus_tsv = "~{base_file_name}.consensus.tsv"
  }

  runtime {
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
    docker: "rocker/tidyverse:4.4.2"
  }
}
