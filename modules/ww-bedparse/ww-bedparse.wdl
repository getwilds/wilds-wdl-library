## WILDS WDL for genomic file format conversions using bedparse.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task gtf2bed {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Convert GTF annotation file to BED12 format using bedparse"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bedparse/ww-bedparse.wdl"
    outputs: {
        bed_file: "BED12-formatted annotation file"
    }
  }

  parameter_meta {
    gtf_file: "GTF annotation file to convert"
    extra_fields: "Comma separated list of extra GTF fields to add after column 12 (e.g., gene_id,gene_name)"
    filter_key: "GTF extra field on which to apply filtering"
    filter_type: "Comma separated list of filterKey field values to retain"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File gtf_file
    String? extra_fields
    String? filter_key
    String? filter_type
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -eo pipefail

    echo "Converting GTF to BED12 format using bedparse..."

    # Convert GTF to BED12 format
    bedparse gtf2bed \
      ~{if defined(extra_fields) then "--extraFields " + extra_fields else ""} \
      ~{if defined(filter_key) then "--filterKey " + filter_key else ""} \
      ~{if defined(filter_type) then "--filterType " + filter_type else ""} \
      "~{gtf_file}" > annotation.bed

    echo "Conversion complete. Total records:"
    wc -l annotation.bed
  >>>

  output {
    File bed_file = "annotation.bed"
  }

  runtime {
    docker: "getwilds/bedparse:0.2.3"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
