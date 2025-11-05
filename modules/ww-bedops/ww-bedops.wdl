## WILDS WDL for genomic file format conversions and operations using bedops.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task gtf_to_bed {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Convert GTF annotation file to BED12 format using bedops"
    outputs: {
        bed_file: "BED12-formatted annotation file"
    }
  }

  parameter_meta {
    gtf_file: "GTF annotation file to convert"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File gtf_file
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -eo pipefail

    echo "Converting GTF to BED12 format using bedops..."
    # Use bedops gtf2bed for reliable GTF to BED12 conversion
    gtf2bed < "~{gtf_file}" > annotation.bed

    echo "Conversion complete. Total records:"
    wc -l annotation.bed
  >>>

  output {
    File bed_file = "annotation.bed"
  }

  runtime {
    docker: "getwilds/bedops:2.4.42"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
