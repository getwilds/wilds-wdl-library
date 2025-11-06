## WILDS WDL for genomic file format conversions using UCSC utilities.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task gtf_to_bed {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Convert GTF annotation file to BED12 format using UCSC tools"
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

    echo "Converting GTF to BED12 format using UCSC tools..."
    # Convert GTF to GenePred format, then to BED12
    # This produces proper BED12 format compatible with tools like RSeQC
    gtfToGenePred -genePredExt "~{gtf_file}" /dev/stdout | \
      genePredToBed stdin annotation.bed

    echo "Conversion complete. Total records:"
    wc -l annotation.bed
  >>>

  output {
    File bed_file = "annotation.bed"
  }

  runtime {
    docker: "quay.io/biocontainers/ucsc-gtftogenepred:447--h2a80c09_3"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
