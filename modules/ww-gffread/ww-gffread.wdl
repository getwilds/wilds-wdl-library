## WILDS WDL Module: ww-gffread
## Wraps the `gffread` tool for GTF/GFF3 manipulation and normalization.
##
## Primary use case: bacterial GTFs from NCBI RefSeq contain mostly CDS rows
## with very few `exon` rows (typically only tRNAs/rRNAs). Downstream tools
## like STAR (GeneCounts mode) and RSeQC only consider `exon` features, so
## these bacterial GTFs silently cause 98% of protein-coding genes to drop
## out of RNA-seq results. This module's `normalize_gtf` task uses gffread
## with `--force-exons` to synthesize the missing exon features, making the
## same pipeline work for both bacterial and eukaryotic references without
## any user configuration.

version 1.0

task normalize_gtf {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Normalizes a GTF file so downstream tools see exon features for every transcript. Synthesizes exon features from CDS records for GTFs that lack them (typical of NCBI bacterial GTFs). Eukaryotic GTFs with proper exon annotations pass through unchanged."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gffread/ww-gffread.wdl"
    outputs: {
        normalized_gtf: "GTF file with exon features synthesized from CDS records where needed",
        feature_counts: "Text report comparing feature counts before and after normalization"
    }
  }

  parameter_meta {
    input_gtf: "Input GTF file to normalize. Can be any GTF - bacterial or eukaryotic."
    output_prefix: "Prefix used for output filenames"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File input_gtf
    String output_prefix = "normalized"
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -eo pipefail

    # Strip `gene`-type rows before feeding to gffread.
    awk -F'\t' '/^#/ || $3 != "gene"' "~{input_gtf}" > stripped.gtf

    # Run gffread with --force-exons to synthesize exon features from CDS
    # -T flag emits GTF (rather than gffread's default GFF3)
    gffread stripped.gtf -T --force-exons -o "~{output_prefix}.gtf"
  >>>

  output {
    File normalized_gtf = "~{output_prefix}.gtf"
    File feature_counts = "~{output_prefix}.feature_counts.txt"
  }

  runtime {
    # TODO: switch to getwilds/gffread:0.12.7 once that image is built and pushed
    docker: "quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task gff3_to_gtf {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Converts a GFF3 annotation file to GTF format using gffread. Useful when an upstream source (e.g. Ensembl Bacteria) only publishes GFF3 but downstream tools expect GTF."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gffread/ww-gffread.wdl"
    outputs: {
        gtf_file: "GTF-format annotation converted from the input GFF3"
    }
  }

  parameter_meta {
    input_gff3: "Input GFF3 file to convert"
    output_prefix: "Prefix used for output filenames"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File input_gff3
    String output_prefix = "converted"
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -eo pipefail

    # -T emits GTF instead of the default GFF3
    gffread "~{input_gff3}" -T -o "~{output_prefix}.gtf"
  >>>

  output {
    File gtf_file = "~{output_prefix}.gtf"
  }

  runtime {
    # TODO: switch to getwilds/gffread:0.12.7 once that image is built and pushed
    docker: "quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
