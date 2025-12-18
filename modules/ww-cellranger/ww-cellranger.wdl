## WILDS WDL for running Cell Ranger pipelines.
## Intended for use alone or as a modular component in the WILDS ecosystem.

version 1.0

task run_count {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Run Cell Ranger count for single-cell RNA-seq with GEX and ADT libraries"
    outputs: {
        results_tar: "Compressed tarball of Cell Ranger count output directory",
        web_summary: "Web summary HTML file",
        metrics_summary: "Metrics summary CSV file"
    }
  }

  parameter_meta {
    ref_gex: "GEX reference transcriptome tarball"
    ref_adt: "ADT feature reference CSV file"
    libraries_csv: "Libraries CSV file defining input data (must include paths to FASTQ files)"
    gex_fastqs: "Array of GEX FASTQ files"
    adt_fastqs: "Optional: Array of ADT FASTQ files"
    sample_id: "Sample ID for output naming"
    force_cells: "Optional: Force number of cells to be called"
    create_bam: "Generate BAM file (default: true)"
    include_introns: "Include intronic reads in count (default: true)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File ref_gex
    File ref_adt
    File libraries_csv
    Array[File] gex_fastqs
    Array[File]? adt_fastqs
    String sample_id = "results"
    Int? force_cells
    Boolean create_bam = true
    Boolean include_introns = true
    Int cpu_cores = 32
    Int memory_gb = 120
  }

  command <<<
    set -eo pipefail

    # Extract GEX reference
    mkdir -p gex_ref
    tar xf "~{ref_gex}" -C gex_ref --strip-components 1

    # Create folders and copy FASTQ files for each library type if provided
    mkdir -p gex_fastqs
    cp ~{sep=' ' gex_fastqs} gex_fastqs/

    ~{if defined(adt_fastqs) then
      "mkdir -p adt_fastqs\n    cp ~{sep=' ', select_first([adt_fastqs, []]))} adt_fastqs/"
      else ""}

    # Run cellranger count
    cellranger count \
      --id="~{sample_id}" \
      --transcriptome=gex_ref \
      ~{if create_bam then "--create-bam true" else ""} \
      --nosecondary \
      ~{if include_introns then "--include-introns true" else ""} \
      --disable-ui \
      ~{if defined(force_cells) then "--force-cells=" + force_cells else ""} \
      --libraries="~{libraries_csv}" \
      --feature-ref="~{ref_adt}" \
      --localcores=~{cpu_cores} \
      --localmem=~{memory_gb}

    # Create tarball of output directory
    tar -czf "~{sample_id}_outs.tar.gz" "~{sample_id}/outs"
  >>>

  output {
    File results_tar = "~{sample_id}_outs.tar.gz"
    File web_summary = "~{sample_id}/outs/web_summary.html"
    File metrics_summary = "~{sample_id}/outs/metrics_summary.csv"
  }

  runtime {
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
    docker: "getwilds/cellranger:8.0.1"
  }
}