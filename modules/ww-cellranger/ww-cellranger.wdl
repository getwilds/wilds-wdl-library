## WILDS WDL for running Cell Ranger pipelines.
## Intended for use alone or as a modular component in the WILDS ecosystem.

version 1.0

task run_count {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Run cellranger count on gene expression reads from one GEM well"
    outputs: {
        results_tar: "Compressed tarball of Cell Ranger count output directory",
        web_summary: "Web summary HTML file",
        metrics_summary: "Metrics summary CSV file"
    }
  }

  parameter_meta {
    gex_fastqs: "Paired GEX FASTQs with naming convention: SampleName_S1_L001_R1_001.fastq.gz"
    ref_gex: "GEX reference transcriptome tarball"
    sample_id: "Sample ID for output naming"
    create_bam: "Generate BAM file (default: true)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
    expect_cells: "Optional: Expected number of recovered cells"
  }

  input {
    Array[File] gex_fastqs
    File ref_gex
    String sample_id
    Boolean create_bam = true
    Int cpu_cores = 8
    Int memory_gb = 64
    Int? expect_cells
  }

  command <<<
    set -eo pipefail

    # Extract GEX reference
    mkdir -p gex_ref
    tar xf "~{ref_gex}" -C gex_ref --strip-components 1

    # Create folder and copy FASTQ files
    mkdir -p gex_fastqs
    cp ~{sep=' ' gex_fastqs} gex_fastqs/
    FASTQS=$(pwd)/gex_fastqs

    mkdir -p "~{sample_id}_outs"

    # Run cellranger count
    cellranger count \
      --transcriptome=gex_ref \
      --fastqs="$FASTQS" \
      --localcores=~{cpu_cores} \
      --localmem=~{memory_gb} \
      --output-dir="~{sample_id}" \
      ~{"--expect-cells=" + expect_cells} \
      --id="~{sample_id}" \
      --create-bam="~{create_bam}"

    # Create tarball of output directory
    tar -czf "~{sample_id}_outs.tar.gz" "~{sample_id}/outs"

    # Move output files to working directory for outputting
    mv "~{sample_id}/outs/web_summary.html" .
    mv "~{sample_id}/outs/metrics_summary.csv" .
  >>>

  output {
    File results_tar = "~{sample_id}_outs.tar.gz"
    File web_summary = "web_summary.html"
    File metrics_summary = "metrics_summary.csv"
  }

  runtime {
    docker: "getwilds/cellranger:10.0.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}


task prepare_fastqs {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Rename a pair of FASTQs to convention: SampleName_S1_L001_R1_001.fastq.gz"
    outputs: {
        renamed_fastqs: "FASTQ files renamed to Cell Ranger convention"
    }
  }

  parameter_meta {
    r1_fastqs: "Array of R1 FASTQ files"
    r2_fastqs: "Array of R2 FASTQ files"
    sample_name: "Sample name for FASTQ naming"
  }

  input {
    Array[File] r1_fastqs
    Array[File] r2_fastqs
    String sample_name
  }

  command <<<
    set -eo pipefail

    # Create folder for FASTQ files
    mkdir -p renamed_fastqs

    # Create arrays from WDL inputs
    R1_FILES=(~{sep=' ' r1_fastqs})
    R2_FILES=(~{sep=' ' r2_fastqs})


    # Process each pair of files
    for i in "${!R1_FILES[@]}"; do
      R1_FILE="${R1_FILES[$i]}"
      R2_FILE="${R2_FILES[$i]}"

      # Define Cell Ranger naming convention for this pair
      TARGET_R1="~{sample_name}_S1_L001_R1_001.fastq.gz"
      TARGET_R2="~{sample_name}_S1_L001_R2_001.fastq.gz"

      cp "$R1_FILE" "renamed_fastqs/$TARGET_R1"
      cp "$R2_FILE" "renamed_fastqs/$TARGET_R2"
    done
  >>>

  output {
    Array[File] renamed_fastqs = glob("renamed_fastqs/*.fastq.gz")
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: 1
    memory: "2 GB"
  }
}
