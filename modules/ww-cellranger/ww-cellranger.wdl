## WILDS WDL for running Cell Ranger pipelines.
## Intended for use alone or as a modular component in the WILDS ecosystem.

version 1.0

task run_count {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Run cellranger count on gene expression reads from one GEM well"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellranger/ww-cellranger.wdl"
    outputs: {
        results_tar: "Compressed tarball of Cell Ranger count output directory",
        web_summary: "Web summary HTML file",
        metrics_summary: "Metrics summary CSV file"
    }
  }

  parameter_meta {
    r1_fastqs: "Array of R1 FASTQ files (contain cell barcodes and UMIs)"
    r2_fastqs: "Array of R2 FASTQ files (contain cDNA sequences)"
    ref_gex: "GEX reference transcriptome tarball"
    sample_id: "Sample ID for output naming"
    create_bam: "Generate BAM file (default: true)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
    expect_cells: "Optional: Expected number of recovered cells"
    chemistry: "Optional: Assay configuration (e.g. SC3Pv2)"
  }

  input {
    Array[File] r1_fastqs
    Array[File] r2_fastqs
    File ref_gex
    String sample_id
    Boolean create_bam = true
    Int cpu_cores = 8
    Int memory_gb = 64
    Int? expect_cells
    String? chemistry
  }

  command <<<
    set -eo pipefail

    # Create arrays from WDL inputs
    R1_FILES=(~{sep=' ' r1_fastqs})
    R2_FILES=(~{sep=' ' r2_fastqs})

    # Validate that R1 and R2 arrays have the same length
    if [ "${#R1_FILES[@]}" -ne "${#R2_FILES[@]}" ]; then
      echo "ERROR: Number of R1 files (${#R1_FILES[@]}) does not match number of R2 files (${#R2_FILES[@]})" >&2
      exit 1
    fi

    # Validate FASTQ naming convention
    # Pattern: SampleName_S[0-9]+_L[0-9]+_R[12]_001.fastq.gz (with lane)
    # Or: SampleName_S[0-9]+_R[12]_001.fastq.gz (without lane, Cell Ranger v4.0+)
    PATTERN_WITH_LANE='_S[0-9]+_L[0-9]+_R[12]_001\.fastq\.gz$'
    PATTERN_WITHOUT_LANE='_S[0-9]+_R[12]_001\.fastq\.gz$'

    validate_fastq_name() {
      local file="$1"
      local file_basename
      file_basename=$(basename "$file")
      if [[ ! "$file_basename" =~ $PATTERN_WITH_LANE ]] && [[ ! "$file_basename" =~ $PATTERN_WITHOUT_LANE ]]; then
        echo "ERROR: FASTQ file '$file_basename' does not follow Cell Ranger naming convention." >&2
        echo "Expected format: SampleName_S1_L001_R1_001.fastq.gz or SampleName_S1_R1_001.fastq.gz" >&2
        exit 1
      fi
    }

    for file in "${R1_FILES[@]}"; do
      validate_fastq_name "$file"
    done
    for file in "${R2_FILES[@]}"; do
      validate_fastq_name "$file"
    done

    echo "FASTQ validation passed"

    # Extract GEX reference
    mkdir -p gex_ref
    tar xf "~{ref_gex}" -C gex_ref --strip-components 1

    # Create folder and copy FASTQ files
    mkdir -p gex_fastqs
    cp ~{sep=' ' r1_fastqs} gex_fastqs/
    cp ~{sep=' ' r2_fastqs} gex_fastqs/
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
      ~{"--chemistry=" + chemistry} \
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
