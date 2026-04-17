## WILDS WDL for running Cell Ranger pipelines.
## Intended for use alone or as a modular component in the WILDS ecosystem.

version 1.0

task run_count {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Run cellranger count on gene expression reads from one GEM well using a private Cell Ranger Docker image. Cell Ranger is not redistributable, so no public WILDS image is provided — see https://github.com/getwilds/wilds-docker-library/blob/main/cellranger/Dockerfile_latest for a recipe to build your own."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellranger/ww-cellranger.wdl"
    outputs: {
        results_tar: "Compressed tarball of Cell Ranger count output directory",
        web_summary: "Web summary HTML file",
        metrics_summary: "Metrics summary CSV file"
    }
    topic: "transcriptomics,gene_expression"
    species: "any"
    operation: "rna_seq_quantification"
    in_sample_req: "r1_fastqs:rna_sequence:fastq,r2_fastqs:rna_sequence:fastq"
    in_sample_opt: "none"
    in_ref_req: "ref_gex:data_index:tar_format"
    in_ref_opt: "none"
    out_sample: "results_tar:gene_expression_matrix:tar_format,web_summary:quality_control_report:html,metrics_summary:quality_control_report:csv"
    out_ref: "none"
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
    docker_image: "Private Cell Ranger Docker image. No public image is provided because Cell Ranger is not redistributable; build your own from the WILDS Dockerfile recipe."
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
    String docker_image = "ghcr.io/getwilds/cellranger:10.0.0"
  }

  # Keep command block in sync with run_count_hpc.
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

    # Clean up Cell Ranger output directory for housekeeping
    # (and to avoid Sprocket symlink validation errors)
    rm -rf "~{sample_id}"
  >>>

  output {
    File results_tar = "~{sample_id}_outs.tar.gz"
    File web_summary = "web_summary.html"
    File metrics_summary = "metrics_summary.csv"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task run_count_hpc {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Run cellranger count using an HPC environment module instead of a Docker image. Intended for institutional HPC backends (e.g. Fred Hutch Cromwell) where Cell Ranger is provided via Lmod/environment modules under an institutional license."
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
    environment_modules: "Space-separated list of HPC environment modules to load (e.g. 'CellRanger/10.0.0'). Honored only by backends configured with a 'modules' runtime attribute."
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
    String environment_modules = "CellRanger/10.0.0"
  }

  # Keep command block in sync with run_count.
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

    # Clean up Cell Ranger output directory for housekeeping
    # (and to avoid Sprocket symlink validation errors)
    rm -rf "~{sample_id}"
  >>>

  output {
    File results_tar = "~{sample_id}_outs.tar.gz"
    File web_summary = "web_summary.html"
    File metrics_summary = "metrics_summary.csv"
  }

  runtime {
    modules: environment_modules
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task rename_fastqs {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Renames FASTQ files to match the Cell Ranger naming convention (SampleName_S1_R1_001.fastq.gz). Useful for preparing SRA downloads or other non-standard FASTQ files for Cell Ranger input."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellranger/ww-cellranger.wdl"
    outputs: {
        r1_renamed: "R1 FASTQ file renamed to Cell Ranger naming convention",
        r2_renamed: "R2 FASTQ file renamed to Cell Ranger naming convention"
    }
  }

  parameter_meta {
    r1_fastq: "R1 FASTQ file to rename"
    r2_fastq: "R2 FASTQ file to rename"
    sample_id: "Sample ID to use as the prefix in the renamed file"
  }

  input {
    File r1_fastq
    File r2_fastq
    String sample_id
  }

  command <<<
    set -eo pipefail
    cp "~{r1_fastq}" "~{sample_id}_S1_R1_001.fastq.gz"
    cp "~{r2_fastq}" "~{sample_id}_S1_R2_001.fastq.gz"
  >>>

  output {
    File r1_renamed = "~{sample_id}_S1_R1_001.fastq.gz"
    File r2_renamed = "~{sample_id}_S1_R2_001.fastq.gz"
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: 1
    memory: "2 GB"
  }
}
