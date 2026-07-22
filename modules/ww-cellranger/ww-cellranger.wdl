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
        chemistry_status: "Marker file: 'ok' or 'skipped_non_single_cell'",
        results_tar: "Compressed tarball of Cell Ranger count output directory (absent if skipped)",
        web_summary: "Web summary HTML file (absent if skipped)",
        metrics_summary: "Metrics summary CSV file (absent if skipped)",
        filtered_h5: "Filtered feature-barcode matrix HDF5 file (absent if skipped)",
        raw_h5: "Raw feature-barcode matrix HDF5 file (absent if skipped)"
    }
    topic: "transcriptomics,gene_expression,single_cell_sequencing"
    species: "human,eukaryote,prokaryote,virus"
    operation: "rna_seq_quantification"
    input_sample_required: "r1_fastqs:rna_sequence:fastq,r2_fastqs:rna_sequence:fastq"
    input_sample_optional: "none"
    input_reference_required: "ref_gex:data_index:tar_format"
    input_reference_optional: "none"
    output_sample: "results_tar:gene_expression_matrix:tar_format,web_summary:quality_control_report:html,metrics_summary:quality_control_report:csv,filtered_h5:gene_expression_matrix:h5,raw_h5:gene_expression_matrix:h5"
    output_reference: "none"
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
    skip_on_chemistry_failure: "When `true`, let task succeed with absent outputs if chemistry can't be auto detected."
    docker_image: "Private Cell Ranger Docker image. No public image is provided because Cell Ranger is not redistributable; build your own from the WILDS Dockerfile recipe. Compatible with Cell Ranger 8.x-10.x; v7 and earlier are not supported."
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
    Boolean skip_on_chemistry_failure = false
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

    # Run cellranger count, capturing both stdout and stderr so we can
    # detect a chemistry-detection failure (typically non-single-cell
    # input) and skip gracefully when skip_on_chemistry_failure=true.
    # Cell Ranger writes most user-facing error text (including the
    # TXRNGR* codes) to stdout, not stderr, so capturing both is
    # required.
    set +e
    cellranger count \
      --transcriptome=gex_ref \
      --fastqs="$FASTQS" \
      --localcores=~{cpu_cores} \
      --localmem=~{memory_gb} \
      --output-dir="~{sample_id}" \
      ~{"--expect-cells=" + expect_cells} \
      ~{"--chemistry=" + chemistry} \
      --id="~{sample_id}" \
      --create-bam="~{create_bam}" > cellranger.log 2>&1
    CR_EXIT=$?
    set -e
    cat cellranger.log

    if [ "$CR_EXIT" -ne 0 ]; then
      # Heuristic chemistry-detection-failure markers (may need tuning
      # across Cell Ranger versions): can't pick a chemistry, read
      # lengths too short, or barcodes don't match the 10x whitelist
      # (all symptoms of non-single-cell input or unusable single-cell
      # input).
      if ~{true="true" false="false" skip_on_chemistry_failure} \
        && grep -qE -i 'could not (auto)?detect|ambiguous chemistry|chemistry .* could not be|NO_INPUT_ANTIBODY_READS|TXRNGR10002|TXRNGR10004|read lengths are incompatible|low rate of correct barcodes' cellranger.log; then
        echo "Cell Ranger could not assign a chemistry to ~{sample_id} (non-single-cell input, or too few/too short reads); skipping (skip_on_chemistry_failure=true)." >&2
        echo "skipped_non_single_cell" > chemistry_status.txt
        rm -rf "~{sample_id}"
        exit 0
      fi
      # Real failure: write a marker so the chemistry_status output is
      # always present (Cromwell's output-collection phase requires it),
      # then re-raise the original exit code.
      echo "failed" > chemistry_status.txt
      exit "$CR_EXIT"
    fi

    echo "ok" > chemistry_status.txt

    # Create tarball of output directory
    tar -czf "~{sample_id}_outs.tar.gz" "~{sample_id}/outs"

    # Move output files to working directory for outputting, prefixing
    # each with the sample ID so basenames stay unique when a downstream
    # workflow flattens multiple samples into one directory (e.g.
    # Cromwell's final_workflow_outputs_dir).
    mv "~{sample_id}/outs/web_summary.html" "~{sample_id}_web_summary.html"
    mv "~{sample_id}/outs/metrics_summary.csv" "~{sample_id}_metrics_summary.csv"
    mv "~{sample_id}/outs/filtered_feature_bc_matrix.h5" "~{sample_id}_filtered_feature_bc_matrix.h5"
    mv "~{sample_id}/outs/raw_feature_bc_matrix.h5" "~{sample_id}_raw_feature_bc_matrix.h5"

    # Clean up Cell Ranger output directory for housekeeping
    # (and to avoid Sprocket symlink validation errors)
    rm -rf "~{sample_id}"
  >>>

  output {
    File chemistry_status = "chemistry_status.txt"
    File? results_tar = "~{sample_id}_outs.tar.gz"
    File? web_summary = "~{sample_id}_web_summary.html"
    File? metrics_summary = "~{sample_id}_metrics_summary.csv"
    File? filtered_h5 = "~{sample_id}_filtered_feature_bc_matrix.h5"
    File? raw_h5 = "~{sample_id}_raw_feature_bc_matrix.h5"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task run_count_hpc_cromwell {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Cromwell-on-HPC variant of cellranger count: omits the docker runtime key so Cromwell runs the task directly on the compute node, where the host's environment-module system can load Cell Ranger."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellranger/ww-cellranger.wdl"
    outputs: {
        chemistry_status: "Marker file: 'ok' or 'skipped_non_single_cell'",
        results_tar: "Compressed tarball of Cell Ranger count output directory (absent if skipped)",
        web_summary: "Web summary HTML file (absent if skipped)",
        metrics_summary: "Metrics summary CSV file (absent if skipped)",
        filtered_h5: "Filtered feature-barcode matrix HDF5 file (absent if skipped)",
        raw_h5: "Raw feature-barcode matrix HDF5 file (absent if skipped)"
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
    skip_on_chemistry_failure: "When `true`, let task succeed with absent outputs if chemistry can't be auto detected."
    cellranger_module: "HPC environment module to load for Cell Ranger (e.g. 'CellRanger/10.0.0'). Compatible with Cell Ranger 8.x-10.x; v7 and earlier are not supported."
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
    Boolean skip_on_chemistry_failure = false
    String cellranger_module = "CellRanger/10.0.0"
  }

  # Keep command block in sync with run_count and run_count_hpc_sprocket.
  command <<<
    set -eo pipefail

    # Load Cell Ranger from the host's environment-module system.
    . /app/lmod/lmod/init/bash
    module load ~{cellranger_module}

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

    # Run cellranger count, capturing both stdout and stderr so we can
    # detect a chemistry-detection failure (typically non-single-cell
    # input) and skip gracefully when skip_on_chemistry_failure=true.
    # Cell Ranger writes most user-facing error text (including the
    # TXRNGR* codes) to stdout, not stderr, so capturing both is
    # required.
    set +e
    cellranger count \
      --transcriptome=gex_ref \
      --fastqs="$FASTQS" \
      --localcores=~{cpu_cores} \
      --localmem=~{memory_gb} \
      --output-dir="~{sample_id}" \
      ~{"--expect-cells=" + expect_cells} \
      ~{"--chemistry=" + chemistry} \
      --id="~{sample_id}" \
      --create-bam="~{create_bam}" > cellranger.log 2>&1
    CR_EXIT=$?
    set -e
    cat cellranger.log

    if [ "$CR_EXIT" -ne 0 ]; then
      # Heuristic chemistry-detection-failure markers (may need tuning
      # across Cell Ranger versions): can't pick a chemistry, read
      # lengths too short, or barcodes don't match the 10x whitelist
      # (all symptoms of non-single-cell input or unusable single-cell
      # input).
      if ~{true="true" false="false" skip_on_chemistry_failure} \
        && grep -qE -i 'could not (auto)?detect|ambiguous chemistry|chemistry .* could not be|NO_INPUT_ANTIBODY_READS|TXRNGR10002|TXRNGR10004|read lengths are incompatible|low rate of correct barcodes' cellranger.log; then
        echo "Cell Ranger could not assign a chemistry to ~{sample_id} (non-single-cell input, or too few/too short reads); skipping (skip_on_chemistry_failure=true)." >&2
        echo "skipped_non_single_cell" > chemistry_status.txt
        rm -rf "~{sample_id}"
        exit 0
      fi
      # Real failure: write a marker so the chemistry_status output is
      # always present (Cromwell's output-collection phase requires it),
      # then re-raise the original exit code.
      echo "failed" > chemistry_status.txt
      exit "$CR_EXIT"
    fi

    echo "ok" > chemistry_status.txt

    # Create tarball of output directory
    tar -czf "~{sample_id}_outs.tar.gz" "~{sample_id}/outs"

    # Move output files to working directory for outputting, prefixing
    # each with the sample ID so basenames stay unique when a downstream
    # workflow flattens multiple samples into one directory (e.g.
    # Cromwell's final_workflow_outputs_dir).
    mv "~{sample_id}/outs/web_summary.html" "~{sample_id}_web_summary.html"
    mv "~{sample_id}/outs/metrics_summary.csv" "~{sample_id}_metrics_summary.csv"
    mv "~{sample_id}/outs/filtered_feature_bc_matrix.h5" "~{sample_id}_filtered_feature_bc_matrix.h5"
    mv "~{sample_id}/outs/raw_feature_bc_matrix.h5" "~{sample_id}_raw_feature_bc_matrix.h5"

    # Clean up Cell Ranger output directory for housekeeping
    # (and to avoid Sprocket symlink validation errors)
    rm -rf "~{sample_id}"
  >>>

  output {
    File chemistry_status = "chemistry_status.txt"
    File? results_tar = "~{sample_id}_outs.tar.gz"
    File? web_summary = "~{sample_id}_web_summary.html"
    File? metrics_summary = "~{sample_id}_metrics_summary.csv"
    File? filtered_h5 = "~{sample_id}_filtered_feature_bc_matrix.h5"
    File? raw_h5 = "~{sample_id}_raw_feature_bc_matrix.h5"
  }

  # No docker/container runtime key: Cromwell's HPC backend runs this
  # task directly on the compute node, where `module load` makes the
  # licensed Cell Ranger binary available on PATH.
  runtime {
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task run_count_hpc_sprocket {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Sprocket-on-HPC variant of cellranger count: runs inside a minimal Lua container so the host's bind-mounted Lmod can execute, while the Cell Ranger binary itself is bind-mounted in from the host."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellranger/ww-cellranger.wdl"
    outputs: {
        chemistry_status: "Marker file: 'ok' or 'skipped_non_single_cell'",
        results_tar: "Compressed tarball of Cell Ranger count output directory (absent if skipped)",
        web_summary: "Web summary HTML file (absent if skipped)",
        metrics_summary: "Metrics summary CSV file (absent if skipped)",
        filtered_h5: "Filtered feature-barcode matrix HDF5 file (absent if skipped)",
        raw_h5: "Raw feature-barcode matrix HDF5 file (absent if skipped)"
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
    skip_on_chemistry_failure: "When `true`, let task succeed with absent outputs if chemistry can't be auto detected."
    cellranger_module: "HPC environment module to load for Cell Ranger (e.g. 'CellRanger/10.0.0'). Compatible with Cell Ranger 8.x-10.x; v7 and earlier are not supported."
    docker_image: "Docker image to use for this task"
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
    Boolean skip_on_chemistry_failure = false
    String cellranger_module = "CellRanger/10.0.0"
    String docker_image = "getwilds/lua:5.3.6"
  }

  # Keep command block in sync with run_count and run_count_hpc_cromwell.
  command <<<
    set -eo pipefail

    # Load Cell Ranger from the host's environment-module system
    # (Lmod tree, modulefiles, and software tree are bind-mounted in
    # via the Sprocket HPC config; see .github/configs/sprocket-hpc.toml).
    . /app/lmod/lmod/init/bash
    module load ~{cellranger_module}

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

    # Run cellranger count, capturing both stdout and stderr so we can
    # detect a chemistry-detection failure (typically non-single-cell
    # input) and skip gracefully when skip_on_chemistry_failure=true.
    # Cell Ranger writes most user-facing error text (including the
    # TXRNGR* codes) to stdout, not stderr, so capturing both is
    # required.
    set +e
    cellranger count \
      --transcriptome=gex_ref \
      --fastqs="$FASTQS" \
      --localcores=~{cpu_cores} \
      --localmem=~{memory_gb} \
      --output-dir="~{sample_id}" \
      ~{"--expect-cells=" + expect_cells} \
      ~{"--chemistry=" + chemistry} \
      --id="~{sample_id}" \
      --create-bam="~{create_bam}" > cellranger.log 2>&1
    CR_EXIT=$?
    set -e
    cat cellranger.log

    if [ "$CR_EXIT" -ne 0 ]; then
      # Heuristic chemistry-detection-failure markers (may need tuning
      # across Cell Ranger versions): can't pick a chemistry, read
      # lengths too short, or barcodes don't match the 10x whitelist
      # (all symptoms of non-single-cell input or unusable single-cell
      # input).
      if ~{true="true" false="false" skip_on_chemistry_failure} \
        && grep -qE -i 'could not (auto)?detect|ambiguous chemistry|chemistry .* could not be|NO_INPUT_ANTIBODY_READS|TXRNGR10002|TXRNGR10004|read lengths are incompatible|low rate of correct barcodes' cellranger.log; then
        echo "Cell Ranger could not assign a chemistry to ~{sample_id} (non-single-cell input, or too few/too short reads); skipping (skip_on_chemistry_failure=true)." >&2
        echo "skipped_non_single_cell" > chemistry_status.txt
        rm -rf "~{sample_id}"
        exit 0
      fi
      # Real failure: write a marker so the chemistry_status output is
      # always present (Cromwell's output-collection phase requires it),
      # then re-raise the original exit code.
      echo "failed" > chemistry_status.txt
      exit "$CR_EXIT"
    fi

    echo "ok" > chemistry_status.txt

    # Create tarball of output directory
    tar -czf "~{sample_id}_outs.tar.gz" "~{sample_id}/outs"

    # Move output files to working directory for outputting, prefixing
    # each with the sample ID so basenames stay unique when a downstream
    # workflow flattens multiple samples into one directory (e.g.
    # Cromwell's final_workflow_outputs_dir).
    mv "~{sample_id}/outs/web_summary.html" "~{sample_id}_web_summary.html"
    mv "~{sample_id}/outs/metrics_summary.csv" "~{sample_id}_metrics_summary.csv"
    mv "~{sample_id}/outs/filtered_feature_bc_matrix.h5" "~{sample_id}_filtered_feature_bc_matrix.h5"
    mv "~{sample_id}/outs/raw_feature_bc_matrix.h5" "~{sample_id}_raw_feature_bc_matrix.h5"

    # Clean up Cell Ranger output directory for housekeeping
    # (and to avoid Sprocket symlink validation errors)
    rm -rf "~{sample_id}"
  >>>

  output {
    File chemistry_status = "chemistry_status.txt"
    File? results_tar = "~{sample_id}_outs.tar.gz"
    File? web_summary = "~{sample_id}_web_summary.html"
    File? metrics_summary = "~{sample_id}_metrics_summary.csv"
    File? filtered_h5 = "~{sample_id}_filtered_feature_bc_matrix.h5"
    File? raw_h5 = "~{sample_id}_raw_feature_bc_matrix.h5"
  }

  # Minimal Lua-having container so the host's bind-mounted Lmod can
  # execute under Apptainer. Cell Ranger itself is not in this image;
  # it comes in via the host bind-mounts in sprocket-hpc.toml.
  runtime {
    docker: docker_image
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
    docker_image: "Docker image to use for this task"
  }

  input {
    File r1_fastq
    File r2_fastq
    String sample_id
    String docker_image = "ubuntu:22.04"
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
    docker: docker_image
    cpu: 1
    memory: "2 GB"
  }
}

