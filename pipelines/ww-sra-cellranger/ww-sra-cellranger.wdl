version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl" as sra_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/organize-sra-cellranger-outputs/modules/ww-cellranger/ww-cellranger.wdl" as cellranger_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-cellbender/modules/ww-cellbender/ww-cellbender.wdl" as cellbender_tasks

workflow sra_cellranger {
  meta {
    author: [
        {
            name: "Taylor Firman",
            email: "tfirman@fredhutch.org"
        },
        {
            name: "Hrishi Venkatesh",
            email: "hvenkate@fredhutch.org"
        }
    ]
    description: "WDL workflow to download single-cell RNA-seq data from SRA and process using Cell Ranger count"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-sra-cellranger/ww-sra-cellranger.wdl"
    outputs: {
        single_cell_sample_list: "Sample IDs that ran Cell Ranger successfully (one per line)",
        skipped_sample_list: "Sample IDs skipped due to chemistry-detection failure (one per line)",
        cellranger_results: "Cell Ranger count output tarballs",
        cellranger_web_summaries: "Cell Ranger web summary HTML files",
        cellranger_metrics: "Cell Ranger metrics summary CSV files",
        cellranger_filtered_h5s: "Cell Ranger filtered feature-barcode matrix HDF5 files",
        cellranger_raw_h5s: "Cell Ranger raw feature-barcode matrix HDF5 files",
        cellbender_output_h5s: "CellBender cleaned count matrix H5 files (all barcodes retained), one per successful sample",
        cellbender_filtered_h5s: "CellBender filtered count matrix H5 files (barcodes with >50% cell probability), one per successful sample",
        organized_results: "Tarball of all Cell Ranger outputs organized into per-sample subdirectories (absent when organize_results is false)"
    }
  }

  parameter_meta {
    sra_id_list: "Array of SRA sample ID's to be pulled down and processed with Cell Ranger. Provide this or `sra_id_file` (the file takes precedence when both are given)."
    sra_id_file: "Optional text file of SRA sample ID's, one per line, as an alternative to `sra_id_list`. This is the format produced by the 'Accession List' button in NCBI's SRA Run Selector."
    ref_gex: "GEX reference transcriptome tarball for Cell Ranger"
    ncpu: "number of CPUs to use for SRA download and Cell Ranger processing"
    memory_gb: "memory allocation in GB for Cell Ranger tasks"
    max_reads: "Optional maximum number of reads to download from SRA (for testing/downsampling). If not specified, downloads all reads."
    ngc_file: "Optional NGC repository key file for downloading controlled-access dbGaP data."
    expect_cells: "Optional expected number of recovered cells per sample"
    chemistry: "Optional assay configuration for Cell Ranger (e.g. SC3Pv2, SC3Pv3)"
    skip_on_chemistry_failure: "When `true`, let task succeed with absent outputs if chemistry can't be auto detected."
    execution_mode: "Which Cell Ranger task to dispatch to: 'docker' (default; private Cell Ranger image), 'hpc_cromwell' (loads Cell Ranger via host env-modules under Cromwell-on-HPC), or 'hpc_sprocket' (same module-load approach inside a Lua container under Sprocket-on-HPC). Any other value fails loudly via select_first."
    docker_image: "Private Cell Ranger Docker image used by run_count. Ignored unless execution_mode is 'docker'."
    cellranger_module: "HPC environment module used by the run_count_hpc_* tasks. Ignored unless execution_mode starts with 'hpc_'."
    organize_results: "When true, package all Cell Ranger outputs into a tarball organized by sample subdirectory."
    output_prefix: "Prefix for the organized results tarball filename (default: 'cellranger_results')."
    cellbender_gpu_enabled: "Enable GPU acceleration for CellBender (default: true); set to false for CPU-only execution."
    cellbender_expected_cells: "Optional expected number of real cells per sample passed to CellBender."
    cellbender_total_droplets_included: "Optional total number of droplets for CellBender to analyze per sample."
    cellbender_epochs: "Number of CellBender training epochs (default: 150)."
    cellbender_cpu_cores: "Number of CPU cores for CellBender (default: 4)."
    cellbender_memory_gb: "Memory in GB for CellBender (default: 32)."
  }

  input {
    Array[String]? sra_id_list
    File? sra_id_file
    File ref_gex
    Int ncpu = 8
    Int memory_gb = 64
    Int? max_reads
    File? ngc_file
    Int? expect_cells
    String? chemistry
    Boolean skip_on_chemistry_failure = false
    String execution_mode = "docker"
    String docker_image = "ghcr.io/getwilds/cellranger:10.0.0"
    String cellranger_module = "CellRanger/10.0.0"
    Boolean organize_results = false
    String output_prefix = "cellranger_results"
    Boolean cellbender_gpu_enabled = true
    Int? cellbender_expected_cells
    Int? cellbender_total_droplets_included
    Int cellbender_epochs = 150
    Int cellbender_cpu_cores = 4
    Int cellbender_memory_gb = 32
  }

  # Resolve the sample list from whichever input was provided. A
  # Run-Selector-style text file (one accession per line) takes
  # precedence over the inline array; supplying neither fails loudly
  # via select_first.
  Array[String] sra_ids = if defined(sra_id_file)
    then read_lines(select_first([sra_id_file]))
    else select_first([sra_id_list])

  scatter (id in sra_ids) {
    # Download FASTQ files from SRA
    call sra_tasks.fastqdump { input:
        sra_id = id,
        ncpu = ncpu,
        max_reads = max_reads,
        ngc_file = ngc_file
    }

    # Rename FASTQs to Cell Ranger naming convention
    call cellranger_tasks.rename_fastqs { input:
        r1_fastq = fastqdump.r1_end,
        r2_fastq = fastqdump.r2_end,
        sample_id = id
    }

    # Dispatch to the Docker, HPC-Cromwell, or HPC-Sprocket variant of
    # cellranger count. Each branch's outputs are Optional; select_first
    # below collapses them to a single Optional File per iteration (the
    # outer Optional comes from skip_on_chemistry_failure skipping a
    # sample inside the chosen branch). An invalid execution_mode causes
    # all branches to skip and the select_first call on chemistry_status
    # to fail loudly.
    if (execution_mode == "docker") {
      call cellranger_tasks.run_count { input:
          r1_fastqs = [rename_fastqs.r1_renamed],
          r2_fastqs = [rename_fastqs.r2_renamed],
          ref_gex = ref_gex,
          sample_id = id,
          create_bam = false,
          cpu_cores = ncpu,
          memory_gb = memory_gb,
          expect_cells = expect_cells,
          chemistry = chemistry,
          skip_on_chemistry_failure = skip_on_chemistry_failure,
          docker_image = docker_image
      }
    }

    if (execution_mode == "hpc_cromwell") {
      call cellranger_tasks.run_count_hpc_cromwell { input:
          r1_fastqs = [rename_fastqs.r1_renamed],
          r2_fastqs = [rename_fastqs.r2_renamed],
          ref_gex = ref_gex,
          sample_id = id,
          create_bam = false,
          cpu_cores = ncpu,
          memory_gb = memory_gb,
          expect_cells = expect_cells,
          chemistry = chemistry,
          skip_on_chemistry_failure = skip_on_chemistry_failure,
          cellranger_module = cellranger_module
      }
    }

    if (execution_mode == "hpc_sprocket") {
      call cellranger_tasks.run_count_hpc_sprocket { input:
          r1_fastqs = [rename_fastqs.r1_renamed],
          r2_fastqs = [rename_fastqs.r2_renamed],
          ref_gex = ref_gex,
          sample_id = id,
          create_bam = false,
          cpu_cores = ncpu,
          memory_gb = memory_gb,
          expect_cells = expect_cells,
          chemistry = chemistry,
          skip_on_chemistry_failure = skip_on_chemistry_failure,
          cellranger_module = cellranger_module
      }
    }

    # chemistry_status is non-optional within each branch, so the
    # outer select_first only needs to pick which branch ran. If
    # execution_mode is invalid, this fails loudly.
    File chemistry_status_file = select_first([run_count.chemistry_status, run_count_hpc_cromwell.chemistry_status, run_count_hpc_sprocket.chemistry_status])

    # Count outputs are File? within each branch (absent on skip).
    # We can't use select_first here: Cromwell errors when every
    # element is None even when the result type is File?. Pick the
    # value from whichever branch ran, defaulting to None.
    File? results_tar = if defined(run_count.results_tar) then run_count.results_tar
                       else if defined(run_count_hpc_cromwell.results_tar) then run_count_hpc_cromwell.results_tar
                       else run_count_hpc_sprocket.results_tar
    File? web_summary = if defined(run_count.web_summary) then run_count.web_summary
                       else if defined(run_count_hpc_cromwell.web_summary) then run_count_hpc_cromwell.web_summary
                       else run_count_hpc_sprocket.web_summary
    File? metrics_summary = if defined(run_count.metrics_summary) then run_count.metrics_summary
                           else if defined(run_count_hpc_cromwell.metrics_summary) then run_count_hpc_cromwell.metrics_summary
                           else run_count_hpc_sprocket.metrics_summary
    File? filtered_h5 = if defined(run_count.filtered_h5) then run_count.filtered_h5
                       else if defined(run_count_hpc_cromwell.filtered_h5) then run_count_hpc_cromwell.filtered_h5
                       else run_count_hpc_sprocket.filtered_h5
    File? raw_h5 = if defined(run_count.raw_h5) then run_count.raw_h5
                  else if defined(run_count_hpc_cromwell.raw_h5) then run_count_hpc_cromwell.raw_h5
                  else run_count_hpc_sprocket.raw_h5

    # Run CellBender on samples that produced a raw h5 (skipped samples won't have one)
    if (defined(raw_h5)) {
      call cellbender_tasks.remove_background { input:
        input_h5 = select_first([raw_h5]),
        sample_name = id,
        expected_cells = cellbender_expected_cells,
        total_droplets_included = cellbender_total_droplets_included,
        epochs = cellbender_epochs,
        gpu_enabled = cellbender_gpu_enabled,
        cpu_cores = cellbender_cpu_cores,
        memory_gb = cellbender_memory_gb
      }
    }
  }

  # Partition sample_ids into "ran" vs "skipped" based on each
  # chemistry_status file's contents.
  call summarize_chemistry_status { input:
      sample_ids = sra_ids,
      chemistry_status_files = chemistry_status_file
  }

  # Optional: package all outputs into per-sample subdirectories. Uses
  # read_lines on the single_cell_sample_list so the sample IDs are in
  # the same order as the select_all'd file arrays (skipped samples are
  # absent from both).
  if (organize_results) {
    call organize_outputs { input:
        sample_ids = read_lines(summarize_chemistry_status.single_cell_sample_list),
        results_tars = select_all(results_tar),
        web_summaries = select_all(web_summary),
        metrics_summaries = select_all(metrics_summary),
        filtered_h5s = select_all(filtered_h5),
        raw_h5s = select_all(raw_h5),
        cellbender_output_h5s = select_all(remove_background.output_h5),
        cellbender_filtered_h5s = select_all(remove_background.filtered_h5),
        output_prefix = output_prefix
    }
  }

  output {
    File single_cell_sample_list = summarize_chemistry_status.single_cell_sample_list
    File skipped_sample_list = summarize_chemistry_status.skipped_sample_list
    Array[File] cellranger_results = select_all(results_tar)
    Array[File] cellranger_web_summaries = select_all(web_summary)
    Array[File] cellranger_metrics = select_all(metrics_summary)
    Array[File] cellranger_filtered_h5s = select_all(filtered_h5)
    Array[File] cellranger_raw_h5s = select_all(raw_h5)
    Array[File] cellbender_output_h5s = select_all(remove_background.output_h5)
    Array[File] cellbender_filtered_h5s = select_all(remove_background.filtered_h5)
    File? organized_results = organize_outputs.results_zip
  }
}

task summarize_chemistry_status {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Partition sample IDs into ran-successfully vs skipped lists based on per-sample chemistry_status marker files from ww-cellranger."
    outputs: {
        single_cell_sample_list: "Newline-delimited sample IDs whose chemistry_status was 'ok'",
        skipped_sample_list: "Newline-delimited sample IDs whose chemistry_status was 'skipped_non_single_cell'"
    }
  }

  parameter_meta {
    sample_ids: "Sample IDs in the same order as chemistry_status_files"
    chemistry_status_files: "Per-sample chemistry_status marker files from ww-cellranger"
  }

  input {
    Array[String] sample_ids
    Array[File] chemistry_status_files
  }

  command <<<
    set -eo pipefail

    IDS=(~{sep=' ' sample_ids})
    FILES=(~{sep=' ' chemistry_status_files})

    if [ "${#IDS[@]}" -ne "${#FILES[@]}" ]; then
      echo "ERROR: sample_ids (${#IDS[@]}) and chemistry_status_files (${#FILES[@]}) have different lengths" >&2
      exit 1
    fi

    : > single_cell_sample_list.txt
    : > skipped_sample_list.txt

    for i in "${!IDS[@]}"; do
      status=$(tr -d '[:space:]' < "${FILES[$i]}")
      case "$status" in
        ok)
          echo "${IDS[$i]}" >> single_cell_sample_list.txt
          ;;
        skipped_non_single_cell)
          echo "${IDS[$i]}" >> skipped_sample_list.txt
          ;;
        *)
          echo "ERROR: unexpected chemistry_status '$status' for sample '${IDS[$i]}'" >&2
          exit 1
          ;;
      esac
    done
  >>>

  output {
    File single_cell_sample_list = "single_cell_sample_list.txt"
    File skipped_sample_list = "skipped_sample_list.txt"
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: 1
    memory: "2 GB"
  }
}

task organize_outputs {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Package Cell Ranger outputs into a tarball organized by sample subdirectory."
    outputs: {
        results_zip: "Tarball containing all Cell Ranger and CellBender outputs organized into per-sample subdirectories"
    }
  }

  parameter_meta {
    sample_ids: "Sample IDs in the same order as the file arrays (skipped samples must be excluded)"
    results_tars: "Cell Ranger output tarballs, one per successful sample"
    web_summaries: "Web summary HTML files, one per successful sample"
    metrics_summaries: "Metrics summary CSV files, one per successful sample"
    filtered_h5s: "Filtered feature-barcode matrix HDF5 files, one per successful sample"
    raw_h5s: "Raw feature-barcode matrix HDF5 files, one per successful sample"
    cellbender_output_h5s: "CellBender cleaned count matrix H5 files, one per successful sample"
    cellbender_filtered_h5s: "CellBender filtered count matrix H5 files, one per successful sample"
    output_prefix: "Prefix for the output tarball filename"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
  }

  input {
    Array[String] sample_ids
    Array[File] results_tars
    Array[File] web_summaries
    Array[File] metrics_summaries
    Array[File] filtered_h5s
    Array[File] raw_h5s
    Array[File] cellbender_output_h5s
    Array[File] cellbender_filtered_h5s
    String output_prefix = "cellranger_results"
    Int memory_gb = 4
    Int cpu_cores = 1
  }

  command <<<
    set -eo pipefail

    OUTDIR="~{output_prefix}"
    SAMPLE_IDS=(~{sep=' ' sample_ids})
    TARS=(~{sep=' ' results_tars})
    WEB=(~{sep=' ' web_summaries})
    METRICS=(~{sep=' ' metrics_summaries})
    FILTERED=(~{sep=' ' filtered_h5s})
    RAW=(~{sep=' ' raw_h5s})
    CB_OUTPUT=(~{sep=' ' cellbender_output_h5s})
    CB_FILTERED=(~{sep=' ' cellbender_filtered_h5s})

    for i in "${!SAMPLE_IDS[@]}"; do
      SAMPLE="${SAMPLE_IDS[$i]}"
      mkdir -p "$OUTDIR/$SAMPLE"
      cp "${TARS[$i]}"        "$OUTDIR/$SAMPLE/"
      cp "${WEB[$i]}"         "$OUTDIR/$SAMPLE/"
      cp "${METRICS[$i]}"     "$OUTDIR/$SAMPLE/"
      cp "${FILTERED[$i]}"    "$OUTDIR/$SAMPLE/"
      cp "${RAW[$i]}"         "$OUTDIR/$SAMPLE/"
      cp "${CB_OUTPUT[$i]}"   "$OUTDIR/$SAMPLE/"
      cp "${CB_FILTERED[$i]}" "$OUTDIR/$SAMPLE/"
    done

    tar -czf "~{output_prefix}.tar.gz" "$OUTDIR"
  >>>

  output {
    File results_zip = "~{output_prefix}.tar.gz"
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
