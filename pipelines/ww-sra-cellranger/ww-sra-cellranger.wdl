version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/skip-non-single-cell/modules/ww-sra/ww-sra.wdl" as sra_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/skip-non-single-cell/modules/ww-cellranger/ww-cellranger.wdl" as cellranger_tasks

workflow sra_cellranger {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow to download single-cell RNA-seq data from SRA and process using Cell Ranger count"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-sra-cellranger/ww-sra-cellranger.wdl"
    outputs: {
        single_cell_sample_list: "Sample IDs that ran Cell Ranger successfully (one per line)",
        skipped_sample_list: "Sample IDs skipped due to chemistry-detection failure (one per line)",
        cellranger_results: "Cell Ranger count output tarballs (successful samples only)",
        cellranger_web_summaries: "Cell Ranger web summary HTML files (successful samples only)",
        cellranger_metrics: "Cell Ranger metrics summary CSV files (successful samples only)",
        cellranger_filtered_h5s: "Cell Ranger filtered feature-barcode matrix HDF5 files (successful samples only)"
    }
  }

  parameter_meta {
    sra_id_list: "list of SRA sample ID's to be pulled down and processed with Cell Ranger"
    ref_gex: "GEX reference transcriptome tarball for Cell Ranger"
    ncpu: "number of CPUs to use for SRA download and Cell Ranger processing"
    memory_gb: "memory allocation in GB for Cell Ranger tasks"
    max_reads: "Optional maximum number of reads to download from SRA (for testing/downsampling). If not specified, downloads all reads."
    ngc_file: "Optional NGC repository key file for downloading controlled-access dbGaP data."
    create_bam: "Whether Cell Ranger should generate a BAM file (default: true)"
    expect_cells: "Optional expected number of recovered cells per sample"
    chemistry: "Optional assay configuration for Cell Ranger (e.g. SC3Pv2, SC3Pv3)"
    skip_on_chemistry_failure: "If true, samples Cell Ranger can't assign a chemistry to (non-single-cell input, or too few/too short reads) are skipped (recorded in skipped_sample_list, absent from cellranger_* arrays) instead of failing the workflow. Default: false."
    execution_mode: "Which Cell Ranger task to dispatch to: 'docker' (default; private Cell Ranger image), 'hpc_cromwell' (loads Cell Ranger via host env-modules under Cromwell-on-HPC), or 'hpc_sprocket' (same module-load approach inside a Lua container under Sprocket-on-HPC). Any other value fails loudly via select_first."
    docker_image: "Private Cell Ranger Docker image used by run_count. Ignored unless execution_mode is 'docker'."
    cellranger_module: "HPC environment module used by the run_count_hpc_* tasks. Ignored unless execution_mode starts with 'hpc_'."
  }

  input {
    Array[String] sra_id_list
    File ref_gex
    Int ncpu = 8
    Int memory_gb = 64
    Int? max_reads
    File? ngc_file
    Boolean create_bam = true
    Int? expect_cells
    String? chemistry
    Boolean skip_on_chemistry_failure = false
    String execution_mode = "docker"
    String docker_image = "ghcr.io/getwilds/cellranger:10.0.0"
    String cellranger_module = "CellRanger/10.0.0"
  }

  scatter (id in sra_id_list) {
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
          create_bam = create_bam,
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
          create_bam = create_bam,
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
          create_bam = create_bam,
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
  }

  # Partition sample_ids into "ran" vs "skipped" based on each
  # chemistry_status file's contents.
  call summarize_chemistry_status { input:
      sample_ids = sra_id_list,
      chemistry_status_files = chemistry_status_file
  }

  output {
    File single_cell_sample_list = summarize_chemistry_status.single_cell_sample_list
    File skipped_sample_list = summarize_chemistry_status.skipped_sample_list
    Array[File] cellranger_results = select_all(results_tar)
    Array[File] cellranger_web_summaries = select_all(web_summary)
    Array[File] cellranger_metrics = select_all(metrics_summary)
    Array[File] cellranger_filtered_h5s = select_all(filtered_h5)
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
