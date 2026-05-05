version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/split-cicd-hpc-testruns/modules/ww-sra/ww-sra.wdl" as sra_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/split-cicd-hpc-testruns/modules/ww-cellranger/ww-cellranger.wdl" as cellranger_tasks

workflow sra_cellranger {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow to download single-cell RNA-seq data from SRA and process using Cell Ranger count"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-sra-cellranger/ww-sra-cellranger.wdl"
    outputs: {
        cellranger_results: "array of compressed tarballs containing Cell Ranger count output directories",
        cellranger_web_summaries: "array of web summary HTML files for each sample",
        cellranger_metrics: "array of metrics summary CSV files for each sample"
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
    execution_mode: "Which Cell Ranger task to dispatch to. One of: 'docker' (default; uses run_count with a private Cell Ranger Docker image), 'hpc_cromwell' (uses run_count_hpc_cromwell, loads Cell Ranger via the host's environment-module system; intended for Cromwell-on-HPC), or 'hpc_sprocket' (uses run_count_hpc_sprocket, same module-load approach inside a Lua container; intended for Sprocket-on-HPC). Any other value will cause select_first to fail at the end of the scatter."
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
    # below collapses them to a single non-optional File per iteration.
    # An invalid execution_mode causes all branches to skip and the
    # select_first call to fail loudly.
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
          cellranger_module = cellranger_module
      }
    }

    File results_tar = select_first([run_count.results_tar, run_count_hpc_cromwell.results_tar, run_count_hpc_sprocket.results_tar])
    File web_summary = select_first([run_count.web_summary, run_count_hpc_cromwell.web_summary, run_count_hpc_sprocket.web_summary])
    File metrics_summary = select_first([run_count.metrics_summary, run_count_hpc_cromwell.metrics_summary, run_count_hpc_sprocket.metrics_summary])
  }

  output {
    Array[File] cellranger_results = results_tar
    Array[File] cellranger_web_summaries = web_summary
    Array[File] cellranger_metrics = metrics_summary
  }
}
