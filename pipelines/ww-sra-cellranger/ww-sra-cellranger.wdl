version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-sra-cellranger/modules/ww-sra/ww-sra.wdl" as sra_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-sra-cellranger/modules/ww-cellranger/ww-cellranger.wdl" as cellranger_tasks

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

    # Run Cell Ranger count
    call cellranger_tasks.run_count { input:
        r1_fastqs = [rename_fastqs.r1_renamed],
        r2_fastqs = [rename_fastqs.r2_renamed],
        ref_gex = ref_gex,
        sample_id = id,
        create_bam = create_bam,
        cpu_cores = ncpu,
        memory_gb = memory_gb,
        expect_cells = expect_cells,
        chemistry = chemistry
    }
  }

  output {
    Array[File] cellranger_results = run_count.results_tar
    Array[File] cellranger_web_summaries = run_count.web_summary
    Array[File] cellranger_metrics = run_count.metrics_summary
  }
}
