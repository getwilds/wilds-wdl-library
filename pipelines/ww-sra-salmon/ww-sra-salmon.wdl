version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl" as sra_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-salmon/ww-salmon.wdl" as salmon_tasks

struct SalmonSample {
    String name
    File r1_fastq
    File r2_fastq
}

workflow sra_salmon {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow to download raw sequencing data from SRA and quantify using Salmon quasi-mapping"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-sra-salmon/ww-sra-salmon.wdl"
    outputs: {
        salmon_quant_dirs: "array of Salmon quantification output directories for each sample",
        tpm_matrix: "matrix of TPM (transcripts per million) values across all samples",
        counts_matrix: "matrix of estimated counts across all samples"
    }
  }

  parameter_meta {
    sra_id_list: "list of SRA sample ID's to be pulled down and quantified"
    transcriptome_fasta: "reference transcriptome FASTA file for Salmon indexing and quantification"
    ncpu: "number of CPUs to use for SRA download and Salmon quantification"
    memory_gb: "memory allocation in GB for Salmon tasks"
    max_reads: "Optional maximum number of reads to download from SRA (for testing/downsampling). If not specified, downloads all reads."
    ngc_file: "Optional NGC repository key file for downloading controlled-access dbGaP data."
  }

  input {
    Array[String] sra_id_list
    File transcriptome_fasta
    Int ncpu = 8
    Int memory_gb = 16
    Int? max_reads
    File? ngc_file
  }

  # Build Salmon index from transcriptome
  call salmon_tasks.build_index { input:
      transcriptome_fasta = transcriptome_fasta,
      cpu_cores = ncpu,
      memory_gb = memory_gb
  }

  # Download FASTQ files from SRA and quantify each sample
  scatter ( id in sra_id_list ){
    call sra_tasks.fastqdump { input:
        sra_id = id,
        ncpu = ncpu,
        max_reads = max_reads,
        ngc_file = ngc_file
    }

    # Paired-end quantification
    if (fastqdump.is_paired_end) {
      call salmon_tasks.quantify as quantify_paired { input:
          salmon_index_dir = build_index.salmon_index,
          sample_name = id,
          fastq_r1 = fastqdump.r1_end,
          fastq_r2 = fastqdump.r2_end,
          cpu_cores = ncpu,
          memory_gb = memory_gb
      }
    }

    # Single-end quantification
    if (!fastqdump.is_paired_end) {
      call salmon_tasks.quantify as quantify_single { input:
          salmon_index_dir = build_index.salmon_index,
          sample_name = id,
          fastq_r1 = fastqdump.r1_end,
          cpu_cores = ncpu,
          memory_gb = memory_gb
      }
    }

    File quant_dir_out = select_first([quantify_paired.salmon_quant_dir, quantify_single.salmon_quant_dir])
    String sample_name_out = select_first([quantify_paired.output_sample_name, quantify_single.output_sample_name])
  }

  # Merge quantification results across all samples
  call salmon_tasks.merge_results { input:
      salmon_quant_dirs = quant_dir_out,
      sample_names = sample_name_out,
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    Array[File] salmon_quant_dirs = quant_dir_out
    File tpm_matrix = merge_results.tpm_matrix
    File counts_matrix = merge_results.counts_matrix
  }
}
