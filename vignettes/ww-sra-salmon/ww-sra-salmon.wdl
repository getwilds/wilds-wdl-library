version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl" as sra_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-salmon/modules/ww-salmon/ww-salmon.wdl" as salmon_tasks

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
    url: "https://github.com/getwilds/wilds-wdl-library/vignettes/ww-sra-salmon"
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
  }

  input {
    Array[String] sra_id_list
    File transcriptome_fasta
    Int ncpu = 8
    Int memory_gb = 16
    Int? max_reads
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
        max_reads = max_reads
    }

    call salmon_tasks.quantify { input:
        salmon_index_dir = build_index.salmon_index,
        sample_name = id,
        fastq_r1 = fastqdump.r1_end,
        fastq_r2 = fastqdump.r2_end,
        cpu_cores = ncpu,
        memory_gb = memory_gb
    }
  }

  # Merge quantification results across all samples
  call salmon_tasks.merge_results { input:
      salmon_quant_dirs = quantify.salmon_quant_dir,
      sample_names = quantify.output_sample_name,
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    Array[File] salmon_quant_dirs = quantify.salmon_quant_dir
    File tpm_matrix = merge_results.tpm_matrix
    File counts_matrix = merge_results.counts_matrix
  }
}
