version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/sra-dbgap/modules/ww-sra/ww-sra.wdl" as sra_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-star/ww-star.wdl" as star_tasks

struct RefGenome {
    String name
    File fasta
    File gtf
}

workflow sra_star {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow to download raw sequencing data from SRA and align using STAR two-pass methodology"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-sra-star/ww-sra-star.wdl"
    outputs: {
        star_bam: "array of aligned bam files for each sample",
        star_bai: "array of corresponding index files for each aligned bam file",
        star_gene_counts: "array of text files containing the number of reads in each gene for each sample",
        star_log_final: "array of text files containing an overarching summary of the analysis performed for each sample",
        star_log_progress: "array of text files containing a detailed progress report for each sample",
        star_log: "array of text files containing STAR's raw command line output for each sample",
        star_sj: "array of text files containing splice junction details for each sample being analyzed"
    }
  }

  parameter_meta {
    sra_id_list: "list of SRA sample ID's to be pulled down and aligned"
    ref_genome: "reference genome object containing name, fasta, and gtf files"
    sjdb_overhang: "Length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database"
    genome_sa_index_nbases: "Length (bases) of the SA pre-indexing string, typically between 10-15 (scales with genome size)"
    ncpu: "number of CPUs to use for SRA download and STAR alignment"
    memory_gb: "memory allocation in GB for STAR tasks"
    max_reads: "Optional maximum number of reads to download from SRA (for testing/downsampling). If not specified, downloads all reads."
    ngc_file: "Optional NGC repository key file for downloading controlled-access dbGaP data."
  }

  input {
    Array[String] sra_id_list
    RefGenome ref_genome
    Int sjdb_overhang = 100
    Int genome_sa_index_nbases = 14
    Int ncpu = 12
    Int memory_gb = 64
    Int? max_reads
    File? ngc_file
  }

  call star_tasks.build_index { input:
      reference_fasta = ref_genome.fasta,
      reference_gtf = ref_genome.gtf,
      sjdb_overhang = sjdb_overhang,
      genome_sa_index_nbases = genome_sa_index_nbases,
      memory_gb = memory_gb,
      cpu_cores = ncpu
  }

  scatter ( id in sra_id_list ){
    call sra_tasks.fasterqdump { input:
        sra_id = id,
        ncpu = ncpu,
        max_reads = max_reads,
        ngc_file = ngc_file
    }

    # Paired-end alignment
    if (fasterqdump.is_paired_end) {
      call star_tasks.align_two_pass as align_paired { input:
          star_genome_tar = build_index.star_index_tar,
          name = id,
          r1 = fasterqdump.r1_end,
          r2 = fasterqdump.r2_end,
          sjdb_overhang = sjdb_overhang,
          memory_gb = memory_gb,
          cpu_cores = ncpu,
          star_threads = ncpu
      }
    }

    # Single-end alignment
    if (!fasterqdump.is_paired_end) {
      call star_tasks.align_two_pass as align_single { input:
          star_genome_tar = build_index.star_index_tar,
          name = id,
          r1 = fasterqdump.r1_end,
          sjdb_overhang = sjdb_overhang,
          memory_gb = memory_gb,
          cpu_cores = ncpu,
          star_threads = ncpu
      }
    }

    File bam_out = select_first([align_paired.bam, align_single.bam])
    File bai_out = select_first([align_paired.bai, align_single.bai])
    File gene_counts_out = select_first([align_paired.gene_counts, align_single.gene_counts])
    File log_final_out = select_first([align_paired.log_final, align_single.log_final])
    File log_progress_out = select_first([align_paired.log_progress, align_single.log_progress])
    File log_out = select_first([align_paired.log, align_single.log])
    File sj_out = select_first([align_paired.sj_out, align_single.sj_out])
  }

  output {
    Array[File] star_bam = bam_out
    Array[File] star_bai = bai_out
    Array[File] star_gene_counts = gene_counts_out
    Array[File] star_log_final = log_final_out
    Array[File] star_log_progress = log_progress_out
    Array[File] star_log = log_out
    Array[File] star_sj = sj_out
  }
}

