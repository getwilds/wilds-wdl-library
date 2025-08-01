version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl" as sra_tasks
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
    url: "https://github.com/getwilds/wilds-wdl-library/vignettes/ww-sra-star"
    outputs: {
        star_bam: "array of aligned bam files for each sample",
        star_bai: "array of corresponding index files for each aligned bam file",
        star_gene_counts: "array of text files containing the number of reads in each gene for each sample",
        star_log_final: "array of text files containing an overarching summary of the analysis performed for each sample",
        star_log_progress: "array of text files containing a detailed progress report for each sample",
        star_log: "array of text files containing STAR's raw command line output for each sample",
        star_sj: "array of text files containing splice junction details for each sample being analyzed",
        validation_report: "validation report confirming all expected outputs were generated correctly"
    }
  }

  parameter_meta {
    sra_id_list: "list of SRA sample ID's to be pulled down and aligned"
    ref_genome: "reference genome object containing name, fasta, and gtf files"
    sjdb_overhang: "Length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database"
    genome_sa_index_nbases: "Length (bases) of the SA pre-indexing string, typically between 10-15 (scales with genome size)"
    ncpu: "number of CPUs to use for SRA download and STAR alignment"
    memory_gb: "memory allocation in GB for STAR tasks"
  }

  input {
    Array[String] sra_id_list
    RefGenome ref_genome
    Int sjdb_overhang = 100
    Int genome_sa_index_nbases = 14
    Int ncpu = 12
    Int memory_gb = 64
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
    call sra_tasks.fastqdump { input:
        sra_id = id,
        ncpu = ncpu
    }

    call star_tasks.align_two_pass { input:
        star_genome_tar = build_index.star_index_tar,
        name = id,
        r1 = fastqdump.r1_end,
        r2 = fastqdump.r2_end,
        sjdb_overhang = sjdb_overhang,
        memory_gb = memory_gb,
        cpu_cores = ncpu,
        star_threads = ncpu
    }
  }

  call star_tasks.validate_outputs { input:
    bam_files = align_two_pass.bam,
    bai_files = align_two_pass.bai,
    gene_count_files = align_two_pass.gene_counts
  }

  output {
    Array[File] star_bam = align_two_pass.bam
    Array[File] star_bai = align_two_pass.bai
    Array[File] star_gene_counts = align_two_pass.gene_counts
    Array[File] star_log_final = align_two_pass.log_final
    Array[File] star_log_progress = align_two_pass.log_progress
    Array[File] star_log = align_two_pass.log
    Array[File] star_sj = align_two_pass.sj_out
    File validation_report = validate_outputs.report
  }
}
