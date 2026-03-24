version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/fix-sra-star-jyoung/modules/ww-ena/ww-ena.wdl" as ena_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/fix-sra-star-jyoung/modules/ww-star/ww-star.wdl" as star_tasks

struct RefGenome {
    String name
    File fasta
    File gtf
}

workflow ena_star {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow to download raw sequencing data from ENA and align using STAR two-pass methodology"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-ena-star/ww-ena-star.wdl"
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
    ena_accession_list: "list of ENA accession ID's to be pulled down and aligned"
    ref_genome: "reference genome object containing name, fasta, and gtf files"
    sjdb_overhang: "Length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database"
    genome_sa_index_nbases: "Length (bases) of the SA pre-indexing string, typically between 10-15 (scales with genome size)"
    ncpu: "number of CPUs to use for ENA download and STAR alignment"
    memory_gb: "memory allocation in GB for STAR tasks"
    ena_file_format: "Format of files to download from ENA (default: READS_FASTQ)"
    ena_protocol: "Transfer protocol for ENA download: FTP or ASPERA (default: FTP)"
  }

  input {
    Array[String] ena_accession_list
    RefGenome ref_genome
    Int sjdb_overhang = 100
    Int genome_sa_index_nbases = 14
    Int ncpu = 12
    Int memory_gb = 64
    String ena_file_format = "READS_FASTQ"
    String ena_protocol = "FTP"
  }

  call star_tasks.build_index { input:
      reference_fasta = ref_genome.fasta,
      reference_gtf = ref_genome.gtf,
      sjdb_overhang = sjdb_overhang,
      genome_sa_index_nbases = genome_sa_index_nbases,
      memory_gb = memory_gb,
      cpu_cores = ncpu
  }

  scatter ( accession in ena_accession_list ){
    call ena_tasks.download_files { input:
        accessions = accession,
        file_format = ena_file_format,
        protocol = ena_protocol,
        output_dir_name = "ena_download_" + accession,
        cpu_cores = ncpu,
        memory_gb = 8
    }

    call ena_tasks.extract_fastq_pairs { input:
        downloaded_files = download_files.downloaded_files
    }

    # Since we scatter by single accession, each extract_fastq_pairs call returns
    # arrays with one element each - access the first element with [0]
    Boolean sample_is_paired = extract_fastq_pairs.is_paired_end_list[0] == "true"

    # Paired-end alignment
    if (sample_is_paired) {
      call star_tasks.align_two_pass as align_paired { input:
          star_genome_tar = build_index.star_index_tar,
          name = extract_fastq_pairs.accessions[0],
          r1 = extract_fastq_pairs.r1_files[0],
          r2 = extract_fastq_pairs.r2_files[0],
          sjdb_overhang = sjdb_overhang,
          memory_gb = memory_gb,
          cpu_cores = ncpu,
          star_threads = ncpu
      }
    }

    # Single-end alignment
    if (!sample_is_paired) {
      call star_tasks.align_two_pass as align_single { input:
          star_genome_tar = build_index.star_index_tar,
          name = extract_fastq_pairs.accessions[0],
          r1 = extract_fastq_pairs.r1_files[0],
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
