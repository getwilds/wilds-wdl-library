version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bwa/ww-bwa.wdl" as bwa_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/cirro-bwa-gatk/modules/ww-gatk/ww-gatk.wdl" as gatk_tasks

struct BwaSample {
    String name
    File reads
    File? mates
}

workflow bwa_gatk {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "WDL workflow to align sequencing data using BWA and perform QC steps with GATK"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-bwa-gatk/ww-bwa-gatk.wdl"
    outputs: {
        duplicate_metrics: "Array of duplicate marking statistics for each sample",
        recalibrated_bam: "Array of BAM files with recalibrated base quality scores",
        recalibrated_bai: "Array of corresponding index files for each recalibrated bam file",
        recalibration_report: "Array of base recalibration report tables"
    }
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index for reference genome FASTA file"
    dbsnp_vcf: "dbSNP VCF file for known variant sites"
    known_indels_vcf: "Array of VCF files with known indel sites"
    samples: "List of BwaSample objects, each containing the sample name, forward reads FASTQ, and optionally reverse reads FASTQs"
    cpu_cores: "Number of CPUs to use for BWA alignment and GATK processing"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File dbsnp_vcf
    Array[File] known_indels_vcf
    Array[BwaSample] samples
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  # Build BWA index once for all samples
  call bwa_tasks.bwa_index { input:
      reference_fasta = reference_fasta,
      cpu_cores = cpu_cores,
      memory_gb = memory_gb
  }

  scatter ( sample in samples ){
    call bwa_tasks.bwa_mem { input:
        bwa_genome_tar = bwa_index.bwa_index_tar,
        reference_fasta = reference_fasta,
        reads = sample.reads,
        name = sample.name,
        mates = sample.mates,
        paired_end = defined(sample.mates),
        cpu_cores = cpu_cores,
        memory_gb = memory_gb
    }

    call gatk_tasks.mark_duplicates { input:
        bam = bwa_mem.sorted_bam,
        bam_index = bwa_mem.sorted_bai,
        base_file_name = sample.name,
        memory_gb = memory_gb,
        cpu_cores = cpu_cores
    }

    call gatk_tasks.create_sequence_dictionary { input:
        reference_fasta = reference_fasta
    }

    call gatk_tasks.base_recalibrator { input:
        bam = mark_duplicates.markdup_bam,
        bam_index = mark_duplicates.markdup_bai,
        dbsnp_vcf = dbsnp_vcf,
        reference_fasta = reference_fasta,
        reference_fasta_index = reference_fasta_index,
        reference_dict = create_sequence_dictionary.sequence_dict,
        known_indels_sites_vcfs = known_indels_vcf,
        base_file_name = sample.name,
        memory_gb = memory_gb,
        cpu_cores = cpu_cores
    }
  }

  output {
    Array[File] duplicate_metrics = mark_duplicates.duplicate_metrics
    Array[File] recalibrated_bam = base_recalibrator.recalibrated_bam
    Array[File] recalibrated_bai = base_recalibrator.recalibrated_bai
    Array[File] recalibration_report = base_recalibrator.recalibration_report
  }
}
