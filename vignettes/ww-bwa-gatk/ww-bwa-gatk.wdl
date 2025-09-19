version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bwa/ww-bwa.wdl" as bwa_tasks
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gatk/ww-gatk.wdl" as gatk_tasks

struct BwaSample {
    String name
    File reads
    File? mates
}

workflow bwa_gatk {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "WDL workflow to align sequencing data using BWA and mark duplicate reads with GATK"
    url: "https://github.com/getwilds/wilds-wdl-library/vignettes/ww-bwa-gatk"
    outputs: {
        markdup_bam: "Array of duplicate-marked bam files for each sample",
        markdup_bai: "Array of corresponding index files for each duplicate-marked bam file",
        duplicate_metrics: "Array of duplicate marking statistics for each sample"
    }
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file"
    samples: "List of BwaSample objects, each containing the sample name, forward reads FASTQ, and optionally reverse reads FASTQs"
    cpu_cores: "Number of CPUs to use for BWA alignment and GATK processing"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File reference_fasta
    Array[BwaSample] samples
    Int cpu_cores = 6
    Int memory_gb = 12
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
  }

  output {
    Array[File] markdup_bam = mark_duplicates.markdup_bam
    Array[File] markdup_bai = mark_duplicates.markdup_bai
    Array[File] duplicate_metrics = mark_duplicates.duplicate_metrics
  }
}
