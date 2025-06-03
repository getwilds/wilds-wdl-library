## WILDS WDL for performing sequence alignment using BWA-MEM.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

struct SampleInfo {
    String name
    File r1
    File r2
}

workflow bwa_example {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "WDL workflow for sequence alignment via BWA-MEM"
    url: "https://github.com/getwilds/wilds-wdl-library/modules/ww-bwa"
    outputs: {
        bwa_bam: "Sorted BWA-MEM alignment output BAM files for each sample",
        bwa_bai: "Index files for the sorted BWA-MEM alignment BAM files"
    }
  }

  parameter_meta {
    samples: "List of sample objects, each containing sample name and R1/R2 FASTQ files"
    reference_fasta: "Reference genome FASTA file"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    memory_gb: "Memory allocated for each task in the workflow in GB"
  }

  input {
    Array[SampleInfo] samples
    File reference_fasta
    Int cpus = 8
    Int memory_gb = 32
  }

  call bwa_index { input:
      reference_fasta = reference_fasta,
      cpu_cores = cpus,
      memory_gb = memory_gb
  }

  scatter (sample in samples) {
    call bwa_mem { input:
        sample_data = sample,
        reference_fasta = bwa_index.fasta,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  output {
    Array[File] bwa_bam = bwa_mem.sorted_bam
    Array[File] bwa_bai = bwa_mem.sorted_bai
  }
}

task bwa_index {
  meta {
    description: "Task for building BWA index files from a reference FASTA"
    outputs: {
        fasta: "Reference genome FASTA file",
        fa_amb: "Text file of ambiguous bases",
        fa_ann: "Text file of reference sequence information, such as name and length",
        fa_bwt: "Binary file of Burrows-Wheeler transformed reference sequence",
        fa_pac: "Binary file of compressed reference sequence",
        fa_sa: "Binary file of the suffix array"
    }
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File reference_fasta
    Int cpu_cores = 8
    Int memory_gb = 32
  }

  command <<<
  set -eo pipefail && \
  bwa index "~{reference_fasta}"
  >>>

  output {
    File fasta = "~{reference_fasta}"
    File fa_amb = "~{reference_fasta}.amb"
    File fa_ann = "~{reference_fasta}.ann"
    File fa_bwt = "~{reference_fasta}.bwt"
    File fa_pac = "~{reference_fasta}.pac"
    File fa_sa = "~{reference_fasta}.sa"
  }

  runtime {
    docker: "getwilds/bwa:0.7.17"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task bwa_mem {
  meta {
    description: "Task for aligning sequence reads using BWA-MEM"
    outputs: {
        sorted_bam: "Sorted BWA-MEM alignment output BAM file",
        sorted_bai: "Index files for the sorted BWA-MEM alignment BAM files"
    }
  }

  parameter_meta {
    sample_data: "Sample object containing sample name and R1/R2 FASTQ files"
    reference_fasta: "BWA indexed reference genome FASTA file"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    SampleInfo sample_data
    File reference_fasta
    Int cpu_cores = 8
    Int memory_gb = 16
  }

  command <<<
  set -eo pipefail && \
  bwa mem -v 3 -t max(1, ~{cpu_cores-1}) -M -R '@RG\tID:~{sample_data.name}\tSM:~{sample_data.name}\tPL:illumina' ~{reference_fasta} "~{sample_data.r1}" "~{sample_data.r2}" - > ~{sample_data.name}.sam && \
  samtools sort -@ max(1, ~{cpu_cores-1}) -o ~{sample_data.name}.sorted_aligned.bam ~{sample_data.name}.sam
  samtools index ~{sample_data.name}.sorted_aligned.bam
  >>>

  output {
    File sorted_bam = "~{sample_data.name}.sorted_aligned.bam"
    File sorted_bai = "~{sample_data.name}.sorted_aligned.bam.bai"
  }

  runtime {
    docker: "getwilds/bwa:0.7.17"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}