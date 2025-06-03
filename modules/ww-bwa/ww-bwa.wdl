## WILDS WDL for performing sequence alignment using BWA-MEM.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

struct SampleInfo {
    String name
    File r1
    File r2
}

struct RefGenome {
    File ref_dict
    File ref_fasta
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_fasta_index
    File ref_pac
    File ref_sa
}
 
workflow bwa_example {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "WDL workflow for sequence alignment via BWA-MEM"
    url: "https://github.com/getwilds/wilds-wdl-library/modules/ww-bwa"
    outputs: {
        bwa_bam: "Sorted BWA-MEM alignment output BAM files for each sample"
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
    RefGenome reference_genome
    Int cpus = 8
    Int memory_gb = 64
  }

  scatter (sample in samples) {
    call bwa_mem { input:
        sample_data = sample,
        reference_fasta = reference_genome.ref_fasta,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  output {
    Array[File] bwa_bam = bwa_mem.sorted_bam
  }
}

task bwa_mem {
  meta {
    description: "Task for aligning sequence reads using BWA-MEM"
    outputs: {
        sorted_bam: "Sorted BWA-MEM alignment output BAM file"
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
  bwa mem -p -v 3 -t ~{cpu_cores} -M -R '@RG\tID:~{sample_data.name}\tSM:~{sample_data.name}\tPL:illumina' ~{reference_fasta} "~{sample_data.r1}" "~{sample_data.r2}" - > ~{sample_data.name}.sam && \
  samtools sort -@ ~{cpu_cores-1} -o ~{sample_data.name}.sorted_aligned.bam ~{sample_data.name}.sam
  >>>

  output {
    File sorted_bam = "~{sample_data.name}.sorted_aligned.bam"
  }

  runtime {
    docker: "getwilds/bwa:0.7.17"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}