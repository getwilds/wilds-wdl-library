## WILDS WDL for processing genomic files with Samtools.
## Intended for use alone or as a modular component in the WILDS ecosystem.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow samtools_example {
  meta {
    author: "Emma Bishop"
    email: "wilds@fredhutch.org"
    description: "WDL workflow for processing genomic files with Samtools"
    url: "https://github.com/getwilds/wilds-wdl-library/modules/ww-samtools"
    outputs: {
        r1_fastqs: "R1 FASTQ files generated from CRAM/BAM/SAM files",
        r2_fastqs: "R2 FASTQ files generated from CRAM/BAM/SAM files"
    }
  }

  # Download test data
  call ww_testdata.download_ref_data { }

  call ww_testdata.download_cram_data { input:
      ref_fasta = download_ref_data.fasta
  }

  # Convert CRAM to FASTQ
  call crams_to_fastq { input:
      cram_files = [download_cram_data.cram],
      ref = download_ref_data.fasta,
      name = "test_sample",
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    File r1_fastqs = crams_to_fastq.r1_fastq
    File r2_fastqs = crams_to_fastq.r2_fastq
  }
}

task crams_to_fastq {
  meta {
    description: "Merge CRAM/BAM/SAM files and convert to FASTQ's using samtools."
    outputs: {
        r1_fastq: "R1 FASTQ file generated from merged CRAM/BAM/SAM file",
        r2_fastq: "R2 FASTQ file generated from merged CRAM/BAM/SAM file",
        sample_name: "Sample name that was processed"
    }
  }

  parameter_meta {
    cram_files: "List of CRAM/BAM/SAM files for a sample"
    ref: "Reference genome FASTA file"
    name: "Name of the sample (used for file naming)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File ref
    Array[File] cram_files
    String name
    Int cpu_cores = 2
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    # Merge CRAM/BAM/SAM files if more than one, then collate and convert to FASTQ
    samtools merge -@ ~{cpu_cores} --reference "~{ref}" -u - ~{sep=" " cram_files} | \
    samtools collate -@ ~{cpu_cores} --reference "~{ref}" -O - | \
    samtools fastq -@ ~{cpu_cores} --reference "~{ref}" -1 "~{name}_R1.fastq.gz" -2 "~{name}_R2.fastq.gz" -0 /dev/null -s /dev/null -
  >>>

  output {
    File r1_fastq = "~{name}_R1.fastq.gz"
    File r2_fastq = "~{name}_R2.fastq.gz"
    String sample_name = "~{name}"
  }

  runtime {
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
    docker: "getwilds/samtools:1.19"
  }
}

