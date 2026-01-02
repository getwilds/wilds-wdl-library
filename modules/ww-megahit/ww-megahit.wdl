## WILDS WDL for de novo metagenome assembly using MEGAHIT.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task megahit {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "De novo metagenome assembly using MEGAHIT"
    outputs: {
        contigs: "Assembled contigs in compressed FASTA format"
    }
  }

  parameter_meta {
    input_fastq: "Interleaved FASTQ file of paired-end reads"
    sample_name: "Sample name for output file naming"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File input_fastq
    String sample_name = basename(input_fastq, ".fastq.gz")
    Int cpu_cores = 2
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    megahit --12 ~{input_fastq} \
            -o megahit_outdir \
            -t ~{cpu_cores}

    # Move and compress contigs
    mv megahit_outdir/final.contigs.fa "~{sample_name}.contigs.fasta"
    gzip "~{sample_name}.contigs.fasta"

    # Clean up MEGAHIT output directory to save space
    rm -rf megahit_outdir
  >>>

  output {
    File contigs = "~{sample_name}.contigs.fasta.gz"
  }

  runtime {
    docker: "megahit:latest"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
