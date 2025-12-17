## WILDS WDL for de novo genome assembly using SPAdes.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task metaspades {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "De novo metagenomic assembly using metaSPAdes"
    outputs: {
        scaffolds_fasta: "Assembled scaffolds in compressed FASTA format",
        contigs_fasta: "Assembled contigs in compressed FASTA format",
        log_file: "SPAdes log file"
    }
  }

  parameter_meta {
    r1_fastq: "Read 1 FASTQ file"
    r2_fastq: "Read 2 FASTQ file"
    interleaved_fastq: "Interleaved FASTQ file of paired-end reads"
    sample_name: "Sample name for output file naming"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File? r1_fastq
    File? r2_fastq
    File? interleaved_fastq
    String sample_name
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Validate input: must have either interleaved OR both R1/R2
    if [[ -z "~{interleaved_fastq}" && ( -z "~{r1_fastq}" || -z "~{r2_fastq}" ) ]]; then
      echo "Error: Must provide either interleaved_fastq OR both r1_fastq and r2_fastq"
      exit 1
    fi

    if [[ -n "~{interleaved_fastq}" && ( -n "~{r1_fastq}" || -n "~{r2_fastq}" ) ]]; then
      echo "Error: Cannot provide both interleaved_fastq and r1_fastq/r2_fastq"
      exit 1
    fi

    # Build SPAdes command
    spades_cmd="spades.py \
      -o spades_outdir \
      --threads ~{cpu_cores} \
      --memory ~{memory_gb} \
      --meta \
      --only-assembler"

    # Add input files
    if [[ -n "~{interleaved_fastq}" ]]; then
      spades_cmd="$spades_cmd --12 ~{interleaved_fastq}"
    else
      spades_cmd="$spades_cmd -1 ~{r1_fastq} -2 ~{r2_fastq}"
    fi

    echo "Running: $spades_cmd"
    eval "$spades_cmd"

    # Move and compress scaffolds and contigs
    mv spades_outdir/scaffolds.fasta "~{sample_name}.scaffolds.fasta"
    gzip "~{sample_name}.scaffolds.fasta"

    mv spades_outdir/contigs.fasta "~{sample_name}.contigs.fasta"
    gzip "~{sample_name}.contigs.fasta"

    # Move log file
    mv spades_outdir/spades.log "~{sample_name}.spades.log"

    # Clean up SPAdes output directory to save space
    rm -rf spades_outdir
  >>>

  output {
    File scaffolds_fasta = "~{sample_name}.scaffolds.fasta.gz"
    File contigs_fasta = "~{sample_name}.contigs.fasta.gz"
    File log_file = "~{sample_name}.spades.log"
  }

  runtime {
    docker: "getwilds/spades:4.2.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
