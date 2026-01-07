## WILDS WDL module for downloading FASTQ data from NCBI's Sequence Read Archive (SRA).
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task fastqdump {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Task for pulling down fastq data from SRA."
    outputs: {
        r1_end: "R1 fastq file downloaded for the sample in question",
        r2_end: "R2 fastq file downloaded for the sample in question (empty file for single-end reads)",
        is_paired_end: "boolean indicating whether the sample used paired-end sequencing"
    }
  }

  parameter_meta {
    sra_id: "SRA ID of the sample to be downloaded via parallel-fastq-dump"
    ncpu: "number of cpus to use during download"
    max_reads: "Optional maximum number of reads to download (for testing/downsampling). If not specified, downloads all reads."
  }

  input {
    String sra_id
    Int ncpu = 8
    Int? max_reads
  }

  command <<<
    set -eo pipefail
    # check if paired ended
    numLines=$(fastq-dump -X 1 -Z --split-spot "~{sra_id}" | wc -l)
    if [ "$numLines" -eq 8 ]; then
      echo true > paired_file
      parallel-fastq-dump \
        --sra-id "~{sra_id}" \
        --threads ~{ncpu} \
        --outdir ./ \
        --split-files \
        --gzip \
        ~{if defined(max_reads) then "--maxSpotId " + max_reads else ""}
    else
      echo false > paired_file
      parallel-fastq-dump \
        --sra-id "~{sra_id}" \
        --threads ~{ncpu} \
        --outdir ./ \
        --gzip \
        ~{if defined(max_reads) then "--maxSpotId " + max_reads else ""}
      # Rename the file to match the expected output format
      mv "~{sra_id}.fastq.gz" "~{sra_id}_1.fastq.gz"
      # Create an empty placeholder for R2
      touch "~{sra_id}_2.fastq.gz"
    fi
  >>>

  output {
    File r1_end = "~{sra_id}_1.fastq.gz"
    File r2_end = "~{sra_id}_2.fastq.gz"
    Boolean is_paired_end = read_boolean("paired_file")
  }

  runtime {
    memory: 2 * ncpu + " GB"
    docker: "getwilds/sra-tools:3.1.1"
    cpu: ncpu
  }
}
