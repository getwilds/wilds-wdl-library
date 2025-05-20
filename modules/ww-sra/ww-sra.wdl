## Pulls down paired fastq's for SRA ID's provided

version 1.0

# TODO: add metadata pull using efetch
# TODO: add checksum for validation purposes
# TODO: add logging for better visibility

#### WORKFLOW DEFINITION

workflow sra_download {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow to download raw sequencing data from SRA in parallel"
    url: "https://github.com/getwilds/ww-sra"
    outputs: {
        r1_fastqs: "array of R1 fastq files for each sample",
        r2_fastqs: "array of R2 fastq files for each sample",
        is_paired_end: "array of booleans indicating whether each sample used paired-end sequencing"
    }
  }

  parameter_meta {
    sra_id_list: "list of SRA sample IDs to be pulled down"
    n_cpu: "number of cpus to use during download"
  }

  input {
    Array[String] sra_id_list
    Int n_cpu = 8
  }

  scatter (id in sra_id_list) {
    call fastqdump { input:
        sra_id = id,
        ncpu = n_cpu,
    }
  }

  output {
    Array[File] r1_fastqs = fastqdump.r1_end
    Array[File] r2_fastqs = fastqdump.r2_end
    Array[Boolean] is_paired_end = fastqdump.is_paired_end
  }
}

#### TASK DEFINITIONS

task fastqdump {
  meta {
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
  }

  input {
    String sra_id
    Int ncpu = 8
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
        --gzip
    else
      echo false > paired_file
      parallel-fastq-dump \
        --sra-id "~{sra_id}" \
        --threads ~{ncpu} \
        --outdir ./ \
        --gzip
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
