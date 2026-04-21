## WILDS WDL module for downloading FASTQ data from NCBI's Sequence Read Archive (SRA).
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.
## Supports both public SRA data and controlled-access dbGaP data via NGC authentication.

version 1.0

task fastqdump {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Task for pulling down fastq data from SRA using prefetch and fasterq-dump. Supports controlled-access dbGaP data via optional NGC repository key file."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl"
    outputs: {
        r1_end: "R1 fastq file downloaded for the sample in question",
        r2_end: "R2 fastq file downloaded for the sample in question (empty file for single-end reads)",
        is_paired_end: "boolean indicating whether the sample used paired-end sequencing"
    }
    topic: "any"
    species: "any"
    operation: "data_retrieval"
    in_sample_req: "none"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "r1_end:dna_sequence:fastq,r2_end:dna_sequence:fastq"
    out_ref: "none"
  }

  parameter_meta {
    sra_id: "SRA ID of the sample to be downloaded via prefetch and fasterq-dump"
    ncpu: "number of cpus to use during download"
    max_reads: "Optional maximum number of reads to download (for testing/downsampling). If not specified, downloads all reads."
    ngc_file: "Optional NGC repository key file for downloading controlled-access dbGaP data. Must be obtained through an approved dbGaP project."
  }

  input {
    String sra_id
    Int ncpu = 8
    Int? max_reads
    File? ngc_file
  }

  command <<<
    set -eo pipefail
    # Prefetch the SRA data (handles both public and dbGaP data)
    prefetch "~{sra_id}" \
      ~{if defined(ngc_file) then "--ngc " + ngc_file else ""}
    # Check if paired ended
    numLines=$(fastq-dump -X 1 -Z --split-spot "~{sra_id}" \
      ~{if defined(ngc_file) then "--ngc " + ngc_file else ""} | wc -l)
    if [ "$numLines" -gt 4 ]; then
      echo true > paired_file
      fasterq-dump "~{sra_id}" \
        --threads ~{ncpu} \
        --outdir ./ \
        --split-files \
        --include-technical \
        ~{if defined(ngc_file) then "--ngc " + ngc_file else ""}
    else
      echo false > paired_file
      fasterq-dump "~{sra_id}" \
        --threads ~{ncpu} \
        --outdir ./ \
        ~{if defined(ngc_file) then "--ngc " + ngc_file else ""}
      # Rename the file to match the expected output format
      mv "~{sra_id}.fastq" "~{sra_id}_1.fastq"
      # Create an empty placeholder for R2
      touch "~{sra_id}_2.fastq"
    fi
    # Handle cases where index reads are included: fasterq-dump --split-files produces 3 files where
    # _1 = R1, _2 = index read, _3 = R2
    # Reassign _3 to _2 so R1/R2 outputs always contain the biological reads.
    if [ -f "~{sra_id}_3.fastq" ]; then
      mv "~{sra_id}_3.fastq" "~{sra_id}_2.fastq"
    fi
    # Truncate to max_reads if specified (fasterq-dump has no built-in read limit)
    MAX_READS="~{select_first([max_reads, 0])}"
    if [ "$MAX_READS" -gt 0 ]; then
      max_lines=$((MAX_READS * 4))
      for f in "~{sra_id}_1.fastq" "~{sra_id}_2.fastq"; do
        if [ -s "$f" ]; then
          head -n $max_lines "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
        fi
      done
    fi
    # Gzip the output files (fasterq-dump does not support --gzip natively)
    gzip "~{sra_id}_1.fastq" "~{sra_id}_2.fastq"
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
