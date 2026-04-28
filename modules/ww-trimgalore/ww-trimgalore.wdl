## WILDS WDL Trim Galore Module
## Wrapper around Cutadapt and FastQC for adapter/quality trimming of FASTQ files.
## Trim Galore auto-detects adapters and performs quality trimming in one step,
## with optional FastQC reports after trimming.

version 1.0

#### TASK DEFINITIONS ####

task trimgalore_paired {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Run Trim Galore on paired-end FASTQ files for adapter and quality trimming"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-trimgalore/ww-trimgalore.wdl"
    outputs: {
        r1_trimmed: "Trimmed and validated R1 FASTQ file",
        r2_trimmed: "Trimmed and validated R2 FASTQ file",
        r1_report: "Trimming report for R1",
        r2_report: "Trimming report for R2"
    }
  }

  parameter_meta {
    sample_name: "Name identifier for the sample"
    r1_fastq: "Read 1 input FASTQ file (gzipped or uncompressed)"
    r2_fastq: "Read 2 input FASTQ file (gzipped or uncompressed)"
    quality: "Phred quality score threshold for trimming low-quality ends (default: 20)"
    length: "Minimum read length after trimming; shorter reads are discarded (default: 20)"
    stringency: "Minimum overlap with adapter sequence required to trim (default: 1)"
    run_fastqc: "Run FastQC on trimmed files (default: false)"
    adapter: "Optional adapter sequence for R1 (auto-detected if not specified)"
    adapter2: "Optional adapter sequence for R2 (auto-detected if not specified)"
    cpu_cores: "Number of CPU cores allocated for the task (default: 4)"
    memory_gb: "Memory allocated for the task in GB (default: 8)"
  }

  input {
    String sample_name
    File r1_fastq
    File r2_fastq
    Int quality = 20
    Int length = 20
    Int stringency = 1
    Boolean run_fastqc = false
    String? adapter
    String? adapter2
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    TG_CMD="trim_galore \
      --paired \
      --quality ~{quality} \
      --length ~{length} \
      --stringency ~{stringency} \
      --cores ~{cpu_cores} \
      --basename ~{sample_name} \
      --gzip \
      --output_dir ."

    if [ "~{run_fastqc}" == "true" ]; then
      TG_CMD="${TG_CMD} --fastqc"
    fi

    if [ ! -z "~{adapter}" ]; then
      TG_CMD="${TG_CMD} --adapter ~{adapter}"
    fi

    if [ ! -z "~{adapter2}" ]; then
      TG_CMD="${TG_CMD} --adapter2 ~{adapter2}"
    fi

    TG_CMD="${TG_CMD} ~{r1_fastq} ~{r2_fastq}"

    echo "Running: ${TG_CMD}"
    ${TG_CMD}
  >>>

  output {
    File r1_trimmed = "~{sample_name}_val_1.fq.gz"
    File r2_trimmed = "~{sample_name}_val_2.fq.gz"
    File r1_report = glob("*trimming_report.txt")[0]
    File r2_report = glob("*trimming_report.txt")[1]
  }

  runtime {
    docker: "getwilds/trim-galore:0.6.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task trimgalore_single {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Run Trim Galore on a single-end FASTQ file for adapter and quality trimming"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-trimgalore/ww-trimgalore.wdl"
    outputs: {
        trimmed_fastq: "Trimmed FASTQ file",
        trimming_report: "Trimming report"
    }
  }

  parameter_meta {
    sample_name: "Name identifier for the sample"
    fastq: "Input FASTQ file (gzipped or uncompressed)"
    quality: "Phred quality score threshold for trimming low-quality ends (default: 20)"
    length: "Minimum read length after trimming; shorter reads are discarded (default: 20)"
    stringency: "Minimum overlap with adapter sequence required to trim (default: 1)"
    run_fastqc: "Run FastQC on trimmed files (default: false)"
    adapter: "Optional adapter sequence (auto-detected if not specified)"
    cpu_cores: "Number of CPU cores allocated for the task (default: 4)"
    memory_gb: "Memory allocated for the task in GB (default: 8)"
  }

  input {
    String sample_name
    File fastq
    Int quality = 20
    Int length = 20
    Int stringency = 1
    Boolean run_fastqc = false
    String? adapter
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    TG_CMD="trim_galore \
      --quality ~{quality} \
      --length ~{length} \
      --stringency ~{stringency} \
      --cores ~{cpu_cores} \
      --basename ~{sample_name} \
      --gzip \
      --output_dir ."

    if [ "~{run_fastqc}" == "true" ]; then
      TG_CMD="${TG_CMD} --fastqc"
    fi

    if [ ! -z "~{adapter}" ]; then
      TG_CMD="${TG_CMD} --adapter ~{adapter}"
    fi

    TG_CMD="${TG_CMD} ~{fastq}"

    echo "Running: ${TG_CMD}"
    ${TG_CMD}
  >>>

  output {
    File trimmed_fastq = "~{sample_name}_trimmed.fq.gz"
    File trimming_report = glob("*_trimming_report.txt")[0]
  }

  runtime {
    docker: "getwilds/trim-galore:0.6.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
