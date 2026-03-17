## WILDS WDL fastp Module
## Ultra-fast all-in-one FASTQ preprocessor for quality control and adapter trimming
## fastp performs quality filtering, adapter trimming, and generates QC reports
## in a single pass over the data

version 1.0

#### TASK DEFINITIONS ####

task fastp_paired {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Run fastp on paired-end FASTQ files for quality filtering, adapter trimming, and QC reporting"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-fastp/ww-fastp.wdl"
    outputs: {
        r1_trimmed: "Trimmed and filtered R1 FASTQ file",
        r2_trimmed: "Trimmed and filtered R2 FASTQ file",
        html_report: "HTML quality control report",
        json_report: "JSON quality control report"
    }
  }

  parameter_meta {
    sample_name: "Name identifier for the sample"
    r1_fastq: "Read 1 input FASTQ file (gzipped or uncompressed)"
    r2_fastq: "Read 2 input FASTQ file (gzipped or uncompressed)"
    qualified_quality_phred: "Minimum base quality score for a base to be qualified (default: 15)"
    length_required: "Minimum read length after trimming (default: 15)"
    detect_adapter_for_pe: "Enable auto-detection of adapters for paired-end data (default: true)"
    adapter_fasta: "Optional FASTA file with custom adapter sequences"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    String sample_name
    File r1_fastq
    File r2_fastq
    Int qualified_quality_phred = 15
    Int length_required = 15
    Boolean detect_adapter_for_pe = true
    File? adapter_fasta
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    FASTP_CMD="fastp \
      --in1 ~{r1_fastq} \
      --in2 ~{r2_fastq} \
      --out1 ~{sample_name}_R1_trimmed.fastq.gz \
      --out2 ~{sample_name}_R2_trimmed.fastq.gz \
      --html ~{sample_name}_fastp.html \
      --json ~{sample_name}_fastp.json \
      --qualified_quality_phred ~{qualified_quality_phred} \
      --length_required ~{length_required} \
      --thread ~{cpu_cores}"

    if [ "~{detect_adapter_for_pe}" == "true" ]; then
      FASTP_CMD="${FASTP_CMD} --detect_adapter_for_pe"
    fi

    if [ ! -z "~{adapter_fasta}" ]; then
      FASTP_CMD="${FASTP_CMD} --adapter_fasta ~{adapter_fasta}"
    fi

    echo "Running: ${FASTP_CMD}"
    ${FASTP_CMD}
  >>>

  output {
    File r1_trimmed = "~{sample_name}_R1_trimmed.fastq.gz"
    File r2_trimmed = "~{sample_name}_R2_trimmed.fastq.gz"
    File html_report = "~{sample_name}_fastp.html"
    File json_report = "~{sample_name}_fastp.json"
  }

  runtime {
    docker: "getwilds/fastp:1.1.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task fastp_single {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Run fastp on single-end FASTQ files for quality filtering, adapter trimming, and QC reporting"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-fastp/ww-fastp.wdl"
    outputs: {
        trimmed_fastq: "Trimmed and filtered FASTQ file",
        html_report: "HTML quality control report",
        json_report: "JSON quality control report"
    }
  }

  parameter_meta {
    sample_name: "Name identifier for the sample"
    fastq: "Input FASTQ file (gzipped or uncompressed)"
    qualified_quality_phred: "Minimum base quality score for a base to be qualified (default: 15)"
    length_required: "Minimum read length after trimming (default: 15)"
    adapter_fasta: "Optional FASTA file with custom adapter sequences"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    String sample_name
    File fastq
    Int qualified_quality_phred = 15
    Int length_required = 15
    File? adapter_fasta
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    FASTP_CMD="fastp \
      --in1 ~{fastq} \
      --out1 ~{sample_name}_trimmed.fastq.gz \
      --html ~{sample_name}_fastp.html \
      --json ~{sample_name}_fastp.json \
      --qualified_quality_phred ~{qualified_quality_phred} \
      --length_required ~{length_required} \
      --thread ~{cpu_cores}"

    if [ ! -z "~{adapter_fasta}" ]; then
      FASTP_CMD="${FASTP_CMD} --adapter_fasta ~{adapter_fasta}"
    fi

    echo "Running: ${FASTP_CMD}"
    ${FASTP_CMD}
  >>>

  output {
    File trimmed_fastq = "~{sample_name}_trimmed.fastq.gz"
    File html_report = "~{sample_name}_fastp.html"
    File json_report = "~{sample_name}_fastp.json"
  }

  runtime {
    docker: "getwilds/fastp:1.1.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
