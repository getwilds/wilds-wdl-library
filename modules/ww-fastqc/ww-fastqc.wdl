## WILDS WDL FastQC Module
## Quality control analysis for high-throughput sequencing data
## FastQC is a tool for quality control checks on raw sequence data
## This module runs FastQC on FASTQ files to generate quality control reports

version 1.0

# Import testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### WORKFLOW DEFINITION ####

workflow fastqc_example {
  meta {
    author: "WILDS Development Team"
    email: "wilds@fredhutch.org"
    description: "FastQC quality control analysis for high-throughput sequencing data"
    url: "https://github.com/getwilds/wilds-wdl-library"
    outputs: {
        fastqc_html_reports: "Array of FastQC HTML quality control reports",
        fastqc_zip_reports: "Array of FastQC ZIP archives containing all report data"
    }
  }

  # Auto-download test data for testing purposes
  call ww_testdata.download_fastq_data as download_demo_data { }

  # Process paired-end sample
  call run_fastqc as run_fastqc_paired { input:
      r1_fastq = download_demo_data.r1_fastq,
      r2_fastq = download_demo_data.r2_fastq,
      cpu_cores = 2,
      memory_gb = 4
  }

  # Process single-end sample (R1 only)
  call run_fastqc as run_fastqc_single { input:
      r1_fastq = download_demo_data.r1_fastq,
      cpu_cores = 2,
      memory_gb = 4
  }

  output {
    Array[File] fastqc_html_reports = flatten([run_fastqc_paired.html_reports, run_fastqc_single.html_reports])
    Array[File] fastqc_zip_reports = flatten([run_fastqc_paired.zip_reports, run_fastqc_single.zip_reports])
  }
}

#### TASK DEFINITIONS ####

task run_fastqc {
  meta {
    description: "Run FastQC quality control analysis on FASTQ files"
    outputs: {
        html_reports: "FastQC HTML quality control reports",
        zip_reports: "FastQC ZIP archives containing all report data"
    }
  }

  parameter_meta {
    r1_fastq: "Read 1 FASTQ file (required)"
    r2_fastq: "Read 2 FASTQ file (optional for paired-end data)"
    cpu_cores: "Number of CPU cores allocated for FastQC"
    memory_gb: "Memory allocated for FastQC in GB"
    adapters: "Optional adapter sequences file for contamination screening"
    limits: "Optional limits file to override default warning/error thresholds"
    contaminants: "Optional contaminants file for contamination screening"
  }

  input {
    File? r1_fastq
    File? r2_fastq
    Int cpu_cores = 2
    Int memory_gb = 4
    File? adapters
    File? limits
    File? contaminants
  }

  command <<<
    set -eo pipefail

    # Create output directory
    mkdir -p fastqc_output

    # Prepare FastQC command with optional parameters
    FASTQC_CMD="fastqc --outdir fastqc_output --threads ~{cpu_cores}"

    # Add optional parameters if provided
    if [ ! -z "~{adapters}" ]; then
        FASTQC_CMD="${FASTQC_CMD} --adapters ~{adapters}"
    fi

    if [ ! -z "~{limits}" ]; then
        FASTQC_CMD="${FASTQC_CMD} --limits ~{limits}"
    fi

    if [ ! -z "~{contaminants}" ]; then
        FASTQC_CMD="${FASTQC_CMD} --contaminants ~{contaminants}"
    fi

    # Run FastQC on input files
    FILES_TO_PROCESS=""
    if [ ! -z "~{r1_fastq}" ]; then
        FILES_TO_PROCESS="${FILES_TO_PROCESS} ~{r1_fastq}"
    fi

    if [ ! -z "~{r2_fastq}" ]; then
        FILES_TO_PROCESS="${FILES_TO_PROCESS} ~{r2_fastq}"
    fi

    if [ -z "${FILES_TO_PROCESS}" ]; then
        echo "Error: No FASTQ files provided"
        exit 1
    fi

    echo "Running FastQC with command: ${FASTQC_CMD} ${FILES_TO_PROCESS}"
    ${FASTQC_CMD} ${FILES_TO_PROCESS}

    # List generated files for verification
    echo "FastQC analysis completed. Generated files:"
    ls -la fastqc_output/
  >>>

  output {
    Array[File] html_reports = glob("fastqc_output/*.html")
    Array[File] zip_reports = glob("fastqc_output/*.zip")
  }

  runtime {
    docker: "getwilds/fastqc:0.12.1"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}