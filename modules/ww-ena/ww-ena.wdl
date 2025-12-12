## WILDS WDL Module for ENA (European Nucleotide Archive)
## Downloads sequencing data files from ENA using the ena-file-downloader tool.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task download_files {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Downloads sequencing data files from the European Nucleotide Archive (ENA) using the ena-file-downloader tool. Supports download by accession numbers with configurable file formats and transfer protocols."
    outputs: {
        downloaded_files: "Array of downloaded files from ENA",
        download_log: "Log file containing download status and details",
        download_summary: "Summary report of the download operation",
        accessions_used: "The accession numbers that were processed"
    }
  }

  parameter_meta {
    accessions: "Comma-separated list of ENA accession numbers (e.g., 'ERR2208926,ERR2208890') or a single accession"
    accessions_file: "Optional file containing accession numbers (one per line or tab-separated with accessions in first column)"
    file_format: "Format of files to download: READS_FASTQ, READS_SUBMITTED, READS_BAM, ANALYSIS_SUBMITTED, or ANALYSIS_GENERATED"
    protocol: "Transfer protocol to use: FTP or ASPERA (requires aspera_location if ASPERA is selected)"
    aspera_location: "Path to Aspera Connect/CLI installation (required if protocol is ASPERA)"
    output_dir_name: "Name for the output directory where files will be downloaded"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    # Either accessions string or accessions_file must be provided
    String? accessions
    File? accessions_file

    # Download configuration
    String file_format = "READS_FASTQ"
    String protocol = "FTP"
    String? aspera_location

    # Output naming
    String output_dir_name = "ena_downloads"

    # Resource parameters
    Int cpu_cores = 2
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Create output directory
    mkdir -p ~{output_dir_name}

    # Create logs directory
    mkdir -p logs

    # Execute download with ena-file-downloader
    java -jar /usr/local/bin/ena-file-downloader.jar \
      ~{if defined(accessions) then "--accessions " + accessions else ""} \
      ~{if defined(accessions_file) then "--accessions " + accessions_file else ""} \
      --format ~{file_format} \
      --protocol ~{protocol} \
      ~{if defined(aspera_location) then "--asperaLocation " + aspera_location else ""} \
      --location ~{output_dir_name}

    # Find all downloaded files
    find ~{output_dir_name} -type f > downloaded_files.txt

    # Create summary
    echo "Download completed at $(date)" > download_summary.txt
    echo "Number of files downloaded: $(cat downloaded_files.txt | wc -l)" >> download_summary.txt
    echo "Files:" >> download_summary.txt
    cat downloaded_files.txt >> download_summary.txt

    # Copy log file if it exists
    if [ -f "logs/app.log" ]; then
      cp logs/app.log download.log
    else
      echo "No log file generated" > download.log
    fi
  >>>

  output {
    Array[File] downloaded_files = read_lines("downloaded_files.txt")
    File download_log = "download.log"
    File download_summary = "download_summary.txt"
    String accessions_used = select_first([accessions, "from_file"])
  }

  runtime {
    docker: "ena-tools:test"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_by_query {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Downloads sequencing data files from ENA using a search query. Allows filtering by result type and query parameters to retrieve multiple files matching specific criteria."
    outputs: {
        downloaded_files: "Array of downloaded files from ENA",
        download_log: "Log file containing download status and details",
        download_summary: "Summary report of the download operation"
    }
  }

  parameter_meta {
    query: "ENA search query string containing result type and query parameters (e.g., 'result=read_run&query=study_accession=PRJEB1234')"
    file_format: "Format of files to download: READS_FASTQ, READS_SUBMITTED, READS_BAM, ANALYSIS_SUBMITTED, or ANALYSIS_GENERATED"
    protocol: "Transfer protocol to use: FTP or ASPERA (requires aspera_location if ASPERA is selected)"
    aspera_location: "Path to Aspera Connect/CLI installation (required if protocol is ASPERA)"
    output_dir_name: "Name for the output directory where files will be downloaded"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    String query

    # Download configuration
    String file_format = "READS_FASTQ"
    String protocol = "FTP"
    String? aspera_location

    # Output naming
    String output_dir_name = "ena_downloads"

    # Resource parameters
    Int cpu_cores = 2
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Create output directory
    mkdir -p ~{output_dir_name}

    # Create logs directory
    mkdir -p logs

    # Execute download with ena-file-downloader
    java -jar /usr/local/bin/ena-file-downloader.jar \
      --query "~{query}" \
      --format ~{file_format} \
      --protocol ~{protocol} \
      ~{if defined(aspera_location) then "--asperaLocation " + aspera_location else ""} \
      --location ~{output_dir_name}

    # Find all downloaded files
    find ~{output_dir_name} -type f > downloaded_files.txt

    # Create summary
    echo "Download completed at $(date)" > download_summary.txt
    echo "Number of files downloaded: $(cat downloaded_files.txt | wc -l)" >> download_summary.txt
    echo "Files:" >> download_summary.txt
    cat downloaded_files.txt >> download_summary.txt

    # Copy log file if it exists
    if [ -f "logs/app.log" ]; then
      cp logs/app.log download.log
    else
      echo "No log file generated" > download.log
    fi
  >>>

  output {
    Array[File] downloaded_files = read_lines("downloaded_files.txt")
    File download_log = "download.log"
    File download_summary = "download_summary.txt"
  }

  runtime {
    docker: "ena-tools:test"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
