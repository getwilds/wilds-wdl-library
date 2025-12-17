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

  # Resolve accessions input - prefer string accessions over file
  String accessions_arg = select_first([accessions, accessions_file])

  command <<<
    set -eo pipefail

    # Execute download with ena-file-downloader
    java -jar /usr/local/bin/ena-file-downloader.jar \
      --accessions=~{accessions_arg} \
      --format=~{file_format} \
      --protocol=~{protocol} \
      ~{if defined(aspera_location) then "--asperaLocation=" + aspera_location else ""} \
      --location=~{output_dir_name}

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
    docker: "getwilds/ena-tools:2.1.1"
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

    # Execute download with ena-file-downloader
    java -jar /usr/local/bin/ena-file-downloader.jar \
      --query="~{query}" \
      --format=~{file_format} \
      --protocol=~{protocol} \
      ~{if defined(aspera_location) then "--asperaLocation=" + aspera_location else ""} \
      --location=~{output_dir_name}

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
    docker: "getwilds/ena-tools:2.1.1"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task extract_fastq_pairs {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Extract R1 and R2 FASTQ files from ENA downloads for downstream processing. This task identifies paired-end FASTQ files from the downloaded files array and extracts the accession ID from the filename."
    outputs: {
        r1: "Read 1 FASTQ file",
        r2: "Read 2 FASTQ file",
        accession: "ENA accession ID extracted from the filename"
    }
  }

  parameter_meta {
    downloaded_files: "Array of files downloaded from ENA"
  }

  input {
    Array[File] downloaded_files
  }

  command <<<
    set -eo pipefail

    # List all downloaded files (using actual localized paths)
    echo "Downloaded files:"
    ls -lh ~{sep=' ' downloaded_files}

    # Find R1 and R2 files (ENA typically names them with _1 and _2 or _R1 and _R2)
    # Look for patterns: *_1.fastq.gz, *_R1.fastq.gz, *_1.fq.gz, *_R1.fq.gz, etc.
    R1_FILE=$(ls ~{sep=' ' downloaded_files} | grep -E "(_1\.fastq|_R1\.fastq|_1\.fq|_R1\.fq)" | head -1 || echo "")
    R2_FILE=$(ls ~{sep=' ' downloaded_files} | grep -E "(_2\.fastq|_R2\.fastq|_2\.fq|_R2\.fq)" | head -1 || echo "")

    if [ -z "$R1_FILE" ] || [ -z "$R2_FILE" ]; then
      echo "ERROR: Could not identify paired FASTQ files"
      echo "Looking for files matching pattern *_1.fastq.gz and *_2.fastq.gz"
      echo "Available files:"
      ls -lh ~{sep=' ' downloaded_files}
      exit 1
    fi

    echo "Identified R1: $R1_FILE"
    echo "Identified R2: $R2_FILE"

    # Extract accession ID from filename (everything before _1 or _R1)
    BASENAME=$(basename "$R1_FILE")
    ACCESSION=$(echo "$BASENAME" | sed -E 's/(_1\.fastq|_R1\.fastq|_1\.fq|_R1\.fq).*//')
    echo "Extracted accession: $ACCESSION"
    echo "$ACCESSION" > accession.txt

    # Copy files to execution directory with standardized names
    cp "$R1_FILE" r1.fastq.gz
    cp "$R2_FILE" r2.fastq.gz
  >>>

  output {
    File r1 = "r1.fastq.gz"
    File r2 = "r2.fastq.gz"
    String accession = read_string("accession.txt")
  }

  runtime {
    docker: "getwilds/ena-tools:2.1.1"
    cpu: 1
    memory: "2 GB"
  }
}
