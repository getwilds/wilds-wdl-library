## WILDS WDL Module for GDC Data Transfer Tool
## This module provides tasks for downloading TCGA and other cancer genomics data
## from the NCI Genomic Data Commons (GDC) using the gdc-client tool.
## The module supports both controlled-access and open-access data downloads.

version 1.0

task download_by_manifest {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Download files from GDC using a manifest file. Supports both controlled-access (with token) and open-access data."
    outputs: {
        downloaded_files: "Array of downloaded data files from GDC",
        download_log: "Log file with download statistics and any errors"
    }
  }

  parameter_meta {
    manifest_file: "GDC manifest file containing file UUIDs and metadata (downloadable from GDC Data Portal)"
    token_file: "Optional authentication token file for controlled-access data (downloadable from GDC user profile)"
    n_processes: "Number of parallel download processes (default: 8)"
    retry_amount: "Number of times to retry failed downloads (default: 5)"
    wait_time: "Seconds to wait between retries (default: 5)"
    output_dir_name: "Name for the output directory where files will be downloaded"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File manifest_file
    File? token_file
    Int n_processes = 8
    Int retry_amount = 5
    Int wait_time = 5
    String output_dir_name = "gdc_download"
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Create output directory
    mkdir -p ~{output_dir_name}

    # Build gdc-client download command
    CMD="gdc-client download \
      --manifest ~{manifest_file} \
      --dir ~{output_dir_name} \
      --n-processes ~{n_processes} \
      --retry-amount ~{retry_amount} \
      --wait-time ~{wait_time} \
      --log-file ~{output_dir_name}/download.log"

    # Add token if provided (for controlled-access data)
    if [ -f "~{token_file}" ]; then
      CMD="$CMD --token-file ~{token_file}"
      echo "Using authentication token for controlled-access data" | tee -a ~{output_dir_name}/download.log
    else
      echo "No token provided - downloading open-access data only" | tee -a ~{output_dir_name}/download.log
    fi

    # Execute download
    echo "Starting GDC download at $(date)" | tee -a ~{output_dir_name}/download.log
    echo "Command: $CMD" | tee -a ~{output_dir_name}/download.log
    eval $CMD
    echo "Download completed at $(date)" | tee -a ~{output_dir_name}/download.log

    # Create a file listing all downloaded files
    find ~{output_dir_name} -type f ! -name "download.log" ! -name "manifest_summary.txt" > ~{output_dir_name}/manifest_summary.txt
    echo "Downloaded $(wc -l < ~{output_dir_name}/manifest_summary.txt) files" | tee -a ~{output_dir_name}/download.log
  >>>

  output {
    Array[File] downloaded_files = read_lines("~{output_dir_name}/manifest_summary.txt")
    File download_log = "~{output_dir_name}/download.log"
  }

  runtime {
    docker: "getwilds/gdc-client:2.3.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_by_uuids {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Download files from GDC using file UUIDs. Supports both controlled-access (with token) and open-access data."
    outputs: {
        downloaded_files: "Array of downloaded data files from GDC",
        download_log: "Log file with download statistics and any errors"
    }
  }

  parameter_meta {
    file_uuids: "Array of GDC file UUIDs to download"
    token_file: "Optional authentication token file for controlled-access data (downloadable from GDC user profile)"
    n_processes: "Number of parallel download processes (default: 8)"
    retry_amount: "Number of times to retry failed downloads (default: 5)"
    wait_time: "Seconds to wait between retries (default: 5)"
    output_dir_name: "Name for the output directory where files will be downloaded"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    Array[String] file_uuids
    File? token_file
    Int n_processes = 8
    Int retry_amount = 5
    Int wait_time = 5
    String output_dir_name = "gdc_download"
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Create output directory
    mkdir -p ~{output_dir_name}

    # Write UUIDs to a file for passing to gdc-client
    cat > uuids.txt <<EOL
~{sep='\n' file_uuids}
EOL

    # Build gdc-client download command
    CMD="gdc-client download \
      $(cat uuids.txt | tr '\n' ' ') \
      --dir ~{output_dir_name} \
      --n-processes ~{n_processes} \
      --retry-amount ~{retry_amount} \
      --wait-time ~{wait_time} \
      --log-file ~{output_dir_name}/download.log"

    # Add token if provided (for controlled-access data)
    if [ -f "~{token_file}" ]; then
      CMD="$CMD --token-file ~{token_file}"
      echo "Using authentication token for controlled-access data" | tee -a ~{output_dir_name}/download.log
    else
      echo "No token provided - downloading open-access data only" | tee -a ~{output_dir_name}/download.log
    fi

    # Execute download
    echo "Starting GDC download at $(date)" | tee -a ~{output_dir_name}/download.log
    echo "Downloading $(wc -l < uuids.txt) files by UUID" | tee -a ~{output_dir_name}/download.log
    echo "Command: $CMD" | tee -a ~{output_dir_name}/download.log
    eval $CMD
    echo "Download completed at $(date)" | tee -a ~{output_dir_name}/download.log

    # Create a file listing all downloaded files
    find ~{output_dir_name} -type f ! -name "download.log" ! -name "manifest_summary.txt" > ~{output_dir_name}/manifest_summary.txt
    echo "Downloaded $(wc -l < ~{output_dir_name}/manifest_summary.txt) files" | tee -a ~{output_dir_name}/download.log
  >>>

  output {
    Array[File] downloaded_files = read_lines("~{output_dir_name}/manifest_summary.txt")
    File download_log = "~{output_dir_name}/download.log"
  }

  runtime {
    docker: "getwilds/gdc-client:2.3.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
