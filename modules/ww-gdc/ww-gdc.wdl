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
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gdc/ww-gdc.wdl"
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
    cpu_cores: "Number of CPU cores allocated for the task (default: 4)"
    memory_gb: "Memory allocated for the task in GB (default: 8)"
  }

  input {
    File manifest_file
    File? token_file
    Int n_processes = 8
    Int retry_amount = 5
    Int wait_time = 5
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Work in /tmp to avoid Cromwell's deep directory structures that cause
    # "AF_UNIX path too long" errors with gdc-client's multiprocessing sockets
    WORK_DIR="/tmp/gdc_$$"
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"

    # Set TMPDIR to short path so Python multiprocessing creates sockets here
    export TMPDIR="$WORK_DIR"
    export TEMP="$WORK_DIR"
    export TMP="$WORK_DIR"

    # Symlink required files to short path
    ln -s ~{manifest_file} manifest.txt

    # Build gdc-client download command
    CMD="gdc-client download \
      --manifest manifest.txt \
      --n-processes ~{n_processes} \
      --retry-amount ~{retry_amount} \
      --wait-time ~{wait_time} \
      --log-file download.log"

    # Add token if provided (for controlled-access data)
    if [ -f "~{token_file}" ]; then
      ln -s ~{token_file} token.txt
      CMD="$CMD --token-file token.txt"
      echo "Using authentication token for controlled-access data" | tee -a download.log
    else
      echo "No token provided - downloading open-access data only" | tee -a download.log
    fi

    # Execute download
    echo "Starting GDC download at $(date)" | tee -a download.log
    echo "Working directory: $WORK_DIR" | tee -a download.log
    echo "TMPDIR: $TMPDIR" | tee -a download.log
    echo "Command: $CMD" | tee -a download.log
    eval $CMD
    echo "Download completed at $(date)" | tee -a download.log

    # Create a file listing all downloaded files
    find . -type f ! -name "download.log" ! -name "manifest_summary.txt" ! -name "manifest.txt" ! -name "token.txt" > manifest_summary.txt
    echo "Downloaded $(wc -l < manifest_summary.txt) files" | tee -a download.log

    # Move results back to original Cromwell execution directory
    EXEC_DIR="$PWD"
    cd -  # Return to previous directory (Cromwell execution dir)

    mv "$EXEC_DIR/download.log" manifest.download.log
    mv "$EXEC_DIR/manifest_summary.txt" manifest.manifest_summary.txt

    # Move downloaded files back, preserving directory structure
    while IFS= read -r file; do
      filepath="$EXEC_DIR/$file"
      if [ -e "$filepath" ]; then
        # Create target directory structure with manifest prefix
        target="manifest/$file"
        target_dir="$(dirname "$target")"
        mkdir -p "$target_dir"
        mv "$filepath" "$target"
      fi
    done < manifest.manifest_summary.txt

    # Update manifest_summary with new paths
    sed -i.bak "s|^|manifest/|" manifest.manifest_summary.txt

    # Clean up
    rm -rf "$EXEC_DIR"
  >>>

  output {
    Array[File] downloaded_files = read_lines("manifest.manifest_summary.txt")
    File download_log = "manifest.download.log"
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
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gdc/ww-gdc.wdl"
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
    cpu_cores: "Number of CPU cores allocated for the task (default: 4)"
    memory_gb: "Memory allocated for the task in GB (default: 8)"
  }

  input {
    Array[String] file_uuids
    File? token_file
    Int n_processes = 8
    Int retry_amount = 5
    Int wait_time = 5
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Work in /tmp to avoid Cromwell's deep directory structures that cause
    # "AF_UNIX path too long" errors with gdc-client's multiprocessing sockets
    WORK_DIR="/tmp/gdc_$$"
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"

    # Set TMPDIR to short path so Python multiprocessing creates sockets here
    export TMPDIR="$WORK_DIR"
    export TEMP="$WORK_DIR"
    export TMP="$WORK_DIR"

    # Write UUIDs to a file for passing to gdc-client
    cat > uuids.txt <<EOL
~{sep='\n' file_uuids}
EOL

    # Build gdc-client download command
    CMD="gdc-client download \
      $(cat uuids.txt | tr '\n' ' ') \
      --n-processes ~{n_processes} \
      --retry-amount ~{retry_amount} \
      --wait-time ~{wait_time} \
      --log-file download.log"

    # Add token if provided (for controlled-access data)
    if [ -f "~{token_file}" ]; then
      ln -s ~{token_file} token.txt
      CMD="$CMD --token-file token.txt"
      echo "Using authentication token for controlled-access data" | tee -a download.log
    else
      echo "No token provided - downloading open-access data only" | tee -a download.log
    fi

    # Execute download
    echo "Starting GDC download at $(date)" | tee -a download.log
    echo "Working directory: $WORK_DIR" | tee -a download.log
    echo "TMPDIR: $TMPDIR" | tee -a download.log
    echo "Downloading $(wc -l < uuids.txt) files by UUID" | tee -a download.log
    echo "Command: $CMD" | tee -a download.log
    eval $CMD
    echo "Download completed at $(date)" | tee -a download.log

    # Create a file listing all downloaded files
    find . -type f ! -name "download.log" ! -name "manifest_summary.txt" ! -name "uuids.txt" ! -name "token.txt" > manifest_summary.txt
    echo "Downloaded $(wc -l < manifest_summary.txt) files" | tee -a download.log

    # Move results back to original Cromwell execution directory
    EXEC_DIR="$PWD"
    cd -  # Return to previous directory (Cromwell execution dir)

    mv "$EXEC_DIR/download.log" uuids.download.log
    mv "$EXEC_DIR/manifest_summary.txt" uuids.manifest_summary.txt

    # Move downloaded files back, preserving directory structure
    while IFS= read -r file; do
      filepath="$EXEC_DIR/$file"
      if [ -e "$filepath" ]; then
        # Create target directory structure with uuids prefix
        target="uuids/$file"
        target_dir="$(dirname "$target")"
        mkdir -p "$target_dir"
        mv "$filepath" "$target"
      fi
    done < uuids.manifest_summary.txt

    # Update manifest_summary with new paths
    sed -i.bak "s|^|uuids/|" uuids.manifest_summary.txt

    # Clean up
    rm -rf "$EXEC_DIR"
  >>>

  output {
    Array[File] downloaded_files = read_lines("uuids.manifest_summary.txt")
    File download_log = "uuids.download.log"
  }

  runtime {
    docker: "getwilds/gdc-client:2.3.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
