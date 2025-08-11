## WILDS WDL module for AWS operations using awscli
##
## This module provides reusable tasks for common AWS operations including:
## - Downloading files from S3 buckets (public and private)
## - Uploading files to S3 buckets
## - Syncing directories with S3
## - Listing S3 bucket contents
##
## All tasks use the getwilds/awscli Docker images from the WILDS Docker Library

version 1.0

workflow aws_example {
  meta {
    description: "Demonstration workflow for AWS operations module"
    author: "WILDS WDL Library"
    outputs: {
        bucket_listing: "List of objects in the specified S3 bucket",
        downloaded_test_file: "Test file downloaded from S3",
        validation_report: "Validation report of AWS operations"
    }
  }

  parameter_meta {
    test_bucket_uri: "S3 bucket/prefix to list (optional)"
    aws_config_file: "Path to AWS config file (optional, uses public access if not provided)"
    aws_credentials_file: "Path to AWS credentials file (optional)"
    cpu_cores: "Number of CPU cores for AWS operations"
    memory_gb: "Memory allocation for AWS operations"
  }

  input {
    String? test_bucket_uri
    String? aws_config_file
    String? aws_credentials_file
    Int cpu_cores = 2
    Int memory_gb = 4
  }

  # Test listing bucket contents
  String bucket_to_list = select_first([test_bucket_uri, "s3://gatk-test-data/wgs_fastq/"])

  call s3_list_bucket { input:
      s3_uri = bucket_to_list,
      aws_config_file = aws_config_file,
      aws_credentials_file = aws_credentials_file,
      cpu_cores = cpu_cores,
      memory_gb = memory_gb
  }

  # Test file download
  call s3_download_file { input:
      s3_uri = bucket_to_list + "NA12878_20k/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq",
      aws_config_file = aws_config_file,
      aws_credentials_file = aws_credentials_file,
      cpu_cores = cpu_cores,
      memory_gb = memory_gb
  }

  # Validate all operations
  call validate_outputs { input:
      downloaded_files = [s3_download_file.downloaded_file],
      bucket_listing = s3_list_bucket.file_list,
      object_count = s3_list_bucket.object_count
  }

  output {
    File bucket_listing = s3_list_bucket.file_list
    File downloaded_test_file = s3_download_file.downloaded_file
    File validation_report = validate_outputs.report
  }
}

task s3_download_file {
  meta {
    description: "Download a single file from an S3 bucket"
    outputs: {
        downloaded_file: "Downloaded file from S3"
    }
  }

  parameter_meta {
    s3_uri: "S3 URI of the file to download (e.g., s3://bucket/path/file.txt)"
    output_filename: "Name for the downloaded file (optional, uses original name if not specified)"
    aws_config_file: "Path to AWS config file (optional, uses --no-sign-request if not provided)"
    aws_credentials_file: "Path to AWS credentials file (optional)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    String s3_uri
    String? output_filename
    String? aws_config_file
    String? aws_credentials_file
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  String final_filename = select_first([output_filename, basename(s3_uri)])
  Boolean use_credentials = defined(aws_config_file)
  String sign_request_flag = if use_credentials then "" else "--no-sign-request"

  command <<<
    set -euo pipefail
    
    # Set AWS config if provided
    ~{if defined(aws_config_file) then 'export AWS_CONFIG_FILE="' + aws_config_file + '"' else ''}
    ~{if defined(aws_credentials_file) then 'export AWS_SHARED_CREDENTIALS_FILE="' + aws_credentials_file + '"' else ''}
    
    echo "=== AWS Configuration Debug ==="
    echo "AWS_CONFIG_FILE: ${AWS_CONFIG_FILE:-not set}"
    echo "AWS_SHARED_CREDENTIALS_FILE: ${AWS_SHARED_CREDENTIALS_FILE:-not set}"
    echo "Using credentials: ~{use_credentials}"
    echo "Sign request flag: '~{sign_request_flag}'"
    echo ""
    
    echo "Downloading ~{s3_uri} to ~{final_filename}"
    
    aws s3 cp ~{sign_request_flag} "~{s3_uri}" "~{final_filename}"
    
    # Verify file was downloaded
    if [[ ! -f "~{final_filename}" ]]; then
      echo "ERROR: Failed to download file from ~{s3_uri}"
      exit 1
    fi
    
    echo "Successfully downloaded: $(ls -lh ~{final_filename})"
  >>>

  output {
    File downloaded_file = final_filename
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task s3_upload_file {
  meta {
    description: "Upload a single file to an S3 bucket"
    outputs: {
        s3_uri: "S3 URI of the uploaded file"
    }
  }

  parameter_meta {
    file_to_upload: "File to upload to S3"
    s3_bucket: "S3 bucket name (without s3:// prefix)"
    s3_key: "S3 key/path for the uploaded file (optional, uses filename if not specified)"
    aws_config_file: "Path to AWS config file (required for uploads)"
    aws_credentials_file: "Path to AWS credentials file (optional)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File file_to_upload
    String s3_bucket
    String? s3_key
    String aws_config_file
    String? aws_credentials_file
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  String final_key = select_first([s3_key, basename(file_to_upload)])
  String final_s3_uri = "s3://~{s3_bucket}/~{final_key}"

  command <<<
    set -euo pipefail
    
    # Set AWS config (required for uploads)
    export AWS_CONFIG_FILE="~{aws_config_file}"
    ~{if defined(aws_credentials_file) then 'export AWS_SHARED_CREDENTIALS_FILE="' + aws_credentials_file + '"' else ''}
    
    echo "=== AWS Configuration Debug ==="
    echo "AWS_CONFIG_FILE: $AWS_CONFIG_FILE"
    echo "AWS_SHARED_CREDENTIALS_FILE: ${AWS_SHARED_CREDENTIALS_FILE:-not set}"
    echo ""
    
    echo "Uploading ~{file_to_upload} to ~{final_s3_uri}"
    
    aws s3 cp "~{file_to_upload}" "~{final_s3_uri}"
    
    echo "Successfully uploaded to: ~{final_s3_uri}"
  >>>

  output {
    String s3_uri = final_s3_uri
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task s3_list_bucket {
  meta {
    description: "List contents of an S3 bucket or prefix"
    outputs: {
        file_list: "Text file containing list of S3 objects",
        object_count: "Number of objects found"
    }
  }

  parameter_meta {
    s3_uri: "S3 URI to list (e.g., s3://bucket/ or s3://bucket/prefix/)"
    aws_config_file: "Path to AWS config file (optional, uses --no-sign-request if not provided)"
    aws_credentials_file: "Path to AWS credentials file (optional)"
    recursive: "List recursively (default: true)"
    human_readable: "Use human-readable file sizes (default: true)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    String s3_uri
    String? aws_config_file
    String? aws_credentials_file
    Boolean recursive = true
    Boolean human_readable = true
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  Boolean use_credentials = defined(aws_config_file)
  String sign_request_flag = if use_credentials then "" else "--no-sign-request"
  String recursive_flag = if recursive then "--recursive" else ""
  String human_readable_flag = if human_readable then "--human-readable" else ""

  command <<<
    set -euo pipefail
    
    # Set AWS config if provided
    ~{if defined(aws_config_file) then 'export AWS_CONFIG_FILE="' + aws_config_file + '"' else ''}
    ~{if defined(aws_credentials_file) then 'export AWS_SHARED_CREDENTIALS_FILE="' + aws_credentials_file + '"' else ''}
    
    echo "=== AWS Configuration Debug ==="
    echo "AWS_CONFIG_FILE: ${AWS_CONFIG_FILE:-not set}"
    echo "AWS_SHARED_CREDENTIALS_FILE: ${AWS_SHARED_CREDENTIALS_FILE:-not set}"
    echo "Using credentials: ~{use_credentials}"
    echo "Sign request flag: '~{sign_request_flag}'"
    echo ""

    echo "Listing contents of ~{s3_uri}"
    
    aws s3 ls ~{sign_request_flag} ~{recursive_flag} ~{human_readable_flag} "~{s3_uri}" > s3_list.txt
    
    # Count objects
    object_count=$(wc -l < s3_list.txt)
    echo "Found $object_count objects"
    echo "$object_count" > object_count.txt
    
    # Show first few entries
    echo "First 10 entries:"
    head -10 s3_list.txt
  >>>

  output {
    File file_list = "s3_list.txt"
    Int object_count = read_int("object_count.txt")
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task validate_outputs {
  meta {
    description: "Validate AWS operation outputs and generate summary report"
    outputs: {
        report: "Validation report with AWS operation statistics"
    }
  }

  parameter_meta {
    downloaded_files: "Files that were downloaded from S3 (optional)"
    bucket_listing: "File containing S3 bucket listing (optional)"
    object_count: "Number of objects found in bucket listing (optional)"
  }

  input {
    Array[File]? downloaded_files
    File? bucket_listing
    Int? object_count
  }

  command <<<
    set -euo pipefail
    
    echo "=== AWS Operations Validation Report ===" > validation_report.txt
    echo "Timestamp: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt
    
    # Validate bucket listing
    if [[ -n "~{bucket_listing}" ]] && [[ -f "~{bucket_listing}" ]]; then
      echo "Bucket Listing Validation:" >> validation_report.txt
      listing_lines=$(wc -l < "~{bucket_listing}")
      echo "  [OK] Bucket listing file exists: ~{basename(select_first([bucket_listing, ""]))}" >> validation_report.txt
      echo "  [OK] Listing contains $listing_lines lines" >> validation_report.txt
      
      if [[ -n "~{object_count}" ]]; then
        reported_count=~{select_first([object_count, 0])}
        if [[ "$listing_lines" -eq "$reported_count" ]]; then
          echo "  [OK] Object count matches: $reported_count objects" >> validation_report.txt
        else
          echo "  [WARNING] Object count mismatch: reported $reported_count, found $listing_lines" >> validation_report.txt
        fi
      fi
      
      echo "  First 5 entries from listing:" >> validation_report.txt
      head -5 "~{bucket_listing}" | sed 's/^/    /' >> validation_report.txt
    fi
    
    # Validate downloaded files
    download_array=~{write_lines(select_first([downloaded_files, []]))}
    if [[ -s "$download_array" ]]; then
      echo "Downloaded Files Validation:" >> validation_report.txt
      download_count=0
      total_size=0
      
      while IFS= read -r file; do
        if [[ -f "$file" ]]; then
          size=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null || echo "0")
          total_size=$((total_size + size))
          download_count=$((download_count + 1))
          echo "  [OK] $(basename "$file"): $size bytes" >> validation_report.txt
        else
          echo "  [ERROR] Missing: $file" >> validation_report.txt
        fi
      done < "$download_array"
      
      echo "  Total downloaded files: $download_count" >> validation_report.txt
      echo "  Total downloaded size: $total_size bytes" >> validation_report.txt
    else
      echo "Downloaded Files Validation: No files to validate" >> validation_report.txt
      echo "" >> validation_report.txt
    fi
    
    echo "All AWS operations completed successfully" >> validation_report.txt
    echo "=== Validation Complete ===" >> validation_report.txt
    
    # Display report
    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: 1
    memory: "2 GB"
  }
}
