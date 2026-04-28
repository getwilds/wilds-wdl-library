## WILDS WDL module for AWS operations using traditional SSO credentials via awscli
##
## This module provides reusable tasks for common AWS operations including:
## - Downloading files from S3 buckets (public and private)
## - Uploading files to S3 buckets
## - Syncing directories with S3
## - Listing S3 bucket contents
##
## All tasks use the getwilds/awscli Docker images from the WILDS Docker Library

version 1.0

task s3_download_file {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Download a single file from an S3 bucket"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-aws-sso/ww-aws-sso.wdl"
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
    File? aws_config_file
    File? aws_credentials_file
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Upload a single file to an S3 bucket"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-aws-sso/ww-aws-sso.wdl"
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
    File aws_config_file
    File? aws_credentials_file
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "List contents of an S3 bucket or prefix"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-aws-sso/ww-aws-sso.wdl"
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
    File? aws_config_file
    File? aws_credentials_file
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
