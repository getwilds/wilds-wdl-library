version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-aws-sso/ww-aws-sso.wdl" as ww_aws_sso

workflow aws_sso_example {
  # Test listing bucket contents
  call ww_aws_sso.s3_list_bucket { input:
      s3_uri = "s3://gatk-test-data/wgs_fastq/",
      cpu_cores = 2,
      memory_gb = 4
  }

  # Test file download
  call ww_aws_sso.s3_download_file { input:
      s3_uri = "s3://gatk-test-data/wgs_fastq/NA12878_20k/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq",
      cpu_cores = 2,
      memory_gb = 4
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
