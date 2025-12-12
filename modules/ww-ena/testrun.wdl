version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-ww-ena/modules/ww-ena/ww-ena.wdl" as ww_ena

workflow ena_example {
  meta {
    description: "Test workflow demonstrating the ENA downloader module. Downloads a small test dataset from ENA using accession numbers."
    outputs: {
        downloaded_files: "Array of files downloaded from ENA",
        download_log: "Log file from the ENA download process",
        download_summary: "Summary report of the ENA download process",
        validation_report: "Report validating the success of the ENA download"
    }
  }

  # Test 1: Download a small FASTQ file using accession number
  call ww_ena.download_files as test_download_accession {
    input:
      accessions = "ERR000001",
      file_format = "READS_FASTQ",
      protocol = "FTP",
      output_dir_name = "test_ena_download",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Validate the download
  call validate_download {
    input:
      downloaded_files = test_download_accession.downloaded_files,
      download_log = test_download_accession.download_log,
      download_summary = test_download_accession.download_summary
  }

  output {
    Array[File] downloaded_files = test_download_accession.downloaded_files
    File download_log = test_download_accession.download_log
    File download_summary = test_download_accession.download_summary
    File validation_report = validate_download.report
  }
}

task validate_download {
  meta {
    description: "Validate that ENA download completed successfully"
    outputs: {
        report: "Validation report confirming successful download of files from ENA"
    }
  }

  input {
    Array[File] downloaded_files
    File download_log
    File download_summary
  }

  command <<<
    set -eo pipefail

    echo "=== ENA Download Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    # Check if files were downloaded
    num_files=~{length(downloaded_files)}
    echo "Number of files downloaded: $num_files" >> validation_report.txt
    echo "" >> validation_report.txt

    if [ "$num_files" -gt 0 ]; then
      echo "✓ SUCCESS: Files were downloaded" >> validation_report.txt
    else
      echo "✗ FAILED: No files were downloaded" >> validation_report.txt
      exit 1
    fi

    echo "" >> validation_report.txt
    echo "=== Download Summary ===" >> validation_report.txt
    cat ~{download_summary} >> validation_report.txt

    echo "" >> validation_report.txt
    echo "=== Download Log ===" >> validation_report.txt
    cat ~{download_log} >> validation_report.txt

    echo "" >> validation_report.txt
    echo "=== Overall Result ===" >> validation_report.txt
    echo "✓ All validations passed successfully" >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/ena-tools:2.1.1"
    memory: "2 GB"
    cpu: 1
  }
}
