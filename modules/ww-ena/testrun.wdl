version 1.0

# import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-ww-ena/modules/ww-ena/ww-ena.wdl" as ww_ena
import "ww-ena.wdl" as ww_ena

workflow ena_example {
  meta {
    description: "Test workflow demonstrating the ENA downloader module. Downloads a small test dataset from ENA using accession numbers and search queries, then extracts FASTQ pairs."
    outputs: {
        downloaded_files: "Array of files downloaded from ENA via accession",
        query_downloaded_files: "Array of files downloaded from ENA via query",
        r1: "Read 1 FASTQ file",
        r2: "Read 2 FASTQ file",
        download_log: "Log file from the ENA download process (accession method)",
        download_summary: "Summary report of the ENA download process (accession method)",
        query_download_log: "Log file from the ENA download process (query method)",
        query_download_summary: "Summary report of the ENA download process (query method)",
        validation_report: "Report validating the success of the ENA download"
    }
  }

  # Test 1: Download a small FASTQ file using accession number
  call ww_ena.download_files as test_download_accession {
    input:
      accessions = "ERR10825982",
      file_format = "READS_FASTQ",
      protocol = "FTP",
      output_dir_name = "test_ena_download",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Test 2: Download the same file using a search query
  # Query format: result=<result_type>&query=<search_criteria>
  # This searches for run_accession matching ERR10825982
  call ww_ena.download_by_query as test_download_query {
    input:
      query = "result=read_run&query=run_accession=\"ERR10825982\"",
      file_format = "READS_FASTQ",
      protocol = "FTP",
      output_dir_name = "test_ena_query_download",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Test 3: Extract FASTQ pairs from downloaded files
  call ww_ena.extract_fastq_pairs {
    input:
      downloaded_files = test_download_accession.downloaded_files
  }

  # Validate the download
  call validate_download {
    input:
      downloaded_files = test_download_accession.downloaded_files,
      download_log = test_download_accession.download_log,
      download_summary = test_download_accession.download_summary,
      query_downloaded_files = test_download_query.downloaded_files,
      query_download_log = test_download_query.download_log,
      query_download_summary = test_download_query.download_summary,
      r1 = extract_fastq_pairs.r1,
      r2 = extract_fastq_pairs.r2
  }

  output {
    Array[File] downloaded_files = test_download_accession.downloaded_files
    Array[File] query_downloaded_files = test_download_query.downloaded_files
    File r1 = extract_fastq_pairs.r1
    File r2 = extract_fastq_pairs.r2
    File download_log = test_download_accession.download_log
    File download_summary = test_download_accession.download_summary
    File query_download_log = test_download_query.download_log
    File query_download_summary = test_download_query.download_summary
    File validation_report = validate_download.report
  }
}

task validate_download {
  meta {
    description: "Validate that ENA download and FASTQ extraction completed successfully"
    outputs: {
        report: "Validation report confirming successful download and extraction of FASTQ files from ENA"
    }
  }

  input {
    Array[File] downloaded_files
    File download_log
    File download_summary
    Array[File] query_downloaded_files
    File query_download_log
    File query_download_summary
    File r1
    File r2
  }

  command <<<
    set -eo pipefail

    echo "=== ENA Download Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    # Check if files were downloaded (accession method)
    echo "=== Accession-Based Download ===" >> validation_report.txt
    num_files=~{length(downloaded_files)}
    echo "Number of files downloaded: $num_files" >> validation_report.txt

    if [ "$num_files" -gt 0 ]; then
      echo "✓ SUCCESS: Files were downloaded via accession" >> validation_report.txt
    else
      echo "✗ FAILED: No files were downloaded via accession" >> validation_report.txt
      exit 1
    fi

    # Check if files were downloaded (query method)
    echo "" >> validation_report.txt
    echo "=== Query-Based Download ===" >> validation_report.txt
    query_num_files=~{length(query_downloaded_files)}
    echo "Number of files downloaded: $query_num_files" >> validation_report.txt

    if [ "$query_num_files" -gt 0 ]; then
      echo "✓ SUCCESS: Files were downloaded via query" >> validation_report.txt
    else
      echo "✗ FAILED: No files were downloaded via query" >> validation_report.txt
      exit 1
    fi

    # Validate R1 and R2 files exist and are non-empty
    echo "" >> validation_report.txt
    echo "=== FASTQ Pair Extraction ===" >> validation_report.txt

    if [ -f "~{r1}" ] && [ -s "~{r1}" ]; then
      echo "✓ SUCCESS: R1 file exists and is non-empty: ~{r1}" >> validation_report.txt
    else
      echo "✗ FAILED: R1 file missing or empty" >> validation_report.txt
      exit 1
    fi

    if [ -f "~{r2}" ] && [ -s "~{r2}" ]; then
      echo "✓ SUCCESS: R2 file exists and is non-empty: ~{r2}" >> validation_report.txt
    else
      echo "✗ FAILED: R2 file missing or empty" >> validation_report.txt
      exit 1
    fi

    echo "" >> validation_report.txt
    echo "=== Accession Download Summary ===" >> validation_report.txt
    cat ~{download_summary} >> validation_report.txt

    echo "" >> validation_report.txt
    echo "=== Accession Download Log ===" >> validation_report.txt
    cat ~{download_log} >> validation_report.txt

    echo "" >> validation_report.txt
    echo "=== Query Download Summary ===" >> validation_report.txt
    cat ~{query_download_summary} >> validation_report.txt

    echo "" >> validation_report.txt
    echo "=== Query Download Log ===" >> validation_report.txt
    cat ~{query_download_log} >> validation_report.txt

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
