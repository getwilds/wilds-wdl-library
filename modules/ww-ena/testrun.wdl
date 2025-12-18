version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-ena/ww-ena.wdl" as ww_ena

workflow ena_example {
  meta {
    description: "Test workflow demonstrating the ENA downloader module. Downloads test datasets from ENA using single accession, multiple accessions, and search queries, then extracts FASTQ pairs."
    outputs: {
        single_downloaded_files: "Array of files downloaded from ENA via single accession",
        single_accessions_used: "Accession string used for single download",
        single_r1_files: "Array of Read 1 FASTQ files from single accession download",
        single_r2_files: "Array of Read 2 FASTQ files from single accession download",
        single_extracted_accessions: "Array of accession IDs extracted from single download filenames",
        multi_downloaded_files: "Array of files downloaded from ENA via multiple accessions",
        multi_accessions_used: "Accession string used for multiple download",
        multi_r1_files: "Array of Read 1 FASTQ files from multiple accession download",
        multi_r2_files: "Array of Read 2 FASTQ files from multiple accession download",
        multi_extracted_accessions: "Array of accession IDs extracted from multiple download filenames",
        query_downloaded_files: "Array of files downloaded from ENA via query",
        validation_report: "Report validating the success of the ENA download"
    }
  }

  # Test 1: Download a single accession
  call ww_ena.download_files as test_single_accession {
    input:
      accessions = "ERR10825982",
      file_format = "READS_FASTQ",
      protocol = "FTP",
      output_dir_name = "test_single_accession",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Test 2: Download multiple accessions
  call ww_ena.download_files as test_multi_accession {
    input:
      accessions = "ERR10825982,ERR10825750",
      file_format = "READS_FASTQ",
      protocol = "FTP",
      output_dir_name = "test_multi_accession",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Test 3: Download using a search query
  # Query format: result=<result_type>&query=<search_criteria>
  call ww_ena.download_by_query as test_download_query {
    input:
      query = "result=read_run&query=run_accession=\"ERR10825982\"",
      file_format = "READS_FASTQ",
      protocol = "FTP",
      output_dir_name = "test_query_download",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Test 4: Extract FASTQ pairs from single accession download
  call ww_ena.extract_fastq_pairs as extract_single {
    input:
      downloaded_files = test_single_accession.downloaded_files
  }

  # Test 5: Extract FASTQ pairs from multiple accession download
  call ww_ena.extract_fastq_pairs as extract_multi {
    input:
      downloaded_files = test_multi_accession.downloaded_files
  }

  # Validate all downloads
  call validate_download {
    input:
      single_accessions_used = test_single_accession.accessions_used,
      single_downloaded_files = test_single_accession.downloaded_files,
      single_r1_files = extract_single.r1_files,
      single_r2_files = extract_single.r2_files,
      single_extracted_accessions = extract_single.accessions,
      multi_accessions_used = test_multi_accession.accessions_used,
      multi_downloaded_files = test_multi_accession.downloaded_files,
      multi_r1_files = extract_multi.r1_files,
      multi_r2_files = extract_multi.r2_files,
      multi_extracted_accessions = extract_multi.accessions,
      query_downloaded_files = test_download_query.downloaded_files
  }

  output {
    # Single accession outputs
    Array[File] single_downloaded_files = test_single_accession.downloaded_files
    String single_accessions_used = test_single_accession.accessions_used
    Array[File] single_r1_files = extract_single.r1_files
    Array[File] single_r2_files = extract_single.r2_files
    Array[String] single_extracted_accessions = extract_single.accessions

    # Multiple accession outputs
    Array[File] multi_downloaded_files = test_multi_accession.downloaded_files
    String multi_accessions_used = test_multi_accession.accessions_used
    Array[File] multi_r1_files = extract_multi.r1_files
    Array[File] multi_r2_files = extract_multi.r2_files
    Array[String] multi_extracted_accessions = extract_multi.accessions

    # Query outputs
    Array[File] query_downloaded_files = test_download_query.downloaded_files

    # Validation
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
    # Single accession inputs
    String single_accessions_used
    Array[File] single_downloaded_files
    Array[File] single_r1_files
    Array[File] single_r2_files
    Array[String] single_extracted_accessions

    # Multiple accession inputs
    String multi_accessions_used
    Array[File] multi_downloaded_files
    Array[File] multi_r1_files
    Array[File] multi_r2_files
    Array[String] multi_extracted_accessions

    # Query inputs
    Array[File] query_downloaded_files
  }

  command <<<
    set -eo pipefail

    echo "=== ENA Download Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    #####################################
    # Test 1: Single Accession Download
    #####################################
    echo "=== Test 1: Single Accession Download ===" >> validation_report.txt
    echo "Accessions used: ~{single_accessions_used}" >> validation_report.txt
    single_num_files=~{length(single_downloaded_files)}
    echo "Number of files downloaded: $single_num_files" >> validation_report.txt

    if [ "$single_num_files" -gt 0 ]; then
      echo "✓ SUCCESS: Files were downloaded via single accession" >> validation_report.txt
    else
      echo "✗ FAILED: No files were downloaded via single accession" >> validation_report.txt
      exit 1
    fi

    # Validate single accession FASTQ extraction
    echo "" >> validation_report.txt
    echo "--- Single Accession FASTQ Extraction ---" >> validation_report.txt
    single_num_r1=~{length(single_r1_files)}
    single_num_r2=~{length(single_r2_files)}
    single_num_acc=~{length(single_extracted_accessions)}
    echo "Number of R1 files: $single_num_r1" >> validation_report.txt
    echo "Number of R2 files: $single_num_r2" >> validation_report.txt
    echo "Extracted accessions: ~{sep=',' single_extracted_accessions}" >> validation_report.txt

    # Verify single accession should produce exactly 1 pair
    if [ "$single_num_r1" -eq 1 ] && [ "$single_num_r2" -eq 1 ]; then
      echo "✓ SUCCESS: Single accession produced exactly 1 FASTQ pair" >> validation_report.txt
    else
      echo "✗ FAILED: Single accession should produce exactly 1 pair, got $single_num_r1 R1 and $single_num_r2 R2" >> validation_report.txt
      exit 1
    fi

    # Verify extracted accession matches input
    if [[ "~{single_accessions_used}" == *"~{single_extracted_accessions[0]}"* ]]; then
      echo "✓ SUCCESS: Extracted accession matches input" >> validation_report.txt
    else
      echo "✗ FAILED: Extracted accession (~{single_extracted_accessions[0]}) doesn't match input (~{single_accessions_used})" >> validation_report.txt
      exit 1
    fi

    #####################################
    # Test 2: Multiple Accession Download
    #####################################
    echo "" >> validation_report.txt
    echo "=== Test 2: Multiple Accession Download ===" >> validation_report.txt
    echo "Accessions used: ~{multi_accessions_used}" >> validation_report.txt
    multi_num_files=~{length(multi_downloaded_files)}
    echo "Number of files downloaded: $multi_num_files" >> validation_report.txt

    if [ "$multi_num_files" -gt 0 ]; then
      echo "✓ SUCCESS: Files were downloaded via multiple accessions" >> validation_report.txt
    else
      echo "✗ FAILED: No files were downloaded via multiple accessions" >> validation_report.txt
      exit 1
    fi

    # Validate multiple accession FASTQ extraction
    echo "" >> validation_report.txt
    echo "--- Multiple Accession FASTQ Extraction ---" >> validation_report.txt
    multi_num_r1=~{length(multi_r1_files)}
    multi_num_r2=~{length(multi_r2_files)}
    multi_num_acc=~{length(multi_extracted_accessions)}
    echo "Number of R1 files: $multi_num_r1" >> validation_report.txt
    echo "Number of R2 files: $multi_num_r2" >> validation_report.txt
    echo "Extracted accessions: ~{sep=',' multi_extracted_accessions}" >> validation_report.txt

    # Verify multiple accessions should produce more than 1 pair
    if [ "$multi_num_r1" -gt 1 ] && [ "$multi_num_r1" -eq "$multi_num_r2" ]; then
      echo "✓ SUCCESS: Multiple accessions produced $multi_num_r1 FASTQ pairs" >> validation_report.txt
    else
      echo "✗ FAILED: Multiple accessions should produce >1 pairs, got $multi_num_r1 R1 and $multi_num_r2 R2" >> validation_report.txt
      exit 1
    fi

    # Verify all extracted accessions are in the input string
    echo "" >> validation_report.txt
    echo "--- Verifying extracted accessions match input ---" >> validation_report.txt
    for acc in ~{sep=' ' multi_extracted_accessions}; do
      if [[ "~{multi_accessions_used}" == *"$acc"* ]]; then
        echo "✓ Accession $acc found in input" >> validation_report.txt
      else
        echo "✗ FAILED: Extracted accession $acc not found in input (~{multi_accessions_used})" >> validation_report.txt
        exit 1
      fi
    done

    #####################################
    # Test 3: Query-Based Download
    #####################################
    echo "" >> validation_report.txt
    echo "=== Test 3: Query-Based Download ===" >> validation_report.txt
    query_num_files=~{length(query_downloaded_files)}
    echo "Number of files downloaded: $query_num_files" >> validation_report.txt

    if [ "$query_num_files" -gt 0 ]; then
      echo "✓ SUCCESS: Files were downloaded via query" >> validation_report.txt
    else
      echo "✗ FAILED: No files were downloaded via query" >> validation_report.txt
      exit 1
    fi

    #####################################
    # Validate all files exist and are non-empty
    #####################################
    echo "" >> validation_report.txt
    echo "=== File Integrity Checks ===" >> validation_report.txt

    for r1_file in ~{sep=' ' single_r1_files} ~{sep=' ' multi_r1_files}; do
      if [ -f "$r1_file" ] && [ -s "$r1_file" ]; then
        echo "✓ R1 file OK: $r1_file" >> validation_report.txt
      else
        echo "✗ FAILED: R1 file missing or empty: $r1_file" >> validation_report.txt
        exit 1
      fi
    done

    for r2_file in ~{sep=' ' single_r2_files} ~{sep=' ' multi_r2_files}; do
      if [ -f "$r2_file" ] && [ -s "$r2_file" ]; then
        echo "✓ R2 file OK: $r2_file" >> validation_report.txt
      else
        echo "✗ FAILED: R2 file missing or empty: $r2_file" >> validation_report.txt
        exit 1
      fi
    done

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
