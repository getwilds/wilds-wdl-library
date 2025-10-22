version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bedtools/ww-bedtools.wdl" as ww_bedtools

workflow bedtools_example {
  # Hard-coded configuration for testing
  Array[String] chromosomes = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
  ]
  String intersect_flags = "-header -wo"
  String tmp_dir = "/tmp"
  Int cpus = 2
  Int memory_gb = 8

  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_bam_data { }

  # Run bedtools analysis on test sample
  call ww_bedtools.coverage { input:
      bed_file = download_ref_data.bed,
      aligned_bam = download_bam_data.bam,
      sample_name = "demo_sample",
      cpu_cores = cpus,
      memory_gb = memory_gb
  }

  call ww_bedtools.intersect { input:
      bed_file = download_ref_data.bed,
      aligned_bam = download_bam_data.bam,
      sample_name = "demo_sample",
      flags = intersect_flags,
      cpu_cores = cpus,
      memory_gb = memory_gb
  }

  call ww_bedtools.makewindows { input:
      bed_file = download_ref_data.bed,
      aligned_bam = download_bam_data.bam,
      bam_index = download_bam_data.bai,
      sample_name = "demo_sample",
      reference_fasta = download_ref_data.fasta,
      reference_index = download_ref_data.fasta_index,
      list_chr = chromosomes,
      tmp_dir = tmp_dir,
      cpu_cores = cpus,
      memory_gb = memory_gb
  }

  call validate_outputs { input:
      intersect_files = [intersect.intersect_output],
      coverage_files = [coverage.mean_coverage],
      window_count_files = [makewindows.counts_bed],
      sample_names = [coverage.name]
  }

  output {
    File intersect_result = intersect.intersect_output
    File coverage_result = coverage.mean_coverage
    File window_count = makewindows.counts_bed
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  # TODO: Do a basic check of the file contents too
  meta {
    description: "Validate that BEDTools output files exist and are non-empty"
    outputs: {
        report: "Validation report summarizing file check results"
    }
  }

  parameter_meta {
    intersect_files: "BEDTools intersect output files for each sample"
    coverage_files: "Coverage analysis results across BED intervals"
    window_count_files: "Tarballs of per-chromosome BED files of read counts"
    sample_names: "Array of sample names that were processed"
  }

  input {
    Array[File] intersect_files
    Array[File] coverage_files
    Array[File] window_count_files
    Array[String] sample_names
  }

  command <<<
    set -eo pipefail

    echo "=== BEDTools Workflow Validation Report ===" > bedtools_validation_report.txt
    echo "" >> bedtools_validation_report.txt

    intersect_files=("~{sep=" " intersect_files}")
    coverage_files=("~{sep=" " coverage_files}")
    window_count_files=("~{sep=" " window_count_files}")
    sample_names=("~{sep=" " sample_names}")

    validation_passed=true

    for i in "${!sample_names[@]}"; do
      sample="${sample_names[$i]}"
      intersect="${intersect_files[$i]}"
      coverage="${coverage_files[$i]}"
      tarball="${window_count_files[$i]}"

      echo "--- Sample: $sample ---" >> bedtools_validation_report.txt

      # Intersect file
      if [[ -f "$intersect" && -s "$intersect" ]]; then
        size=$(stat -c%s "$intersect")
        lines=$(wc -l < "$intersect")
        echo "  Intersect: PASS (${size} bytes, ${lines} lines)" >> bedtools_validation_report.txt
      else
        echo "  Intersect: FAIL - MISSING OR EMPTY" >> bedtools_validation_report.txt
        validation_passed=false
      fi

      # Coverage file
      if [[ -f "$coverage" && -s "$coverage" ]]; then
        size=$(stat -c%s "$coverage")
        lines=$(wc -l < "$coverage")
        echo "  Coverage: PASS (${size} bytes, ${lines} lines)" >> bedtools_validation_report.txt
      else
        echo "  Coverage: FAIL - MISSING OR EMPTY" >> bedtools_validation_report.txt
        validation_passed=false
      fi

      # Window count tarball
      if [[ -f "$tarball" && -s "$tarball" ]]; then
        size=$(stat -c%s "$tarball")
        echo "  Windows BED file Tarball: PASS (${size} bytes)" >> bedtools_validation_report.txt
      else
        echo "  Window BED file Tarball: FAIL - MISSING OR EMPTY" >> bedtools_validation_report.txt
        validation_passed=false
      fi

      echo "" >> bedtools_validation_report.txt
    done

    echo "=== Validation Summary ===" >> bedtools_validation_report.txt
    echo "Total samples processed: ${#sample_names[@]}" >> bedtools_validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> bedtools_validation_report.txt
    else
      echo "Overall Status: FAILED" >> bedtools_validation_report.txt
      exit 1
    fi

    # Also output to stdout for immediate feedback
    cat bedtools_validation_report.txt
  >>>

  output {
    File report = "bedtools_validation_report.txt"
  }

  runtime {
    docker: "getwilds/bedtools:2.31.1"
    cpu: 1
    memory: "2 GB"
  }
}
