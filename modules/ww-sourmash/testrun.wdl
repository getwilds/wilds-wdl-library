version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-sourmash/modules/ww-sourmash/ww-sourmash.wdl" as ww_sourmash

workflow sourmash_example {
  # Pull down reference genome and index files for chr1
  call ww_testdata.download_ref_data { input:
      chromo = "chr1",
      version = "hg38",
      region = "1-10000000"
  }

  call ww_testdata.download_bam_data { }

  # Create sketch from NA12878 BAM
  call ww_sourmash.sketch as sketch_from_bam { input:
      infile = download_bam_data.bam,
      bam_as_input = true,
      k_value = 31,
      scaled = 1000,
      track_abundance = true,
      cpu_cores = 2,
      output_name = "NA12878"
  }

  # Create sketch from chr1 FASTA
  call ww_sourmash.sketch as sketch_from_fasta { input:
      infile = download_ref_data.fasta,
      bam_as_input = false,
      k_value = 31,
      scaled = 1000,
      track_abundance = true,
      cpu_cores = 2,
      output_name = "chr1"
  }

  # Run gather
  call ww_sourmash.gather as gather_sample1 { input:
      query_sig = sketch_from_bam.sig,
      reference_databases_sigs = [sketch_from_fasta.sig],
      output_name = "gather"
  }

  # Compare the two signatures
  call ww_sourmash.compare as compare_samples { input:
      sig_inputs = [sketch_from_bam.sig, sketch_from_fasta.sig],
      save_name = "compare",
      k_value = 31
  }

  call validate_outputs { input:
      sig_files = [sketch_from_bam.sig, sketch_from_fasta.sig],
      gather_results = [gather_sample1.result],
      compare_npy = compare_samples.npy,
      compare_csv = compare_samples.csv
  }

  output {
    File sample1_sig = sketch_from_bam.sig
    File sample2_sig = sketch_from_fasta.sig
    File sample1_gather_result = gather_sample1.result
    File comparison_matrix_npy = compare_samples.npy
    File comparison_matrix_csv = compare_samples.csv
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate all sourmash output files"
    outputs: {
        report: "Validation report summarizing file checks and statistics"
    }
  }

  parameter_meta {
    sig_files: "Array of signature files to validate"
    gather_results: "Array of gather result CSV files to validate"
    compare_npy: "Comparison matrix numpy file to validate"
    compare_csv: "Comparison matrix CSV file to validate"
  }

  input {
    Array[File] sig_files
    Array[File] gather_results
    File compare_npy
    File compare_csv
  }

  command <<<
    set -eo pipefail

    echo "=== Sourmash Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Validate signature files
    echo "--- Sourmash Sketch Signature Files ---" >> validation_report.txt
    for sig_file in ~{sep=" " sig_files}; do
      if [[ -f "$sig_file" && -s "$sig_file" ]]; then
        sig_size=$(wc -c < "$sig_file")
        echo "Signature file: $sig_file (${sig_size} bytes)" >> validation_report.txt
      else
        echo "Signature file: $sig_file - MISSING OR EMPTY" >> validation_report.txt
        validation_passed=false
      fi
    done
    echo "" >> validation_report.txt

    # Validate gather results
    echo "--- Sourmash Gather Results ---" >> validation_report.txt
    for gather_file in ~{sep=" " gather_results}; do
      if [[ -f "$gather_file" && -s "$gather_file" ]]; then
        gather_size=$(wc -c < "$gather_file")
        echo "Gather result: $gather_file (${gather_size} bytes)" >> validation_report.txt
      else
        echo "Gather result: $gather_file - MISSING OR EMPTY" >> validation_report.txt
        validation_passed=false
      fi
    done
    echo "" >> validation_report.txt

    # Validate compare outputs
    echo "--- Sourmash Compare Results ---" >> validation_report.txt
    if [[ -f "~{compare_npy}" && -s "~{compare_npy}" ]]; then
      npy_size=$(wc -c < "~{compare_npy}")
      echo "Compare matrix (NPY): ~{compare_npy} (${npy_size} bytes)" >> validation_report.txt
    else
      echo "Compare matrix (NPY): ~{compare_npy} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{compare_csv}" && -s "~{compare_csv}" ]]; then
      csv_size=$(wc -c < "~{compare_csv}")
      echo "Compare matrix (CSV): ~{compare_csv} (${csv_size} bytes)" >> validation_report.txt
    else
      echo "Compare matrix (CSV): ~{compare_csv} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Overall summary
    echo "=== Validation Summary ===" >> validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> validation_report.txt
    else
      echo "Overall Status: FAILED" >> validation_report.txt
      exit 1
    fi

    # Also output to stdout for immediate feedback
    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/sourmash:4.8.2"
    memory: "2 GB"
    cpu: 1
  }
}
