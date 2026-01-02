version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-diamond/ww-diamond.wdl" as ww_diamond
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow diamond_example {
  call ww_testdata.download_diamond_data { }

  # Build DIAMOND database from reference proteome
  call ww_diamond.make_database { input:
    fasta = download_diamond_data.reference,
    memory_gb = 4,
    cpu_cores = 2
  }

  # Run DIAMOND BLASTP alignment
  call ww_diamond.diamond_blastp { input:
    diamond_db = make_database.diamond_db,
    query = download_diamond_data.query,
    memory_gb = 4,
    cpu_cores = 2
  }

  # Validate outputs
  call validate_outputs { input:
    database = make_database.diamond_db,
    alignment = diamond_blastp.aln
  }

  output {
    File diamond_database = make_database.diamond_db
    File alignment_results = diamond_blastp.aln
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that DIAMOND output files were generated correctly"
    outputs: {
        report: "Validation report summarizing file checks"
    }
  }

  parameter_meta {
    database: "DIAMOND database file to validate"
    alignment: "DIAMOND alignment file to validate"
    cpu_cores: "Number of CPU cores to use for validation"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    File database
    File alignment
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -euo pipefail

    # Function to validate that a file exists and is non-empty
    # Function to validate a file exists and is non-empty
    validate_file() {
      local file_path="$1"
      local file_label="$2"

      if [[ -f "$file_path" && -s "$file_path" ]]; then
        echo "$file_label: $file_path - PASSED" >> validation_report.txt
      else
        echo "$file_label: $file_path - MISSING OR EMPTY" >> validation_report.txt
      fi
    }

    echo "=== DIAMOND Alignment Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    validate_file "~{database}" "DIAMOND database" || validation_passed=false
    validate_file "~{alignment}" "DIAMOND alignment" || validation_passed=false

    {
      echo ""
      echo "=== Validation Summary ==="
      echo "Total files validated: 2"
    } >> validation_report.txt

    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> validation_report.txt
    else
      echo "Overall Status: FAILED" >> validation_report.txt
      exit 1
    fi

    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
