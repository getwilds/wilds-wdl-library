version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-annotsv/ww-annotsv.wdl" as ww_annotsv

workflow annotsv_example {
  # Download test data for annotation
  call ww_testdata.download_annotsv_vcf { }

  # Use test VCF data
  Array[File] vcfs_to_process = [download_annotsv_vcf.test_vcf]

  scatter (vcf in vcfs_to_process) {
    call ww_annotsv.annotsv_annotate { input:
        raw_vcf = vcf,
        genome_build = "GRCh38",
        sv_min_size = 50,
        annotation_mode = "full",
        include_ci = true,
        overlap_threshold = 70,
        exclude_benign = false,
        cpu_cores = 4,
        memory_gb = 8
    }
  }

  call validate_outputs { input:
      annotated_tsv_files = annotsv_annotate.annotated_tsv
  }

  output {
    Array[File] annotated_tsv = annotsv_annotate.annotated_tsv
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate AnnotSV outputs and generate comprehensive statistics"
    outputs: {
        report: "Validation summary with structural variant annotation statistics"
    }
  }

  parameter_meta {
    annotated_tsv_files: "Array of annotated TSV files to validate"
  }

  input {
    Array[File] annotated_tsv_files
  }

  command <<<
    set -eo pipefail
    
    echo "AnnotSV Validation Report" > validation_report.txt
    echo "=========================" >> validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt
    
    # Validate TSV files
    echo "TSV File Validation:" >> validation_report.txt
    echo "-------------------" >> validation_report.txt
    
    TSV_COUNT=0
    TOTAL_ANNOTATIONS=0
    for tsv_file in ~{sep=" " annotated_tsv_files}; do
      TSV_COUNT=$((TSV_COUNT + 1))
      echo "TSV $TSV_COUNT: $(basename $tsv_file)" >> validation_report.txt
      
      # Count annotations in the TSV file
      ANNOTATION_COUNT=$(wc -l < "$tsv_file")
      TOTAL_ANNOTATIONS=$((TOTAL_ANNOTATIONS + ANNOTATION_COUNT))

      # Check if file exists and is not empty
      if [ -f "$tsv_file" ] && [ -s "$tsv_file" ]; then
        echo "  File exists and is not empty" >> validation_report.txt
      else
        echo "  File missing or empty" >> validation_report.txt
      fi
    done
    
    echo "" >> validation_report.txt
    echo "Summary Statistics:" >> validation_report.txt
    echo "------------------" >> validation_report.txt
    echo "Total TSV files processed: $TSV_COUNT" >> validation_report.txt
    echo "Total annotations generated: $TOTAL_ANNOTATIONS" >> validation_report.txt
    
    echo "" >> validation_report.txt
    echo "Validation completed successfully." >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/annotsv:3.4.4"
    memory: "4GB"
    cpu: 1
  }
}
