version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/move-consensus/modules/ww-annovar/ww-annovar.wdl" as ww_annovar

workflow annovar_example {
  # Download test data for annotation
  call ww_testdata.download_gnomad_vcf { input:
      region = "chr1:1-1000000",
      filter_name = "chr1_1Mb"
  }

  # Use test VCF data
  Array[File] vcfs_to_process = [download_gnomad_vcf.gnomad_vcf]

  scatter (vcf in vcfs_to_process) {
    call ww_annovar.annovar_annotate { input:
        vcf_to_annotate = vcf,
        ref_name = "hg38",
        annovar_protocols = "refGene,knownGene,cosmic70,esp6500siv2_all,clinvar_20180603,gnomad211_exome",
        annovar_operation = "g,f,f,f,f,f"
    }
  }

  call validate_outputs { input:
      annotated_vcf_files = annovar_annotate.annotated_vcf,
      annotated_table_files = annovar_annotate.annotated_table
  }

  output {
    Array[File] annotated_vcf = annovar_annotate.annotated_vcf
    Array[File] annotated_table = annovar_annotate.annotated_table
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate Annovar outputs and generate summary statistics"
    outputs: {
        report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    annotated_vcf_files: "Array of annotated VCF files to validate"
    annotated_table_files: "Array of annotated table files to validate"
  }

  input {
    Array[File] annotated_vcf_files
    Array[File] annotated_table_files
  }

  command <<<
    set -euo pipefail

    echo "=== ANNOVAR OUTPUT VALIDATION REPORT ===" > validation_report.txt
    echo "Generated: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    # Validate VCF files
    echo "VCF Files Validation:" >> validation_report.txt
    vcf_count=0
    for vcf in ~{sep=" " annotated_vcf_files}; do
      vcf_count=$((vcf_count + 1))
      echo "  File $vcf_count: $(basename $vcf)" >> validation_report.txt

      # Check if file exists and is not empty
      if [ -f "$vcf" ] && [ -s "$vcf" ]; then
        echo "    Status: File exists and is not empty" >> validation_report.txt

        # Count variants
        variant_count=$(zcat "$vcf" | grep -v "^#" | wc -l || echo "0")
        echo "    Variants: $variant_count" >> validation_report.txt

        # Check VCF format
        if zcat "$vcf" | head -1 | grep -q "^##fileformat=VCF"; then
          echo "    Format: Valid VCF header" >> validation_report.txt
        else
          echo "    Format: Invalid VCF header" >> validation_report.txt
        fi
      else
        echo "    Status: File missing or empty" >> validation_report.txt
      fi
      echo "" >> validation_report.txt
    done

    # Validate annotation table files
    echo "Annotation Table Files Validation:" >> validation_report.txt
    table_count=0
    for table in ~{sep=" " annotated_table_files}; do
      table_count=$((table_count + 1))
      echo "  File $table_count: $(basename $table)" >> validation_report.txt

      if [ -f "$table" ] && [ -s "$table" ]; then
        echo "    Status: File exists and is not empty" >> validation_report.txt

        # Count lines (excluding header)
        line_count=$(($(wc -l < "$table") - 1))
        echo "    Annotated variants: $line_count" >> validation_report.txt

        # Check for key annotation columns
        if head -1 "$table" | grep -q "Chr.*Start.*End.*Ref.*Alt"; then
          echo "    Format: Contains expected annotation columns" >> validation_report.txt
        else
          echo "    Format: Missing expected annotation columns" >> validation_report.txt
        fi
      else
        echo "    Status: File missing or empty" >> validation_report.txt
      fi
      echo "" >> validation_report.txt
    done

    echo "=== VALIDATION SUMMARY ===" >> validation_report.txt
    echo "Total VCF files processed: $vcf_count" >> validation_report.txt
    echo "Total annotation tables generated: $table_count" >> validation_report.txt
    echo "Validation completed: $(date)" >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/bcftools:1.19"
    cpu: 1
    memory: "2 GB"
  }
}
