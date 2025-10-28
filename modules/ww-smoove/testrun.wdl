version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-smoove/ww-smoove.wdl" as ww_smoove

struct SmooveSample {
    String name
    File bam
    File bai
}

workflow smoove_example {
  # Download test reference data
  call ww_testdata.download_ref_data { }

  # Download test sample data
  call ww_testdata.download_bam_data { }

  # Create test samples array
  Array[SmooveSample] final_samples = [
    {
      "name": "demo_sample",
      "bam": download_bam_data.bam,
      "bai": download_bam_data.bai
    }
  ]

  scatter (sample in final_samples) {
    call ww_smoove.smoove_call { input:
        aligned_bam = sample.bam,
        aligned_bam_index = sample.bai,
        sample_name = sample.name,
        reference_fasta = download_ref_data.fasta,
        reference_fasta_index = download_ref_data.fasta_index,
        cpu_cores = 2,
        memory_gb = 8
    }
  }

  call validate_outputs { input:
      vcf_files = smoove_call.vcf,
      vcf_index_files = smoove_call.vcf_index
  }

  output {
    Array[File] smoove_vcfs = smoove_call.vcf
    Array[File] smoove_vcf_indexes = smoove_call.vcf_index
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate Smoove outputs and generate a comprehensive report"
    outputs: {
        report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    vcf_files: "Array of VCF files to validate"
    vcf_index_files: "Array of VCF index files to validate"
  }

  input {
    Array[File] vcf_files
    Array[File] vcf_index_files
  }

  command <<<
    set -euo pipefail

    echo "=== Smoove Output Validation Report ===" > validation_report.txt
    echo "Generated: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    # Validate VCF files
    echo "VCF File Validation:" >> validation_report.txt
    vcf_count=0
    for vcf in ~{sep=" " vcf_files}; do
      vcf_count=$((vcf_count + 1))
      echo "  File $vcf_count: $(basename $vcf)" >> validation_report.txt
      
      # Check file exists and is not empty
      if [[ -f "$vcf" && -s "$vcf" ]]; then
        echo "    Status: PASS - File exists and is not empty" >> validation_report.txt
        
        # Count variants
        variant_count=$(zcat "$vcf" | grep -v "^#" | wc -l || echo "0")
        echo "    Variants: $variant_count" >> validation_report.txt
        
        # Check VCF format validity
        if zcat "$vcf" | head -n 50 | grep -q "^#CHROM"; then
          echo "    Format: PASS - Valid VCF header detected" >> validation_report.txt
        else
          echo "    Format: FAIL - Invalid VCF header" >> validation_report.txt
        fi
      else
        echo "    Status: FAIL - File missing or empty" >> validation_report.txt
      fi
      echo "" >> validation_report.txt
    done

    # Validate index files
    echo "VCF Index File Validation:" >> validation_report.txt
    index_count=0
    for index in ~{sep=" " vcf_index_files}; do
      index_count=$((index_count + 1))
      echo "  Index $index_count: $(basename $index)" >> validation_report.txt
      
      if [[ -f "$index" && -s "$index" ]]; then
        echo "    Status: PASS - Index file exists" >> validation_report.txt
      else
        echo "    Status: FAIL - Index file missing" >> validation_report.txt
      fi
      echo "" >> validation_report.txt
    done

    # Summary
    echo "Summary:" >> validation_report.txt
    echo "  Total VCF files processed: $vcf_count" >> validation_report.txt
    echo "  Total index files processed: $index_count" >> validation_report.txt
    echo "  Validation completed: $(date)" >> validation_report.txt

    echo "Validation completed successfully"
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "ubuntu:20.04"
    cpu: 1
    memory: "2GB"
  }
}
