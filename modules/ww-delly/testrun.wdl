version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-delly/ww-delly.wdl" as ww_delly

workflow delly_example {
  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_bam_data { }

  # Call structural variants on test sample
  call ww_delly.delly_call { input:
      aligned_bam = download_bam_data.bam,
      aligned_bam_index = download_bam_data.bai,
      reference_fasta = download_ref_data.fasta,
      reference_fasta_index = download_ref_data.fasta_index,
      sv_type = "",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Validate all outputs at the end
  call validate_outputs { input:
      delly_vcfs = [delly_call.vcf],
      delly_vcf_indices = [delly_call.vcf_index]
  }

  output {
    File delly_vcf = delly_call.vcf
    File delly_vcf_index = delly_call.vcf_index
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validates Delly outputs and generates a comprehensive report"
    outputs: {
        report: "Validation report summarizing file checks and variant statistics"
    }
  }

  parameter_meta {
    delly_vcfs: "Array of VCF files from Delly calling"
    delly_vcf_indices: "Array of index files for the VCFs"
  }

  input {
    Array[File] delly_vcfs
    Array[File] delly_vcf_indices
  }

  command <<<
    set -euo pipefail

    {
      echo "=== Delly Workflow Validation Report ==="
      echo "Generated on: $(date)"
      echo "Total samples processed: ~{length(delly_vcfs)}"
      echo ""
    } >> validation_report.txt

    # Validate each VCF file
    VCF_FILES=(~{sep=" " delly_vcfs})
    INDEX_FILES=(~{sep=" " delly_vcf_indices})

    for i in "${!VCF_FILES[@]}"; do
      VCF="${VCF_FILES[$i]}"
      INDEX="${INDEX_FILES[$i]}"
      
      {
        echo "====================="
        echo "- VCF file: $VCF"
        echo "  Size: $(ls -lh "$VCF" | awk '{print $5}')"
        echo "- VCF index: $INDEX"
        echo "  Size: $(ls -lh "$INDEX" | awk '{print $5}')"
      } >> validation_report.txt
      
      # Validate VCF file format
      if bcftools view -h "$VCF" > /dev/null 2>&1; then
        echo "- VCF format: VALID" >> validation_report.txt
      else
        echo "- VCF format: INVALID" >> validation_report.txt
      fi
      
      # Count variants
      VARIANT_COUNT=$(bcftools view -H "$VCF" | wc -l)
      echo "- Total variants: $VARIANT_COUNT" >> validation_report.txt
      
      echo "" >> validation_report.txt
    done

    # Overall summary
    echo "=== Overall Summary ===" >> validation_report.txt
    TOTAL_VARIANTS=0
    for vcf in "${VCF_FILES[@]}"; do
      COUNT=$(bcftools view -H "$vcf" | wc -l)
      TOTAL_VARIANTS=$((TOTAL_VARIANTS + COUNT))
    done
    {
      echo "Total variants across all samples: $TOTAL_VARIANTS"
      echo ""
      echo "Validation completed successfully!"
    } >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/delly:1.2.9"
    cpu: 1
    memory: "4GB"
  }
}
