version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-manta/ww-manta.wdl" as ww_manta

workflow manta_example {
  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_bam_data { }

  # Call structural variants on test sample
  call ww_manta.manta_call { input:
      aligned_bam = download_bam_data.bam,
      aligned_bam_index = download_bam_data.bai,
      sample_name = "demo_sample",
      reference_fasta = download_ref_data.fasta,
      reference_fasta_index = download_ref_data.fasta_index,
      is_rna = false,
      cpu_cores = 2,
      memory_gb = 8
  }

  call validate_outputs { input:
      vcf_files = [manta_call.vcf],
      vcf_index_files = [manta_call.vcf_index]
  }

  output {
    File manta_vcf = manta_call.vcf
    File manta_vcf_index = manta_call.vcf_index
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate Manta outputs and generate a comprehensive report"
    outputs: {
        report: "Validation summary with structural variant calling statistics"
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
    set -eo pipefail
    
    {
      echo "Manta Structural Variant Calling Validation Report"
      echo "================================================="
      echo "Generated on: $(date)"
      echo ""
      echo "Sample Summary:"
      echo "Total samples processed: ~{length(vcf_files)}"
      echo ""
    } >> validation_report.txt
    
    # Validate each sample's outputs
    vcf_files=(~{sep=" " vcf_files})
    vcf_index_files=(~{sep=" " vcf_index_files})
    
    for i in "${!vcf_files[@]}"; do
        vcf="${vcf_files[$i]}"
        vcf_index="${vcf_index_files[$i]}"
        
        echo "Sample: $vcf" >> validation_report.txt
        
        # Check if VCF file exists and is not empty
        if [[ -f "$vcf" && -s "$vcf" ]]; then
            echo "VCF file present and non-empty" >> validation_report.txt
            variant_count=$(zcat "$vcf" | grep -c -v '^#' || true)
            echo "Variants called: $variant_count" >> validation_report.txt
        else
            echo "VCF file missing or empty" >> validation_report.txt
        fi
        
        # Check if VCF index exists
        if [[ -f "$vcf_index" ]]; then
            echo "VCF index file present" >> validation_report.txt
        else
            echo "VCF index file missing" >> validation_report.txt
        fi
        
        echo "" >> validation_report.txt
    done
    
    echo "Validation completed successfully" >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/manta:1.6.0"
    memory: "4GB"
    cpu: 1
  }
}
