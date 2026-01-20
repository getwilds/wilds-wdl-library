version 1.0

import "ww-bcftools.wdl" as ww_bcftools
import "../../modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow bcftools_example {
  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_bam_data { }

  # Call variants on test sample
  call ww_bcftools.mpileup_call { input:
      bam_file = download_bam_data.bam,
      bam_index = download_bam_data.bai,
      reference_fasta = download_ref_data.fasta,
      reference_fasta_index = download_ref_data.fasta_index,
      max_depth = 10000,
      max_idepth = 10000,
      memory_gb = 8,
      cpu_cores = 2
  }

  # Test concat task by concatenating the same VCF twice (with allow_overlaps)
  call ww_bcftools.concat { input:
      vcf_files = [mpileup_call.mpileup_vcf, mpileup_call.mpileup_vcf],
      vcf_indices = [mpileup_call.mpileup_vcf_index, mpileup_call.mpileup_vcf_index],
      output_prefix = "concat_test",
      output_format = "vcf.gz",
      allow_overlaps = true
  }

  call validate_outputs { input:
      vcf_files = [mpileup_call.mpileup_vcf, concat.concatenated_vcf]
  }

  output {
    File variant_vcf = mpileup_call.mpileup_vcf
    File concatenated_vcf = concat.concatenated_vcf
    File concatenated_vcf_index = concat.concatenated_vcf_index
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that all expected bcftools output files were generated correctly"
    outputs: {
        report: "Validation report summarizing file checks and variant statistics"
    }
  }

  parameter_meta {
    vcf_files: "Array of VCF files to validate"
  }

  input {
    Array[File] vcf_files
  }

  command <<<
    set -eo pipefail

    echo "=== bcftools Variant Calling Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    # Arrays for bash processing
    vcf_files=~{sep=" " vcf_files}

    validation_passed=true
    total_variants=0

    # Check each sample
    for i in "${!vcf_files[@]}"; do
      vcf_file="${vcf_files[$i]}"

      echo "--- Sample: $vcf_file ---" >> validation_report.txt

      # Check VCF file exists and is not empty
      if [[ -f "$vcf_file" && -s "$vcf_file" ]]; then
        vcf_size=$(stat -c%s "$vcf_file")
        echo "VCF file: $vcf_file (${vcf_size} bytes)" >> validation_report.txt

        # Try to get variant counts if bcftools is available
        if command -v bcftools &> /dev/null; then
          variant_count=$(bcftools view -H "$vcf_file" | wc -l 2>/dev/null || echo "N/A")
          echo "  Total variants: $variant_count" >> validation_report.txt

          if [[ "$variant_count" =~ ^[0-9]+$ ]]; then
            total_variants=$((total_variants + variant_count))
          fi

          # Get basic statistics
          snp_count=$(bcftools view -H -v snps "$vcf_file" | \
                      wc -l 2>/dev/null || echo "N/A")
          indel_count=$(bcftools view -H -v indels "$vcf_file" | \
                      wc -l 2>/dev/null || echo "N/A")
          echo "  SNPs: $snp_count" >> validation_report.txt
          echo "  Indels: $indel_count" >> validation_report.txt
        fi
      else
        echo "VCF file: $vcf_file - MISSING OR EMPTY" >> validation_report.txt
        validation_passed=false
      fi

      echo "" >> validation_report.txt
    done

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
    docker: "getwilds/bcftools:1.19"
    memory: "2 GB"
    cpu: 1
  }
}
