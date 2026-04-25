version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-samtools/ww-samtools.wdl" as ww_samtools
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-varscan/ww-varscan.wdl" as ww_varscan

workflow varscan_example {
  # Download test data
  call ww_testdata.download_ref_data { }

  # Download BAM files to simulate tumor and normal samples
  call ww_testdata.download_bam_data as download_sample_bam { }

  call ww_testdata.download_bam_data as download_normal_bam { input:
    filename = "NA12878_chr1_normal.bam"
  }

  # Make normal mpileup for VarScan somatic
  call ww_samtools.mpileup as normal_mpileup {
    input:
      bamfile = download_normal_bam.bam,
      ref_fasta = download_ref_data.fasta,
      sample_name = "test_normal",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Make normal mpileup for VarScan mpileup2cns
  call ww_samtools.mpileup as normal_mpileup_cns {
    input:
      bamfile = download_normal_bam.bam,
      ref_fasta = download_ref_data.fasta,
      sample_name = "test_normal_nobaq",
      disable_baq = true,
      cpu_cores = 2,
      memory_gb = 8
  }

  # Make 'tumor' mpileup
  call ww_samtools.mpileup as tumor_mpileup {
    input:
      bamfile = download_sample_bam.bam,
      ref_fasta = download_ref_data.fasta,
      sample_name = "test_tumor",
      cpu_cores = 2,
      memory_gb = 8
  }

  # VarScan somatic
  call ww_varscan.somatic {
    input:
      sample_name = "test_sample",
      normal_pileup = normal_mpileup.pileup,
      tumor_pileup = tumor_mpileup.pileup,
      cpu_cores = 2,
      memory_gb = 8
  }

  # VarScan mpileup2cns
  call ww_varscan.mpileup2cns {
    input:
      pileup = normal_mpileup_cns.pileup,
      sample_name = "test_normal_nobaq",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Validate outputs
  call validate_outputs {
    input:
      somatic_snvs_vcf = somatic.somatic_snvs_vcf,
      somatic_indels_vcf = somatic.somatic_indels_vcf,
      germline_vcf = mpileup2cns.vcf
  }

  output {
    File somatic_snvs = somatic.somatic_snvs_vcf
    File somatic_indels = somatic.somatic_indels_vcf
    File germline_vcf = mpileup2cns.vcf
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate all VarScan output files"
    outputs: {
        report: "Validation report summarizing file checks and variant statistics"
    }
  }

  parameter_meta {
    somatic_snvs_vcf: "Somatic SNVs VCF file to validate"
    somatic_indels_vcf: "Somatic indels VCF file to validate"
    germline_vcf: "Germline VCF file to validate"
  }

  input {
    File somatic_snvs_vcf
    File somatic_indels_vcf
    File germline_vcf
  }

  command <<<
    set -eo pipefail

    echo "=== VarScan Somatic Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Check somatic SNVs VCF file
    echo "--- Somatic SNVs VCF File ---" >> validation_report.txt
    if [[ -f "~{somatic_snvs_vcf}" && -s "~{somatic_snvs_vcf}" ]]; then
      snvs_size=$(wc -c < "~{somatic_snvs_vcf}")
      echo "Somatic SNVs VCF: ~{somatic_snvs_vcf} (${snvs_size} bytes)" >> validation_report.txt

      # Count variants (non-header lines)
      snv_count=$(grep -cv '^#' "~{somatic_snvs_vcf}" || echo "0")
      echo "SNV variants called: ${snv_count}" >> validation_report.txt
    else
      echo "Somatic SNVs VCF: ~{somatic_snvs_vcf} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Check somatic indels VCF file
    echo "--- Somatic Indels VCF File ---" >> validation_report.txt
    if [[ -f "~{somatic_indels_vcf}" && -s "~{somatic_indels_vcf}" ]]; then
      indels_size=$(wc -c < "~{somatic_indels_vcf}")
      echo "Somatic Indels VCF: ~{somatic_indels_vcf} (${indels_size} bytes)" >> validation_report.txt

      # Count variants (non-header lines)
      indel_count=$(grep -cv '^#' "~{somatic_indels_vcf}" || echo "0")
      echo "Indel variants called: ${indel_count}" >> validation_report.txt
    else
      echo "Somatic Indels VCF: ~{somatic_indels_vcf} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Check germline VCF file
    echo "--- Germline VCF File ---" >> validation_report.txt
    if [[ -f "~{germline_vcf}" && -s "~{germline_vcf}" ]]; then
      germ_size=$(wc -c < "~{germline_vcf}")
      echo "Germline variants VCF: ~{germline_vcf} (${germ_size} bytes)" >> validation_report.txt

      # Count variants (non-header lines)
      germ_count=$(grep -cv '^#' "~{germline_vcf}" || echo "0")
      echo "Germline variants called: ${germ_count}" >> validation_report.txt
    else
      echo "Germline variants VCF: ~{germline_vcf} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi
    echo "" >> validation_report.txt

    # Overall summary
    echo "=== Validation Summary ===" >> validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> validation_report.txt
      echo "All VarScan outputs were generated successfully." >> validation_report.txt
    else
      echo "Overall Status: FAILED" >> validation_report.txt
      echo "One or more output files are missing or empty." >> validation_report.txt
      exit 1
    fi

    # Also output to stdout for immediate feedback
    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/varscan:2.4.6"
    memory: "2 GB"
    cpu: 1
  }
}
