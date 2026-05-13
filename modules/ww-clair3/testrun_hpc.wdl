version 1.2

# HPC test workflow for ww-clair3. Mirrors testrun.wdl coverage but exercises
# the GPU code path on realistic input data — intended to validate end-to-end
# behavior on Fred Hutch HPC infrastructure (Sprocket or PROOF).

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/wdl-v1.2-support/modules/ww-clair3/ww-clair3.wdl" as ww_clair3
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/hpc-testruns/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####

workflow clair3_example {
  # Download full chr1 reference and BAM test data for a realistic call set
  call ww_testdata.download_ref_data as download_reference { input:
    chromo = "chr1"
  }

  call ww_testdata.download_bam_data as download_bam { input:
    filename = "sample1.bam"
  }

  # Pileup-only mode with GPU
  call ww_clair3.run_clair3 { input:
    sample_name = "sample1",
    input_bam = download_bam.bam,
    input_bam_index = download_bam.bai,
    ref_fasta = download_reference.fasta,
    ref_fasta_index = download_reference.fasta_index,
    platform = "ilmn",
    model_path = "/opt/models/ilmn",
    bed_file = download_reference.bed,
    pileup_only = true,
    gpu_enabled = true,
    cpu_cores = 4,
    memory_gb = 16
  }

  # Full-alignment mode with gVCF output and GPU
  call ww_clair3.run_clair3 as run_clair3_gvcf { input:
    sample_name = "sample1_gvcf",
    input_bam = download_bam.bam,
    input_bam_index = download_bam.bai,
    ref_fasta = download_reference.fasta,
    ref_fasta_index = download_reference.fasta_index,
    platform = "ilmn",
    model_path = "/opt/models/ilmn",
    bed_file = download_reference.bed,
    gvcf_enabled = true,
    gpu_enabled = true,
    cpu_cores = 4,
    memory_gb = 16
  }

  call validate_outputs { input:
      pileup_vcf = run_clair3.output_vcf,
      full_vcf = run_clair3_gvcf.output_vcf,
      gvcfs = run_clair3_gvcf.output_gvcf
  }

  output {
    File output_vcf = run_clair3.output_vcf
    File output_vcf_index = run_clair3.output_vcf_index
    File output_gvcf_vcf = run_clair3_gvcf.output_vcf
    File output_gvcf_vcf_index = run_clair3_gvcf.output_vcf_index
    Array[File] output_gvcf = run_clair3_gvcf.output_gvcf
    Array[File] output_gvcf_index = run_clair3_gvcf.output_gvcf_index
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that Clair3 produced non-empty VCF and gVCF outputs"
    outputs: {
        report: "Validation report summarizing record counts and file checks"
    }
  }

  parameter_meta {
    pileup_vcf: "VCF file from pileup-only run_clair3 invocation"
    full_vcf: "VCF file from full-alignment run_clair3 invocation"
    gvcfs: "Array of gVCF files from full-alignment run_clair3 invocation"
  }

  input {
    File pileup_vcf
    File full_vcf
    Array[File] gvcfs
  }

  command <<<
    set -eo pipefail

    echo "=== Clair3 HPC Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    pileup_records=$(zcat ~{pileup_vcf} | grep -vc "^#" || true)
    echo "Pileup VCF variant records: ${pileup_records}" >> validation_report.txt
    if [ "${pileup_records}" -eq 0 ]; then
      echo "WARNING: pileup VCF has no variant records" >> validation_report.txt
      validation_passed=false
    fi

    full_records=$(zcat ~{full_vcf} | grep -vc "^#" || true)
    echo "Full-alignment VCF variant records: ${full_records}" >> validation_report.txt
    if [ "${full_records}" -eq 0 ]; then
      echo "WARNING: full-alignment VCF has no variant records" >> validation_report.txt
      validation_passed=false
    fi

    gvcf_count=$(echo "~{sep=' ' gvcfs}" | wc -w)
    echo "gVCF files produced: ${gvcf_count}" >> validation_report.txt
    if [ "${gvcf_count}" -eq 0 ]; then
      echo "WARNING: no gVCF files produced" >> validation_report.txt
      validation_passed=false
    fi

    echo "" >> validation_report.txt
    echo "=== Validation Summary ===" >> validation_report.txt
    if [ "${validation_passed}" = "true" ]; then
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
    docker: "ubuntu:22.04"
    cpu: 1
    memory: "2 GB"
  }
}
