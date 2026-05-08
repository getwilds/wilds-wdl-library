version 1.2

# HPC test workflow for ww-deepvariant. Mirrors testrun.wdl coverage with
# GPU acceleration enabled — intended to validate end-to-end behavior on
# Fred Hutch HPC infrastructure (Sprocket or PROOF).

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/split-cicd-hpc-testruns/modules/ww-deepvariant/ww-deepvariant.wdl" as ww_deepvariant
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/split-cicd-hpc-testruns/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct DeepVariantSample {
    String name
    File bam
    File bai
}

workflow deepvariant_example {
  # Download chr1 reference and BAM test data
  call ww_testdata.download_ref_data as download_reference { input:
    chromo = "chr1"
  }

  call ww_testdata.download_bam_data as download_bam_1 { input:
    filename = "sample1.bam"
  }

  call ww_testdata.download_bam_data as download_bam_2 { input:
    filename = "sample2.bam"
  }

  Array[DeepVariantSample] samples = [
    {
      "name": "sample1",
      "bam": download_bam_1.bam,
      "bai": download_bam_1.bai
    },
    {
      "name": "sample2",
      "bam": download_bam_2.bam,
      "bai": download_bam_2.bai
    }
  ]

  # Scatter VCF-only DeepVariant calls across samples on GPU
  scatter (sample in samples) {
    call ww_deepvariant.run_deepvariant { input:
      sample_name = sample.name,
      input_bam = sample.bam,
      input_bam_index = sample.bai,
      ref_fasta = download_reference.fasta,
      ref_fasta_index = download_reference.fasta_index,
      model_type = "WGS",
      regions = "chr1",
      gpu_enabled = true,
      cpu_cores = 8,
      memory_gb = 32
    }
  }

  # gVCF call on the first sample, also on GPU
  call ww_deepvariant.run_deepvariant as run_deepvariant_gvcf { input:
    sample_name = "sample1_gvcf",
    input_bam = download_bam_1.bam,
    input_bam_index = download_bam_1.bai,
    ref_fasta = download_reference.fasta,
    ref_fasta_index = download_reference.fasta_index,
    model_type = "WGS",
    output_gvcf_enabled = true,
    regions = "chr1",
    gpu_enabled = true,
    cpu_cores = 8,
    memory_gb = 32
  }

  call validate_outputs { input:
      vcfs = run_deepvariant.output_vcf,
      gvcf_vcf = run_deepvariant_gvcf.output_vcf,
      gvcfs = run_deepvariant_gvcf.output_gvcf
  }

  output {
    Array[File] output_vcfs = run_deepvariant.output_vcf
    Array[File] output_vcf_indices = run_deepvariant.output_vcf_index
    File output_gvcf_vcf = run_deepvariant_gvcf.output_vcf
    Array[File] output_gvcfs = run_deepvariant_gvcf.output_gvcf
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate that DeepVariant produced non-empty VCF and gVCF outputs"
    outputs: {
        report: "Validation report summarizing record counts and file checks"
    }
  }

  parameter_meta {
    vcfs: "Array of per-sample VCF files from scattered run_deepvariant calls"
    gvcf_vcf: "VCF file from the gVCF-enabled run_deepvariant invocation"
    gvcfs: "Array of gVCF files from the gVCF-enabled run_deepvariant invocation"
  }

  input {
    Array[File] vcfs
    File gvcf_vcf
    Array[File] gvcfs
  }

  command <<<
    set -eo pipefail

    echo "=== DeepVariant HPC Validation Report ===" > validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    for vcf in ~{sep=' ' vcfs}; do
      records=$(zcat "${vcf}" | grep -vc "^#" || true)
      echo "Sample VCF ${vcf}: ${records} variant records" >> validation_report.txt
      if [ "${records}" -eq 0 ]; then
        echo "WARNING: VCF ${vcf} has no variant records" >> validation_report.txt
        validation_passed=false
      fi
    done

    gvcf_vcf_records=$(zcat ~{gvcf_vcf} | grep -vc "^#" || true)
    echo "gVCF-mode VCF variant records: ${gvcf_vcf_records}" >> validation_report.txt
    if [ "${gvcf_vcf_records}" -eq 0 ]; then
      echo "WARNING: gVCF-mode VCF has no variant records" >> validation_report.txt
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
