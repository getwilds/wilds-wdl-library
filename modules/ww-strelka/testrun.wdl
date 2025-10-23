version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-strelka/ww-strelka.wdl" as ww_strelka

struct StrelkaSample {
    String name
    File bam
    File bai
}

workflow strelka_example {
  # Download test reference data
  call ww_testdata.download_ref_data { }

  # Download test sample data
  call ww_testdata.download_bam_data as sample_data { }

  # Create test samples array
  Array[StrelkaSample] final_samples = [
    {
      "name": "demo_sample",
      "bam": sample_data.bam,
      "bai": sample_data.bai
    }
  ]

  # Germline variant calling
  scatter (sample in final_samples) {
    call ww_strelka.strelka_germline { input:
        sample_name = sample.name,
        bam = sample.bam,
        bai = sample.bai,
        ref_fasta = download_ref_data.fasta,
        ref_fasta_index = download_ref_data.fasta_index,
        is_exome = false,
        cpus = 4,
        memory_gb = 8
    }
  }

  # Validate germline outputs
  call validate_germline_outputs { input:
      germline_vcfs = strelka_germline.variants_vcf,
      germline_indices = strelka_germline.variants_vcf_index
  }

  # Somatic variant calling (tumor/normal pairs)
  # Download normal test data
  call ww_testdata.download_bam_data as normal_data { input:
    filename = "NA12878_chr1_normal.bam"
  }

  # Create normal samples array
  Array[StrelkaSample] final_normal_samples = [
    {
      "name": "normal_sample",
      "bam": normal_data.bam,
      "bai": normal_data.bai
    }
  ]

  # Somatic variant calling for tumor/normal pairs
  scatter (i in range(length(final_samples))) {
    call ww_strelka.strelka_somatic { input:
        tumor_sample_name = final_samples[i].name,
        tumor_bam = final_samples[i].bam,
        tumor_bai = final_samples[i].bai,
        normal_sample_name = final_normal_samples[i].name,
        normal_bam = final_normal_samples[i].bam,
        normal_bai = final_normal_samples[i].bai,
        ref_fasta = download_ref_data.fasta,
        ref_fasta_index = download_ref_data.fasta_index,
        is_exome = false,
        cpus = 4,
        memory_gb = 8
    }
  }

  # Validate somatic outputs
  call validate_somatic_outputs { input:
      somatic_snvs_vcfs = strelka_somatic.somatic_snvs_vcf,
      somatic_indels_vcfs = strelka_somatic.somatic_indels_vcf,
      somatic_snvs_indices = strelka_somatic.somatic_snvs_vcf_index,
      somatic_indels_indices = strelka_somatic.somatic_indels_vcf_index
  }

  call combine_validation_reports { input:
      germline_report = validate_germline_outputs.report,
      somatic_report = validate_somatic_outputs.report
  }

  output {
    Array[File]? germline_vcfs = strelka_germline.variants_vcf
    Array[File]? germline_vcf_indices = strelka_germline.variants_vcf_index
    Array[File]? somatic_snvs_vcfs = strelka_somatic.somatic_snvs_vcf
    Array[File]? somatic_indels_vcfs = strelka_somatic.somatic_indels_vcf
    Array[File]? somatic_snvs_vcf_indices = strelka_somatic.somatic_snvs_vcf_index
    Array[File]? somatic_indels_vcf_indices = strelka_somatic.somatic_indels_vcf_index
    File validation_report = combine_validation_reports.report  # Single combined report
  }
}

task validate_germline_outputs {
  meta {
    description: "Validate Strelka germline outputs and generate summary report"
    outputs: {
        report: "Validation summary with file checks and basic statistics for germline outputs"
    }
  }

  parameter_meta {
    germline_vcfs: "Array of germline VCF files to validate"
    germline_indices: "Array of germline VCF index files to validate"
  }

  input {
    Array[File] germline_vcfs
    Array[File] germline_indices
  }

  command <<<
    set -euo pipefail

    echo "=== Strelka Germline Validation Report ===" > validation_report.txt
    echo "Generated: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Check germline VCF files
    echo "Germline VCF Files:" >> validation_report.txt
    for vcf in ~{sep=" " germline_vcfs}; do
        if [[ -f "$vcf" && -s "$vcf" ]]; then
            # Count variants
            variant_count=$(zcat "$vcf" | grep -v "^#" | wc -l || echo "0")
            echo "  $vcf - PASSED ($variant_count variants)" >> validation_report.txt
        else
            echo "  $vcf - MISSING OR EMPTY" >> validation_report.txt
            validation_passed=false
        fi
    done

    # Check germline VCF index files  
    echo "Germline VCF Index Files:" >> validation_report.txt
    for idx in ~{sep=" " germline_indices}; do
        if [[ -f "$idx" && -s "$idx" ]]; then
            echo "  $idx - PASSED" >> validation_report.txt
        else
            echo "  $idx - MISSING OR EMPTY" >> validation_report.txt
            validation_passed=false
        fi
    done

    echo "" >> validation_report.txt
    echo "=== Summary ===" >> validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
        echo "Validation Status: PASSED" >> validation_report.txt
        echo "All germline output files are present and non-empty." >> validation_report.txt
    else
        echo "Validation Status: FAILED" >> validation_report.txt
        echo "One or more germline output files are missing or empty." >> validation_report.txt
        exit 1
    fi

    echo "Germline validation completed successfully!" >> validation_report.txt
    >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/strelka:2.9.10"
    memory: "2 GB"
    cpu: 1
  }
}

task validate_somatic_outputs {
  meta {
    description: "Validate Strelka somatic outputs and generate summary report"
    outputs: {
        report: "Validation summary with file checks and basic statistics for somatic outputs"
    }
  }

  parameter_meta {
    somatic_snvs_vcfs: "Array of somatic SNV VCF files to validate"
    somatic_indels_vcfs: "Array of somatic indel VCF files to validate"
    somatic_snvs_indices: "Array of somatic SNV VCF index files to validate"
    somatic_indels_indices: "Array of somatic indel VCF index files to validate"
  }

  input {
    Array[File] somatic_snvs_vcfs
    Array[File] somatic_indels_vcfs
    Array[File] somatic_snvs_indices
    Array[File] somatic_indels_indices
  }

  command <<<
    set -euo pipefail

    echo "=== Strelka Somatic Validation Report ===" > validation_report.txt
    echo "Generated: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Check somatic SNV VCF files
    echo "Somatic SNV VCF Files:" >> validation_report.txt
    for vcf in ~{sep=" " somatic_snvs_vcfs}; do
        if [[ -f "$vcf" && -s "$vcf" ]]; then
            # Count variants
            variant_count=$(zcat "$vcf" | grep -v "^#" | wc -l || echo "0")
            echo "  $vcf - PASSED ($variant_count SNVs)" >> validation_report.txt
        else
            echo "  $vcf - MISSING OR EMPTY" >> validation_report.txt
            validation_passed=false
        fi
    done

    # Check somatic indel VCF files
    echo "Somatic Indel VCF Files:" >> validation_report.txt
    for vcf in ~{sep=" " somatic_indels_vcfs}; do
        if [[ -f "$vcf" && -s "$vcf" ]]; then
            # Count variants
            variant_count=$(zcat "$vcf" | grep -v "^#" | wc -l || echo "0")
            echo "  $vcf - PASSED ($variant_count indels)" >> validation_report.txt
        else
            echo "  $vcf - MISSING OR EMPTY" >> validation_report.txt
            validation_passed=false
        fi
    done

    # Check somatic SNV index files
    echo "Somatic SNV Index Files:" >> validation_report.txt
    for idx in ~{sep=" " somatic_snvs_indices}; do
        if [[ -f "$idx" && -s "$idx" ]]; then
            echo "  $idx - PASSED" >> validation_report.txt
        else
            echo "  $idx - MISSING OR EMPTY" >> validation_report.txt
            validation_passed=false
        fi
    done

    # Check somatic indel index files
    echo "Somatic Indel Index Files:" >> validation_report.txt
    for idx in ~{sep=" " somatic_indels_indices}; do
        if [[ -f "$idx" && -s "$idx" ]]; then
            echo "  $idx - PASSED" >> validation_report.txt
        else
            echo "  $idx - MISSING OR EMPTY" >> validation_report.txt
            validation_passed=false
        fi
    done

    echo "" >> validation_report.txt
    echo "=== Summary ===" >> validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
        echo "Validation Status: PASSED" >> validation_report.txt
        echo "All somatic output files are present and non-empty." >> validation_report.txt
    else
        echo "Validation Status: FAILED" >> validation_report.txt
        echo "One or more somatic output files are missing or empty." >> validation_report.txt
        exit 1
    fi

    echo "Somatic validation completed successfully!" >> validation_report.txt
    >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/strelka:2.9.10"
    memory: "2 GB"
    cpu: 1
  }
}

task combine_validation_reports {
  meta {
    description: "Combine germline and somatic validation reports into a single report"
    outputs: {
        report: "Combined validation summary with both germline and somatic results"
    }
  }

  parameter_meta {
    germline_report: "Optional germline validation report"
    somatic_report: "Optional somatic validation report"
  }

  input {
    File? germline_report
    File? somatic_report
  }

  command <<<
    set -euo pipefail

    echo "=== Strelka Combined Validation Report ===" > validation_report.txt
    echo "Generated: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    # Include germline report if provided
    if [[ "~{defined(germline_report)}" == "true" ]]; then
        echo "=== GERMLINE RESULTS ===" >> validation_report.txt
        # Skip the header from the germline report and include the rest
        tail -n +4 "~{germline_report}" >> validation_report.txt
        echo "" >> validation_report.txt
    fi

    # Include somatic report if provided  
    if [[ "~{defined(somatic_report)}" == "true" ]]; then
        echo "=== SOMATIC RESULTS ===" >> validation_report.txt
        # Skip the header from the somatic report and include the rest
        tail -n +4 "~{somatic_report}" >> validation_report.txt
        echo "" >> validation_report.txt
    fi

    # If neither report is provided, add a note
    if [[ "~{defined(germline_report)}" == "false" && "~{defined(somatic_report)}" == "false" ]]; then
        echo "No validation reports provided." >> validation_report.txt
        echo "This may indicate that neither germline nor somatic calling was performed." >> validation_report.txt
    fi

    echo "=== COMBINED SUMMARY ===" >> validation_report.txt
    
    # Check if any validation failed by looking for "FAILED" in the reports
    overall_status="PASSED"
    
    if [[ "~{defined(germline_report)}" == "true" ]]; then
        if grep -q "Validation Status: FAILED" "~{germline_report}"; then
            overall_status="FAILED"
        fi
    fi
    
    if [[ "~{defined(somatic_report)}" == "true" ]]; then
        if grep -q "Validation Status: FAILED" "~{somatic_report}"; then
            overall_status="FAILED"
        fi
    fi
    
    echo "Overall Validation Status: $overall_status" >> validation_report.txt
    
    if [[ "$overall_status" == "FAILED" ]]; then
        echo "One or more validation checks failed. Please review the detailed results above." >> validation_report.txt
        exit 1
    else
        echo "All validation checks passed successfully." >> validation_report.txt
    fi
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/strelka:2.9.10"
    memory: "2 GB"
    cpu: 1
  }
}
