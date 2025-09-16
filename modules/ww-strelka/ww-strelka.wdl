## WILDS WDL module for Strelka variant calling.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-strelka/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct StrelkaSample {
    String name
    File bam
    File bai
}

workflow strelka_example {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow for germline and somatic variant calling using Strelka"
    url: "https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-strelka"
    outputs: {
        germline_vcfs: "Germline variant calls in compressed VCF format with index files",
        germline_vcf_indices: "Index files for germline VCFs",
        somatic_snvs_vcfs: "Somatic SNV variant calls in compressed VCF format with index files",
        somatic_indels_vcfs: "Somatic indel variant calls in compressed VCF format with index files",
        somatic_snvs_vcf_indices: "Index files for somatic SNV VCFs",
        somatic_indels_vcf_indices: "Index files for somatic indel VCFs",
        validation_report: "Combined validation report for germline and somatic outputs"
    }
  }


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
    call strelka_germline { input:
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
    call strelka_somatic { input:
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

task strelka_germline {
  meta {
    description: "Run Strelka germline variant calling on a single sample"
    outputs: {
        variants_vcf: "Compressed VCF file containing germline variant calls",
        variants_vcf_index: "Index file for the variants VCF"
    }
  }

  parameter_meta {
    sample_name: "Name of the sample being analyzed"
    bam: "Input BAM file for variant calling"
    bai: "Index file for the input BAM"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Reference genome FASTA index file"
    target_regions_bed: "Optional BED file specifying regions to analyze"
    is_exome: "Whether this is exome sequencing data (enables exome-specific settings)"
    cpus: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    String sample_name
    File bam
    File bai
    File ref_fasta
    File ref_fasta_index
    File? target_regions_bed
    Boolean is_exome = false
    Int cpus = 4
    Int memory_gb = 8
  }

  String exome_flag = if is_exome then "--exome" else ""
  String regions_flag = if defined(target_regions_bed) then "--callRegions ~{target_regions_bed}" else ""

  command <<<
    set -euo pipefail

    # Create output directory
    mkdir -p strelka_germline_output

    # Configure Strelka germline workflow
    configureStrelkaGermlineWorkflow.py \
      --bam "~{bam}" \
      --referenceFasta "~{ref_fasta}" \
      --runDir strelka_germline_output \
      ~{exome_flag} \
      ~{regions_flag}

    # Execute the workflow
    strelka_germline_output/runWorkflow.py \
      -m local \
      -j ~{cpus}

    # Copy and rename output files
    cp strelka_germline_output/results/variants/variants.vcf.gz "~{sample_name}.strelka.germline.vcf.gz"
    cp strelka_germline_output/results/variants/variants.vcf.gz.tbi "~{sample_name}.strelka.germline.vcf.gz.tbi"

    # Basic validation
    echo "Validating output files..."
    if [[ ! -f ~{sample_name}.strelka.germline.vcf.gz ]]; then
      echo "ERROR: Variants VCF file not found"
      exit 1
    fi

    if [[ ! -f ~{sample_name}.strelka.germline.vcf.gz.tbi ]]; then
      echo "ERROR: Variants VCF index file not found"
      exit 1
    fi

    # Count variants
    variant_count=$(zcat ~{sample_name}.strelka.germline.vcf.gz | grep -v "^#" | wc -l || echo "0")
    echo "Total variants called: $variant_count"
  >>>

  output {
    File variants_vcf = "~{sample_name}.strelka.germline.vcf.gz"
    File variants_vcf_index = "~{sample_name}.strelka.germline.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/strelka:2.9.10"
    memory: "~{memory_gb} GB"
    cpu: cpus
  }
}

task strelka_somatic {
  meta {
    description: "Run Strelka somatic variant calling on tumor/normal pair"
    outputs: {
        somatic_snvs_vcf: "Compressed VCF file containing somatic SNV calls",
        somatic_indels_vcf: "Compressed VCF file containing somatic indel calls",
        somatic_snvs_vcf_index: "Index file for the somatic SNVs VCF",
        somatic_indels_vcf_index: "Index file for the somatic indels VCF"
    }
  }

  parameter_meta {
    tumor_sample_name: "Name of the tumor sample"
    tumor_bam: "Tumor sample BAM file"
    tumor_bai: "Tumor sample BAM index file"
    normal_sample_name: "Name of the normal sample"
    normal_bam: "Normal sample BAM file"
    normal_bai: "Normal sample BAM index file"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Reference genome FASTA index file"
    target_regions_bed: "Optional BED file specifying regions to analyze"
    is_exome: "Whether this is exome sequencing data (enables exome-specific settings)"
    cpus: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    String tumor_sample_name
    File tumor_bam
    File tumor_bai
    String normal_sample_name
    File normal_bam
    File normal_bai
    File ref_fasta
    File ref_fasta_index
    File? target_regions_bed
    Boolean is_exome = false
    Int cpus = 4
    Int memory_gb = 8
  }

  String exome_flag = if is_exome then "--exome" else ""
  String regions_flag = if defined(target_regions_bed) then "--callRegions ~{target_regions_bed}" else ""

  command <<<
    set -euo pipefail

    # Create output directory
    mkdir -p strelka_somatic_output

    # Configure Strelka somatic workflow
    configureStrelkaSomaticWorkflow.py \
      --tumorBam "~{tumor_bam}" \
      --normalBam "~{normal_bam}" \
      --referenceFasta "~{ref_fasta}" \
      --runDir strelka_somatic_output \
      ~{exome_flag} \
      ~{regions_flag}

    # Execute the workflow
    strelka_somatic_output/runWorkflow.py \
      -m local \
      -j ~{cpus}

    # Copy and rename output files
    cp strelka_somatic_output/results/variants/somatic.snvs.vcf.gz "~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.snvs.vcf.gz"
    cp strelka_somatic_output/results/variants/somatic.snvs.vcf.gz.tbi "~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.snvs.vcf.gz.tbi"
    cp strelka_somatic_output/results/variants/somatic.indels.vcf.gz "~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.indels.vcf.gz"
    cp strelka_somatic_output/results/variants/somatic.indels.vcf.gz.tbi "~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.indels.vcf.gz.tbi"

    # Basic validation
    echo "Validating output files..."
    for file in ~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.*.vcf.gz; do
      if [[ ! -f "$file" ]]; then
        echo "ERROR: Output file $file not found"
        exit 1
      fi
    done

    for file in ~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.*.vcf.gz.tbi; do
      if [[ ! -f "$file" ]]; then
        echo "ERROR: Index file $file not found"
        exit 1
      fi
    done

    # Count variants
    snv_count=$(zcat ~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.snvs.vcf.gz | grep -v "^#" | wc -l || echo "0")
    indel_count=$(zcat ~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.indels.vcf.gz | grep -v "^#" | wc -l || echo "0")
    echo "Somatic SNVs called: $snv_count"
    echo "Somatic indels called: $indel_count"
  >>>

  output {
    File somatic_snvs_vcf = "~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.snvs.vcf.gz"
    File somatic_indels_vcf = "~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.indels.vcf.gz"
    File somatic_snvs_vcf_index = "~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.snvs.vcf.gz.tbi"
    File somatic_indels_vcf_index = "~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.indels.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/strelka:2.9.10"
    memory: "~{memory_gb} GB"
    cpu: cpus
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
