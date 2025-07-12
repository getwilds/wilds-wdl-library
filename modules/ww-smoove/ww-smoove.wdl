## WILDS WDL module for structural variant calling using Smoove.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct SmooveSample {
    String name
    File bam
    File bai
}

workflow smoove_example {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow for structural variant calling via Smoove"
    url: "https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-smoove"
    outputs: {
        smoove_vcfs: "Structural variant calls in VCF format from Smoove",
        smoove_vcf_indexes: "Index file for the VCF output",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    samples: "List of sample objects, each containing name, BAM file, and BAM index"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Reference genome FASTA index file"
    exclude_bed: "Optional BED file defining regions to exclude from calling"
    include_bed: "Optional BED file defining regions to include for calling"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    memory_gb: "Memory allocated for each task in the workflow in GB"
  }

  input {
    Array[SmooveSample]? samples
    File? ref_fasta
    File? ref_fasta_index
    File? exclude_bed
    File? include_bed
    Int cpus = 2
    Int memory_gb = 8
  }

  # Determine which genome files to use
  if (!defined(ref_fasta) || !defined(ref_fasta_index)) {
    call ww_testdata.download_ref_data { }
  }
  File genome_fasta = select_first([ref_fasta, download_ref_data.fasta])
  File genome_fasta_index = select_first([ref_fasta_index, download_ref_data.fasta_index])

  # If no samples provided, download demonstration data
  if (!defined(samples)) {
    call ww_testdata.download_bam_data { }
  }

  # Create samples array - either from input or from BWA alignment
  Array[SmooveSample] final_samples = if defined(samples) then select_first([samples]) else [
    {
      "name": "demo_sample",
      "bam": select_first([download_bam_data.bam]),
      "bai": select_first([download_bam_data.bai])
    }
  ]

  scatter (sample in final_samples) {
    call smoove_call { input:
        aligned_bam = sample.bam,
        aligned_bam_index = sample.bai,
        sample_name = sample.name,
        reference_fasta = genome_fasta,
        reference_fasta_index = genome_fasta_index,
        target_regions_bed = include_bed,
        exclude_bed = exclude_bed,
        cpu_cores = cpus,
        memory_gb = memory_gb
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

task smoove_call {
  meta {
    description: "Call structural variants using Smoove for a single sample"
    outputs: {
        vcf: "Structural variant calls in compressed VCF format",
        vcf_index: "Index file for the VCF output"
    }
  }

  parameter_meta {
    aligned_bam: "Input aligned BAM file containing reads for variant calling"
    aligned_bam_index: "Index file for the aligned BAM above"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    sample_name: "Name of the sample provided for output files"
    target_regions_bed: "Optional BED file defining target regions to include in final output"
    exclude_bed: "Optional BED file defining regions to exclude from calling"
    exclude_chroms: "Optional comma-separated list of chromosomes to exclude"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File aligned_bam
    File aligned_bam_index
    File reference_fasta
    File reference_fasta_index
    String sample_name
    File? target_regions_bed
    File? exclude_bed
    String? exclude_chroms
    Int cpu_cores = 8
    Int memory_gb = 16
  }

  String vcf_filename = "${sample_name}.smoove.vcf.gz"
  String vcf_index_filename = "${vcf_filename}.tbi"

  command <<<
    set -euo pipefail

    # Create output directory
    mkdir -p results

    # Build smoove command
    smoove call \
      --outdir results \
      --name "~{sample_name}" \
      --fasta "~{reference_fasta}" \
      --processes ~{cpu_cores} \
      ~{if defined(exclude_bed) then "--exclude " + exclude_bed else ""} \
      ~{if defined(exclude_chroms) then "--excludechroms " + exclude_chroms else ""} \
      "~{aligned_bam}"
    
    # Index the raw output
    tabix -p vcf "results/~{sample_name}-smoove.vcf.gz"

    # Filter to target regions if provided, otherwise use raw output
    if [ -n "~{target_regions_bed}" ]; then
      # Filter VCF to only include variants overlapping target regions
      bcftools view \
        -R "~{target_regions_bed}" \
        -Oz \
        -o "~{vcf_filename}" \
        "results/~{sample_name}-smoove.vcf.gz"
      
      # Index the filtered VCF
      tabix -p vcf "~{vcf_filename}"
    else
      # Move and rename output files for consistency
      mv "results/~{sample_name}-smoove.vcf.gz" "~{vcf_filename}"
      mv "results/~{sample_name}-smoove.vcf.gz.tbi" "~{vcf_index_filename}"
    fi
  >>>

  output {
    File vcf = vcf_filename
    File vcf_index = vcf_index_filename
  }

  runtime {
    docker: "brentp/smoove:latest"
    cpu: cpu_cores
    memory: "${memory_gb}GB"
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
