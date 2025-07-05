## WILDS WDL module for bcftools variant calling and manipulation.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct SampleInfo {
    String name
    File bam
    File bai
}

struct RefGenome {
    String name
    File fasta
    File fasta_index
}

workflow bcftools_example {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow for variant calling using bcftools"
    url: "https://github.com/getwilds/wilds-wdl-library"
    outputs: {
        variant_vcfs: "VCF files containing variants called by bcftools mpileup for each sample",
        validation_report: "validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    samples: "List of sample objects, each containing name, BAM file, and BAM index"
    regions_bed: "Optional BED file specifying regions to analyze"
    reference_genome: "Optional reference genome object containing name, fasta, and fasta index files. If not provided, test data will be used."
    max_depth: "Maximum read depth for mpileup (default: 10000)"
    max_idepth: "Maximum per-sample depth for indel calling (default: 10000)"
    memory_gb: "Memory allocated for each task in the workflow in GB"
    cpu_cores: "Number of CPU cores allocated for each task in the workflow"
  }

  input {
    Array[SampleInfo] samples
    File? regions_bed
    RefGenome? reference_genome
    Int max_depth = 10000
    Int max_idepth = 10000
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  # If no reference genome provided, download test data
  if (!defined(reference_genome)) {
    call ww_testdata.download_ref_data { }
  }

  # Determine which genome files to use
  File genome_fasta = select_first([reference_genome.fasta, download_ref_data.fasta])
  File genome_fasta_index = select_first([reference_genome.fasta_index, download_ref_data.fasta_index])

  scatter (sample in samples) {
    call mpileup_call { input:
        bam_file = sample.bam,
        bam_index = sample.bai,
        reference_fasta = genome_fasta,
        reference_fasta_index = genome_fasta_index,
        regions_bed = regions_bed,
        max_depth = max_depth,
        max_idepth = max_idepth,
        memory_gb = memory_gb,
        cpu_cores = cpu_cores
    }
  }

  call validate_outputs { input:
      vcf_files = mpileup_call.mpileup_vcf
  }

  output {
    Array[File] variant_vcfs = mpileup_call.mpileup_vcf
    File validation_report = validate_outputs.report
  }
}

task mpileup_call {
  meta {
    description: "Call variants using bcftools mpileup and call"
    outputs: {
        mpileup_vcf: "Compressed VCF file containing variants called by mpileup",
        mpileup_vcf_index: "Index file for the mpileup VCF"
    }
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Reference genome FASTA index file"
    bam_file: "Input BAM file for the sample"
    bam_index: "Index file for the input BAM"
    regions_bed: "Optional BED file specifying regions to analyze"
    annotate_format: "FORMAT annotations to add (default: AD,DP)"
    ignore_rg: "Ignore read groups during analysis"
    disable_baq: "Disable BAQ (per-Base Alignment Quality) computation"
    max_depth: "Maximum read depth for mpileup"
    max_idepth: "Maximum per-sample depth for indel calling"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File bam_file
    File bam_index
    File? regions_bed
    String annotate_format = "FORMAT/AD,FORMAT/DP"
    Boolean ignore_rg = true
    Boolean disable_baq = true
    Int max_depth = 10000
    Int max_idepth = 10000
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail
    
    # Get the basename of the BAM file (without path and extension)
    sample_name=$(basename "~{bam_file}" .bam)
    
    # Build the bcftools mpileup command
    bcftools_cmd="bcftools mpileup \
      --max-depth ~{max_depth} \
      --max-idepth ~{max_idepth} \
      --annotate \"~{annotate_format}\" \
      --fasta-ref \"~{reference_fasta}\""
    
    # Add regions file if provided
    if [ -n "~{regions_bed}" ]; then
      bcftools_cmd="$bcftools_cmd --regions-file \"~{regions_bed}\""
    fi
    
    # Add optional flags
    if [ "~{ignore_rg}" == "true" ]; then
      bcftools_cmd="$bcftools_cmd --ignore-RG"
    fi
    
    if [ "~{disable_baq}" == "true" ]; then
      bcftools_cmd="$bcftools_cmd --no-BAQ"
    fi
    
    # Complete the command with input BAM and piping to bcftools call
    bcftools_cmd="$bcftools_cmd \"~{bam_file}\" | \
      bcftools call -Oz -mv -o \"${sample_name}.bcftools.vcf.gz\""
    
    echo "Running: $bcftools_cmd"
    eval "$bcftools_cmd"
    
    # Create index for the output VCF
    bcftools index "${sample_name}.bcftools.vcf.gz"
  >>>

  output {
    File mpileup_vcf = "~{basename(bam_file, '.bam')}.bcftools.vcf.gz"
    File mpileup_vcf_index = "~{basename(bam_file, '.bam')}.bcftools.vcf.gz.csi"
  }

  runtime {
    docker: "getwilds/bcftools:1.19"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
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
