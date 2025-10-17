## WILDS WDL module for bcftools variant calling and manipulation.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

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
