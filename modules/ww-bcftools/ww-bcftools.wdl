## WILDS WDL module for bcftools variant calling and manipulation.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task mpileup_call {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Call variants using bcftools mpileup and call"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bcftools/ww-bcftools.wdl"
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
