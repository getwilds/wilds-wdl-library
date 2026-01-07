## WILDS WDL module for structural variant calling using Smoove.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task smoove_call {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Call structural variants using Smoove for a single sample"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-smoove/ww-smoove.wdl"
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
    
    # Sort the resulting vcf (SV's are tricky)
    bcftools sort -Oz "results/~{sample_name}-smoove.vcf.gz" \
      > "results/~{sample_name}-smoove.sorted.vcf.gz"

    # Index the raw output
    tabix -p vcf "results/~{sample_name}-smoove.sorted.vcf.gz"

    # Filter to target regions if provided, otherwise use raw output
    if [ -n "~{target_regions_bed}" ]; then
      # Filter VCF to only include variants overlapping target regions
      bcftools view -R "~{target_regions_bed}" \
        "results/~{sample_name}-smoove.sorted.vcf.gz" | \
        bcftools sort -Oz -o "~{vcf_filename}"
      
      # Index the filtered VCF
      tabix -p vcf "~{vcf_filename}"
    else
      # Move and rename output files for consistency
      mv "results/~{sample_name}-smoove.sorted.vcf.gz" "~{vcf_filename}"
      mv "results/~{sample_name}-smoove.sorted.vcf.gz.tbi" "~{vcf_index_filename}"
    fi
  >>>

  output {
    File vcf = vcf_filename
    File vcf_index = vcf_index_filename
  }

  runtime {
    docker: "getwilds/smoove:0.2.8"
    cpu: cpu_cores
    memory: "${memory_gb}GB"
  }
}
