## WILDS WDL module for structural variant calling using Delly.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task delly_call {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Calls structural variants using Delly for a single sample"
    outputs: {
        vcf: "VCF file containing structural variant calls",
        vcf_index: "Index file for the VCF output",
        summary: "Summary text file with Delly run details and statistics"
    }
  }

  parameter_meta {
    aligned_bam: "Input aligned BAM file containing reads for SV calling"
    aligned_bam_index: "Index file for the aligned BAM above"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    target_regions_bed: "Optional BED file of regions to target for SV calling"
    exclude_regions_bed: "Optional BED file to exclude problematic regions"
    sv_type: "Structural variant type to call (DEL, DUP, INV, TRA, INS) or empty for all"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File aligned_bam
    File aligned_bam_index
    File reference_fasta
    File reference_fasta_index
    File? target_regions_bed
    File? exclude_regions_bed
    String sv_type = ""
    Int cpu_cores = 8
    Int memory_gb = 16
  }

  String sample_name = basename(aligned_bam, ".bam")
  String output_prefix = "~{sample_name}.delly"
  String sv_type_arg = if sv_type != "" then "-t " + sv_type else ""
  String exclude_arg = if defined(exclude_regions_bed) then "-x " + select_first([exclude_regions_bed]) else ""

  command <<<
    set -euo pipefail

    # Set OpenMP threads for Delly parallelization
    export OMP_NUM_THREADS=~{cpu_cores}

    # Run Delly structural variant calling
    delly call \
      ~{sv_type_arg} \
      ~{exclude_arg} \
      -g "~{reference_fasta}" \
      -o "~{output_prefix}.bcf" \
      "~{aligned_bam}"

    # Filter and convert BCF to VCF
    bcftools view \
        ~{if defined(target_regions_bed) then "-R " + target_regions_bed else ""} \
        -Oz \
        -o "~{output_prefix}.vcf.gz" \
        "~{output_prefix}.bcf"
    bcftools index "~{output_prefix}.vcf.gz"

    # Generate summary statistics
    echo "Delly SV calling completed for sample: ~{sample_name}" > "~{output_prefix}.summary.txt"
    echo "SV type filter: ~{if sv_type != "" then sv_type else "ALL"}" >> "~{output_prefix}.summary.txt"
    echo "Reference genome: ~{basename(reference_fasta)}" >> "~{output_prefix}.summary.txt"
    echo "Exclude regions: ~{if defined(exclude_regions_bed) then "YES" else "NO"}" >> "~{output_prefix}.summary.txt"
    VARIANT_COUNT=$(bcftools view -H "~{output_prefix}.vcf.gz" | wc -l)
    echo "Total variants called: $VARIANT_COUNT" >> "~{output_prefix}.summary.txt"
  >>>

  output {
    File vcf = "~{output_prefix}.vcf.gz"
    File vcf_index = "~{output_prefix}.vcf.gz.csi"
    File summary = "~{output_prefix}.summary.txt"
  }

  runtime {
    docker: "getwilds/delly:1.2.9"
    cpu: cpu_cores
    memory: "~{memory_gb}GB"
  }
}
