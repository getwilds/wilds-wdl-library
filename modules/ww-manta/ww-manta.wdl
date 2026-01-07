## WILDS WDL module for structural variant calling using Manta.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task manta_call {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Call structural variants using Manta on a single sample"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-manta/ww-manta.wdl"
    outputs: {
        vcf: "Structural variant calls in compressed VCF format",
        vcf_index: "Index file for the VCF output"
    }
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    aligned_bam: "Input aligned BAM file containing reads for variant calling"
    aligned_bam_index: "Index file for the aligned BAM"
    sample_name: "Name of the sample provided for output files"
    call_regions_bed: "Optional BED file to restrict variant calling to specific regions"
    call_regions_index: "Index file for the optional BED file"
    is_rna: "Boolean flag for RNA-seq mode (enables RNA-specific settings)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File aligned_bam
    File aligned_bam_index
    String sample_name
    File? call_regions_bed
    File? call_regions_index
    Boolean is_rna = false
    Int cpu_cores = 8
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail
    
    # Create working directory
    mkdir -p manta_work
    
    # Configure Manta workflow
    configManta.py \
      --bam "~{aligned_bam}" \
      --referenceFasta "~{reference_fasta}" \
      ~{if defined(call_regions_bed) then "--callRegions " + call_regions_bed else ""} \
      ~{true="--rna" false="" is_rna} \
      --runDir manta_work
    
    # Execute Manta workflow
    cd manta_work
    ./runWorkflow.py -m local -j ~{cpu_cores}
    
    # Copy outputs to working directory
    cp results/variants/diploidSV.vcf.gz "../~{sample_name}.manta.vcf.gz"
    cp results/variants/diploidSV.vcf.gz.tbi "../~{sample_name}.manta.vcf.gz.tbi"
  >>>

  output {
    File vcf = "~{sample_name}.manta.vcf.gz"
    File vcf_index = "~{sample_name}.manta.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/manta:1.6.0"
    memory: "~{memory_gb}GB"
    cpu: cpu_cores
  }
}
