## WILDS WDL module for Clair3 variant calling.
## A deep learning-based germline small variant caller for long-read
## sequencing data. Uses a two-stage pileup and full-alignment approach
## for high-accuracy SNP and indel calling from ONT, PacBio, and Illumina reads.

version 1.0

#### TASK DEFINITIONS ####

task run_clair3 {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Runs Clair3 to call germline variants from aligned reads using deep learning pileup and full-alignment models"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-clair3/ww-clair3.wdl"
    outputs: {
        output_vcf: "VCF file containing merged variant calls from pileup and full-alignment models",
        output_vcf_index: "Index file for the output VCF",
        output_gvcf: "Array containing the gVCF file for downstream joint genotyping (empty when gvcf_enabled is false)",
        output_gvcf_index: "Array containing the gVCF index file (empty when gvcf_enabled is false)"
    }
  }

  parameter_meta {
    sample_name: "Name identifier for the sample"
    input_bam: "Aligned reads in BAM format (must be sorted and indexed)"
    input_bam_index: "Index file for the input BAM"
    ref_fasta: "Reference genome FASTA file (must be indexed with .fai)"
    ref_fasta_index: "Index file for the reference FASTA"
    platform: "Sequencing platform: ont, hifi, or ilmn"
    model_path: "Path to Clair3 model directory (default uses bundled ONT model)"
    gvcf_enabled: "Whether to also produce a gVCF file for downstream joint genotyping"
    bed_file: "Optional BED file to restrict variant calling to specific regions"
    ctg_name: "Optional contig/chromosome name to restrict calling"
    include_all_ctgs: "Call variants on all contigs (required for non-human species, default calls chr1-22,X,Y only)"
    pileup_only: "Use only the pileup model for faster variant calling (skips full-alignment stage)"
    gpu_enabled: "Enable GPU acceleration for Clair3 inference (requests GPU in runtime)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    String sample_name
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    String platform = "ont"
    String model_path = "/opt/models/ont"
    Boolean gvcf_enabled = false
    File? bed_file
    String? ctg_name
    Boolean include_all_ctgs = false
    Boolean pileup_only = false
    Boolean gpu_enabled = false
    Int cpu_cores = 8
    Int memory_gb = 32
  }

  command <<<
    set -eo pipefail

    /opt/bin/run_clair3.sh \
      --bam_fn=~{input_bam} \
      --ref_fn=~{ref_fasta} \
      --output=clair3_output \
      --threads=~{cpu_cores} \
      --platform=~{platform} \
      --model_path=~{model_path} \
      --sample_name=~{sample_name} \
      ~{if gvcf_enabled then "--gvcf" else ""} \
      ~{if defined(bed_file) then "--bed_fn=~{bed_file}" else ""} \
      ~{if defined(ctg_name) then "--ctg_name=~{ctg_name}" else ""} \
      ~{if include_all_ctgs then "--include_all_ctgs" else ""} \
      ~{if pileup_only then "--pileup_only" else ""} \
      ~{if gpu_enabled then "--use_gpu" else ""}

    # Move final outputs to working directory with sample-specific names
    # pileup_only mode outputs pileup.vcf.gz; full mode outputs merge_output.vcf.gz
    if [ "~{pileup_only}" = "true" ]; then
      cp clair3_output/pileup.vcf.gz "~{sample_name}.clair3.vcf.gz"
      cp clair3_output/pileup.vcf.gz.tbi "~{sample_name}.clair3.vcf.gz.tbi"
    else
      cp clair3_output/merge_output.vcf.gz "~{sample_name}.clair3.vcf.gz"
      cp clair3_output/merge_output.vcf.gz.tbi "~{sample_name}.clair3.vcf.gz.tbi"
    fi

    if [ "~{gvcf_enabled}" = "true" ]; then
      cp clair3_output/merge_output.gvcf.gz "~{sample_name}.clair3.g.vcf.gz"
      cp clair3_output/merge_output.gvcf.gz.tbi "~{sample_name}.clair3.g.vcf.gz.tbi"
    fi
  >>>

  output {
    File output_vcf = "~{sample_name}.clair3.vcf.gz"
    File output_vcf_index = "~{sample_name}.clair3.vcf.gz.tbi"
    Array[File] output_gvcf = glob("*.g.vcf.gz")
    Array[File] output_gvcf_index = glob("*.g.vcf.gz.tbi")
  }

  runtime {
    docker: "getwilds/clair3:2.0.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
    gpus: if gpu_enabled then "1" else "0"
  }
}
