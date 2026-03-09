## WILDS WDL module for DeepVariant variant calling.
## Google's deep-learning-based germline variant caller for next-generation
## DNA sequencing data. Supports WGS, WES, PacBio, and ONT sequencing technologies.

version 1.0

#### TASK DEFINITIONS ####

task run_deepvariant {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Runs DeepVariant to call germline variants from aligned reads"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deepvariant/ww-deepvariant.wdl"
    outputs: {
        output_vcf: "VCF file containing variant calls",
        output_vcf_index: "Index file for the output VCF"
    }
  }

  parameter_meta {
    sample_name: "Name identifier for the sample"
    input_bam: "Aligned reads in BAM format"
    input_bam_index: "Index file for the input BAM"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Index file for the reference FASTA"
    model_type: "Sequencing model type: WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA, or MASSEQ"
    regions: "Optional genomic regions to restrict variant calling (e.g., chr1)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    String sample_name
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    String model_type = "WGS"
    String? regions
    Int cpu_cores = 8
    Int memory_gb = 32
  }

  command <<<
    set -eo pipefail

    /opt/deepvariant/bin/run_deepvariant \
      --model_type=~{model_type} \
      --ref=~{ref_fasta} \
      --reads=~{input_bam} \
      --output_vcf="~{sample_name}.deepvariant.vcf.gz" \
      --num_shards=~{cpu_cores} \
      ~{if defined(regions) then "--regions=~{regions}" else ""}
  >>>

  output {
    File output_vcf = "~{sample_name}.deepvariant.vcf.gz"
    File output_vcf_index = "~{sample_name}.deepvariant.vcf.gz.tbi"
  }

  runtime {
    docker: "google/deepvariant:1.10.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task run_deepvariant_gvcf {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Runs DeepVariant to call germline variants and produce a gVCF for joint calling"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deepvariant/ww-deepvariant.wdl"
    outputs: {
        output_vcf: "VCF file containing variant calls",
        output_vcf_index: "Index file for the output VCF",
        output_gvcf: "gVCF file for downstream joint genotyping",
        output_gvcf_index: "Index file for the output gVCF"
    }
  }

  parameter_meta {
    sample_name: "Name identifier for the sample"
    input_bam: "Aligned reads in BAM format"
    input_bam_index: "Index file for the input BAM"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Index file for the reference FASTA"
    model_type: "Sequencing model type: WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA, or MASSEQ"
    regions: "Optional genomic regions to restrict variant calling (e.g., chr1)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    String sample_name
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    String model_type = "WGS"
    String? regions
    Int cpu_cores = 8
    Int memory_gb = 32
  }

  command <<<
    set -eo pipefail

    /opt/deepvariant/bin/run_deepvariant \
      --model_type=~{model_type} \
      --ref=~{ref_fasta} \
      --reads=~{input_bam} \
      --output_vcf="~{sample_name}.deepvariant.vcf.gz" \
      --output_gvcf="~{sample_name}.deepvariant.g.vcf.gz" \
      --num_shards=~{cpu_cores} \
      ~{if defined(regions) then "--regions=~{regions}" else ""}
  >>>

  output {
    File output_vcf = "~{sample_name}.deepvariant.vcf.gz"
    File output_vcf_index = "~{sample_name}.deepvariant.vcf.gz.tbi"
    File output_gvcf = "~{sample_name}.deepvariant.g.vcf.gz"
    File output_gvcf_index = "~{sample_name}.deepvariant.g.vcf.gz.tbi"
  }

  runtime {
    docker: "google/deepvariant:1.10.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
