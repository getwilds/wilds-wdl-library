## WILDS WDL module for Strelka variant calling.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task strelka_germline {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Run Strelka germline variant calling on a single sample"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-strelka/ww-strelka.wdl"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Run Strelka somatic variant calling on tumor/normal pair"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-strelka/ww-strelka.wdl"
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
