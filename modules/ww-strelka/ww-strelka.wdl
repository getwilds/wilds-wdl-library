## WILDS WDL module for Strelka variant calling.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

# import "../ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct StrelkaSample {
    String name
    File bam
    File bai
}

workflow strelka_example {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "WDL workflow for germline and somatic variant calling using Strelka"
    url: "https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-strelka"
    outputs: {
        germline_vcfs: "Germline variant calls in VCF format with index files",
        somatic_vcfs: "Somatic variant calls in VCF format with index files (tumor/normal pairs)",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    samples: "Optional array of StrelkaSample objects, each containing sample name, BAM file, and BAM index. If not provided, workflow will download test data."
    normal_samples: "Optional array of normal samples for somatic variant calling. Should match with tumor samples by array index."
    ref_fasta: "Optional reference genome FASTA file. If not provided, test data will be used."
    ref_fasta_index: "Optional reference genome FASTA index file. If not provided, test data will be used."
    target_regions_bed: "Optional BED file of regions to target for calling (e.g., exons, specific genomic regions)"
    is_exome: "Set to true for exome sequencing data (enables exome-specific optimizations)"
    call_germline: "Whether to perform germline variant calling (default: true)"
    call_somatic: "Whether to perform somatic variant calling (requires paired samples, default: false)"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    memory_gb: "Memory allocated for each task in the workflow in GB"
  }

  input {
    Array[StrelkaSample]? samples
    Array[StrelkaSample]? normal_samples
    File? ref_fasta
    File? ref_fasta_index
    File? target_regions_bed
    Boolean is_exome = false
    Boolean call_germline = true
    Boolean call_somatic = false
    Int cpus = 4
    Int memory_gb = 8
  }

  # Determine which genome files to use
  if (!defined(ref_fasta) || !defined(ref_fasta_index)) {
    call ww_testdata.download_ref_data { }
  }

  # Define final inputs
  File final_ref_fasta = select_first([ref_fasta, download_ref_data.fasta])
  File final_ref_fasta_index = select_first([ref_fasta_index, download_ref_data.fasta_index])

  # Download test data if necessary
  if (!defined(samples)) {
    call ww_testdata.download_bam_data as sample_data { }
  }

  # Create samples array - either from input or from test data download
  Array[StrelkaSample] final_samples = if defined(samples) then select_first([samples]) else [
    {
      "name": "demo_sample",
      "bam": select_first([sample_data.bam]),
      "bai": select_first([sample_data.bai])
    }
  ]

  # Germline variant calling
  if (call_germline) {
    scatter (sample in final_samples) {
      call strelka_germline { input:
          sample_name = sample.name,
          bam_file = sample.bam,
          bai_file = sample.bai,
          ref_fasta = final_ref_fasta,
          ref_fasta_index = final_ref_fasta_index,
          target_regions_bed = target_regions_bed,
          is_exome = is_exome,
          cpus = cpus,
          memory_gb = memory_gb
      }
    }
  }

  # Download test data if necessary
  if (!defined(normal_samples)) {
    call ww_testdata.download_bam_data as normal_data { }
  }

  # Create normals array - either from input or from test data download
  Array[StrelkaSample] final_normal_samples = if defined(normal_samples) then select_first([normal_samples]) else [
    {
      "name": "normal_sample",
      "bam": select_first([normal_data.bam]),
      "bai": select_first([normal_data.bai])
    }
  ]

  # Somatic variant calling (tumor/normal pairs)
  if (call_somatic) {
    # Ensure we have matching numbers of tumor and normal samples
    Int sample_count = length(final_samples)
    Int normal_count = length(final_normal_samples)
    
    if (sample_count == normal_count) {
      scatter (i in range(sample_count)) {
        call strelka_somatic { input:
            tumor_sample_name = final_samples[i].name,
            tumor_bam = final_samples[i].bam,
            tumor_bai = final_samples[i].bai,
            normal_sample_name = final_normal_samples[i].name,
            normal_bam = final_normal_samples[i].bam,
            normal_bai = final_normal_samples[i].bai,
            ref_fasta = final_ref_fasta,
            ref_fasta_index = final_ref_fasta_index,
            target_regions_bed = target_regions_bed,
            is_exome = is_exome,
            cpus = cpus,
            memory_gb = memory_gb
        }
      }
    }
  }

  # Validate outputs
  call validate_outputs { input:
      germline_vcfs = if call_germline && defined(strelka_germline.variants_vcf) then select_first([strelka_germline.variants_vcf]) else [],
      germline_indices = if call_germline && defined(strelka_germline.variants_vcf_index) then select_first([strelka_germline.variants_vcf_index]) else [],
      somatic_snvs_vcfs = if call_somatic && defined(strelka_somatic.somatic_snvs_vcf) then select_first([strelka_somatic.somatic_snvs_vcf]) else [],
      somatic_indels_vcfs = if call_somatic && defined(strelka_somatic.somatic_indels_vcf) then select_first([strelka_somatic.somatic_indels_vcf]) else [],
      somatic_snvs_indices = if call_somatic && defined(strelka_somatic.somatic_snvs_vcf_index) then select_first([strelka_somatic.somatic_snvs_vcf_index]) else [],
      somatic_indels_indices = if call_somatic && defined(strelka_somatic.somatic_indels_vcf_index) then select_first([strelka_somatic.somatic_indels_vcf_index]) else []
  }

  output {
    Array[File]? germline_vcfs = strelka_germline.variants_vcf
    Array[File]? germline_vcf_indices = strelka_germline.variants_vcf_index
    Array[File]? somatic_snvs_vcfs = strelka_somatic.somatic_snvs_vcf
    Array[File]? somatic_indels_vcfs = strelka_somatic.somatic_indels_vcf
    Array[File]? somatic_snvs_vcf_indices = strelka_somatic.somatic_snvs_vcf_index
    Array[File]? somatic_indels_vcf_indices = strelka_somatic.somatic_indels_vcf_index
    File validation_report = validate_outputs.report
  }
}

task strelka_germline {
  meta {
    description: "Run Strelka germline variant calling on a single sample"
    outputs: {
        variants_vcf: "Compressed VCF file containing germline variant calls",
        variants_vcf_index: "Index file for the variants VCF"
    }
  }

  parameter_meta {
    sample_name: "Name of the sample being analyzed"
    bam_file: "Input BAM file for variant calling"
    bai_file: "Index file for the input BAM"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Reference genome FASTA index file"
    target_regions_bed: "Optional BED file specifying regions to analyze"
    is_exome: "Whether this is exome sequencing data (enables exome-specific settings)"
    cpus: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    String sample_name
    File bam_file
    File bai_file
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
      --bam ~{bam_file} \
      --referenceFasta ~{ref_fasta} \
      --runDir strelka_germline_output \
      ~{exome_flag} \
      ~{regions_flag}

    # Execute the workflow
    strelka_germline_output/runWorkflow.py \
      -m local \
      -j ~{cpus}

    # Copy and rename output files
    cp strelka_germline_output/results/variants/variants.vcf.gz ~{sample_name}.strelka.germline.vcf.gz
    cp strelka_germline_output/results/variants/variants.vcf.gz.tbi ~{sample_name}.strelka.germline.vcf.gz.tbi

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
    description: "Run Strelka somatic variant calling on tumor/normal pair"
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
      --tumorBam ~{tumor_bam} \
      --normalBam ~{normal_bam} \
      --referenceFasta ~{ref_fasta} \
      --runDir strelka_somatic_output \
      ~{exome_flag} \
      ~{regions_flag}

    # Execute the workflow
    strelka_somatic_output/runWorkflow.py \
      -m local \
      -j ~{cpus}

    # Copy and rename output files
    cp strelka_somatic_output/results/variants/somatic.snvs.vcf.gz ~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.snvs.vcf.gz
    cp strelka_somatic_output/results/variants/somatic.snvs.vcf.gz.tbi ~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.snvs.vcf.gz.tbi
    cp strelka_somatic_output/results/variants/somatic.indels.vcf.gz ~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.indels.vcf.gz
    cp strelka_somatic_output/results/variants/somatic.indels.vcf.gz.tbi ~{tumor_sample_name}_vs_~{normal_sample_name}.strelka.somatic.indels.vcf.gz.tbi

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

task validate_outputs {
  meta {
    description: "Validate Strelka outputs and generate summary report"
    outputs: {
        report: "Validation summary with file checks and basic statistics"
    }
  }

  parameter_meta {
    germline_vcfs: "Array of germline VCF files to validate"
    germline_indices: "Array of germline VCF index files to validate"
    somatic_snvs_vcfs: "Array of somatic SNV VCF files to validate"
    somatic_indels_vcfs: "Array of somatic indel VCF files to validate"
    somatic_snvs_indices: "Array of somatic SNV VCF index files to validate"
    somatic_indels_indices: "Array of somatic indel VCF index files to validate"
  }

  input {
    Array[File] germline_vcfs
    Array[File] germline_indices
    Array[File] somatic_snvs_vcfs
    Array[File] somatic_indels_vcfs
    Array[File] somatic_snvs_indices
    Array[File] somatic_indels_indices
  }

  command <<<
    set -euo pipefail

    echo "=== Strelka Validation Report ===" > validation_report.txt
    echo "Generated: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    # Validate germline outputs
    if [ ~{length(germline_vcfs)} -gt 0 ]; then
      echo "Germline Variant Files:" >> validation_report.txt
      for vcf in ~{sep=' ' germline_vcfs}; do
        echo "  - $vcf" >> validation_report.txt
        variant_count=$(zcat "$vcf" | grep -v "^#" | wc -l || echo "0")
        echo "    Variants: $variant_count" >> validation_report.txt
      done
      echo "" >> validation_report.txt
    fi

    # Validate somatic outputs
    if [ ~{length(somatic_snvs_vcfs)} -gt 0 ]; then
      echo "Somatic Variant Files:" >> validation_report.txt
      for i in $(seq 0 $((~{length(somatic_snvs_vcfs)} - 1))); do
        snv_vcf=$(echo "~{sep=' ' somatic_snvs_vcfs}" | cut -d' ' -f$((i+1)))
        indel_vcf=$(echo "~{sep=' ' somatic_indels_vcfs}" | cut -d' ' -f$((i+1)))
        echo "  SNVs: $snv_vcf" >> validation_report.txt
        echo "  Indels: $indel_vcf" >> validation_report.txt
        
        snv_count=$(zcat "$snv_vcf" | grep -v "^#" | wc -l || echo "0")
        indel_count=$(zcat "$indel_vcf" | grep -v "^#" | wc -l || echo "0")
        echo "    SNV count: $snv_count" >> validation_report.txt
        echo "    Indel count: $indel_count" >> validation_report.txt
        echo "" >> validation_report.txt
      done
    fi

    # Validate file integrity
    echo "File Integrity Checks:" >> validation_report.txt
    all_files=(~{sep=' ' germline_vcfs} ~{sep=' ' germline_indices} ~{sep=' ' somatic_snvs_vcfs} ~{sep=' ' somatic_indels_vcfs} ~{sep=' ' somatic_snvs_indices} ~{sep=' ' somatic_indels_indices})
    
    for file in "${all_files[@]}"; do
      if [ -f "$file" ]; then
        size=$(stat -c%s "$file")
        echo "  $file ($size bytes)" >> validation_report.txt
      else
        echo "  Missing: $file" >> validation_report.txt
      fi
    done

    echo "" >> validation_report.txt
    echo "Validation completed successfully!" >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/strelka:2.9.10"
    memory: "2 GB"
    cpu: 1
  }
}
