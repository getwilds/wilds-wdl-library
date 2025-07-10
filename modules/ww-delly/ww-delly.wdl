## WILDS WDL module for structural variant calling using Delly.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/switch-test-data/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct DellySample {
    String name
    File bam
    File bai
}

workflow delly_example {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "WDL workflow for structural variant calling via Delly"
    url: "https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-delly"
    outputs: {
        delly_vcfs: "Structural variant calls in BCF/VCF format",
        delly_vcf_indices: "Index files for the BCF/VCF output",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    samples: "Optional list of DellySample objects, each containing sample name, BAM file, and BAM index. If not provided, workflow will download and align a demonstration sample from SRA"
    ref_fasta: "Optional reference genome FASTA file. If not provided, test data will be used."
    ref_fasta_index: "Optional reference genome FASTA index file. If not provided, test data will be used."
    target_regions_bed: "Optional BED file of regions to target for calling (e.g., exons, specific genomic regions)"
    exclude_regions_bed: "Optional BED file of regions to exclude from calling (e.g., centromeres, telomeres)"
    sv_type: "Type of structural variant to call (DEL, DUP, INV, TRA, INS, or leave empty for all types)"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    memory_gb: "Memory allocated for each task in the workflow in GB"
  }

  input {
    Array[DellySample]? samples
    File? ref_fasta
    File? ref_fasta_index
    File? target_regions_bed
    File? exclude_regions_bed
    String sv_type = ""
    Int cpus = 2
    Int memory_gb = 8
  }

  # If no reference genome provided, download test data
  if (!defined(ref_fasta) || !defined(ref_fasta_index)) {
    call ww_testdata.download_ref_data { }
  }

  # Determine which genome files to use
  File genome_fasta = select_first([ref_fasta, download_ref_data.fasta])
  File genome_fasta_index = select_first([ref_fasta_index, download_ref_data.fasta_index])

  # If no samples provided, download demonstration data
  if (!defined(samples)) {
    call ww_testdata.download_bam_data { }
  }

  # Create samples array - either from input or from test data download
  Array[DellySample] final_samples = if defined(samples) then select_first([samples]) else [
    {
      "name": "demo_sample",
      "bam": select_first([download_bam_data.bam]),
      "bai": select_first([download_bam_data.bai])
    }
  ]

  # Call structural variants using Delly for each sample
  scatter (sample in final_samples) {
    call delly_call { input:
        aligned_bam = sample.bam,
        aligned_bam_index = sample.bai,
        reference_fasta = genome_fasta,
        reference_fasta_index = genome_fasta_index,
        target_regions_bed = target_regions_bed,
        exclude_regions_bed = exclude_regions_bed,
        sv_type = sv_type,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  # Validate all outputs at the end
  call validate_outputs { input:
      delly_vcfs = delly_call.vcf,
      delly_vcf_indices = delly_call.vcf_index
  }

  output {
    Array[File] delly_vcfs = delly_call.vcf
    Array[File] delly_vcf_indices = delly_call.vcf_index
    File validation_report = validate_outputs.report
  }
}

task delly_call {
  meta {
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

task validate_outputs {
  meta {
    description: "Validates Delly outputs and generates a comprehensive report"
    outputs: {
        report: "Validation report summarizing file checks and variant statistics"
    }
  }

  parameter_meta {
    delly_vcfs: "Array of VCF files from Delly calling"
    delly_vcf_indices: "Array of index files for the VCFs"
  }

  input {
    Array[File] delly_vcfs
    Array[File] delly_vcf_indices
  }

  command <<<
    set -euo pipefail

    echo "=== Delly Workflow Validation Report ===" > validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "Total samples processed: ~{length(delly_vcfs)}" >> validation_report.txt
    echo "" >> validation_report.txt

    # Validate each VCF file
    VCF_FILES=(~{sep=" " delly_vcfs})
    INDEX_FILES=(~{sep=" " delly_vcf_indices})

    for i in "${!VCF_FILES[@]}"; do
      VCF="${VCF_FILES[$i]}"
      INDEX="${INDEX_FILES[$i]}"
      
      echo "=====================" >> validation_report.txt
      echo "- VCF file: $VCF" >> validation_report.txt
      echo "  Size: $(ls -lh "$VCF" | awk '{print $5}')" >> validation_report.txt
      echo "- VCF index: $INDEX" >> validation_report.txt
      echo "  Size: $(ls -lh "$INDEX" | awk '{print $5}')" >> validation_report.txt
      
      # Validate VCF file format
      if bcftools view -h "$VCF" > /dev/null 2>&1; then
        echo "- VCF format: VALID" >> validation_report.txt
      else
        echo "- VCF format: INVALID" >> validation_report.txt
      fi
      
      # Count variants
      VARIANT_COUNT=$(bcftools view -H "$VCF" | wc -l)
      echo "- Total variants: $VARIANT_COUNT" >> validation_report.txt
      
      echo "" >> validation_report.txt
    done

    # Overall summary
    echo "=== Overall Summary ===" >> validation_report.txt
    TOTAL_VARIANTS=0
    for vcf in "${VCF_FILES[@]}"; do
      COUNT=$(bcftools view -H "$vcf" | wc -l)
      TOTAL_VARIANTS=$((TOTAL_VARIANTS + COUNT))
    done
    echo "Total variants across all samples: $TOTAL_VARIANTS" >> validation_report.txt
    echo "" >> validation_report.txt
    echo "Validation completed successfully!" >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/delly:1.2.9"
    cpu: 1
    memory: "4GB"
  }
}
