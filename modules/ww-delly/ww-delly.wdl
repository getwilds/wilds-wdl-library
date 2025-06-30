## WILDS WDL module for structural variant calling using Delly.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

struct SampleInfo {
    String name
    File bam
    File bai
}

struct RefGenome {
    String name
    File fasta
    File fasta_index
}

workflow delly_example {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "WDL workflow for structural variant calling via Delly"
    url: "https://github.com/getwilds/ww-delly"
    outputs: {
        delly_vcf: "Structural variant calls in BCF/VCF format",
        delly_vcf_index: "Index file for the BCF/VCF output",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    samples: "List of sample objects, each containing name, BAM file, and BAM index"
    reference_genome: "Reference genome object containing name, fasta, and fasta index files"
    exclude_regions_bed: "Optional BED file of regions to exclude from calling (e.g., centromeres, telomeres)"
    sv_type: "Type of structural variant to call (DEL, DUP, INV, TRA, INS, or leave empty for all types)"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    memory_gb: "Memory allocated for each task in the workflow in GB"
  }

  input {
    Array[SampleInfo] samples
    RefGenome reference_genome
    File? exclude_regions_bed
    String sv_type = ""
    Int cpus = 8
    Int memory_gb = 16
  }

  # Call structural variants using Delly for each sample
  scatter (sample in samples) {
    call delly_call {
      input:
        aligned_bam = sample.bam,
        aligned_bam_index = sample.bai,
        reference_fasta = reference_genome.fasta,
        reference_fasta_index = reference_genome.fasta_index,
        exclude_regions_bed = exclude_regions_bed,
        sv_type = sv_type,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  # Validate all outputs at the end
  call validate_outputs {
    input:
      delly_bcfs = delly_call.bcf,
      delly_bcf_indices = delly_call.bcf_index,
      sample_names = delly_call.sample_name
  }

  output {
    Array[File] delly_vcfs = delly_call.bcf
    Array[File] delly_vcf_indices = delly_call.bcf_index
    File validation_report = validate_outputs.report
  }
}

task delly_call {
  meta {
    description: "Calls structural variants using Delly for a single sample"
  }

  parameter_meta {
    aligned_bam: "Input aligned BAM file containing reads for variant calling"
    aligned_bam_index: "Index file for the aligned BAM above"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
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
    File? exclude_regions_bed
    String sv_type = ""
    Int cpu_cores = 8
    Int memory_gb = 16
  }

  String sample_name = basename(aligned_bam, ".bam")
  String output_prefix = "~{sample_name}_delly"
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
      -g ~{reference_fasta} \
      -o ~{output_prefix}.bcf \
      ~{aligned_bam}

    # Index the output BCF file
    bcftools index ~{output_prefix}.bcf

    # Generate summary statistics
    echo "Delly SV calling completed for sample: ~{sample_name}" > ~{output_prefix}.summary.txt
    echo "SV type filter: ~{if sv_type != "" then sv_type else "ALL"}" >> ~{output_prefix}.summary.txt
    echo "Reference genome: ~{basename(reference_fasta)}" >> ~{output_prefix}.summary.txt
    echo "Exclude regions: ~{if defined(exclude_regions_bed) then "YES" else "NO"}" >> ~{output_prefix}.summary.txt
    
    # Count variants if bcftools is available
    if command -v bcftools &> /dev/null; then
      VARIANT_COUNT=$(bcftools view -H ~{output_prefix}.bcf | wc -l)
      echo "Total variants called: $VARIANT_COUNT" >> ~{output_prefix}.summary.txt
    fi
  >>>

  output {
    File bcf = "~{output_prefix}.bcf"
    File bcf_index = "~{output_prefix}.bcf.csi"
    File summary = "~{output_prefix}.summary.txt"
    String sample_name = sample_name
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
  }

  parameter_meta {
    delly_bcfs: "Array of BCF files from Delly calling"
    delly_bcf_indices: "Array of index files for the BCFs"
    sample_names: "Array of sample names for validation"
  }

  input {
    Array[File] delly_bcfs
    Array[File] delly_bcf_indices
    Array[String] sample_names
  }

  command <<<
    set -euo pipefail

    echo "=== Delly Workflow Validation Report ===" > validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "Total samples processed: ~{length(sample_names)}" >> validation_report.txt
    echo "" >> validation_report.txt

    # Validate each BCF file
    SAMPLE_NAMES=(~{sep=' ' sample_names})
    BCF_FILES=(~{sep=' ' delly_bcfs})
    INDEX_FILES=(~{sep=' ' delly_bcf_indices})

    for i in "${!SAMPLE_NAMES[@]}"; do
      SAMPLE="${SAMPLE_NAMES[$i]}"
      BCF="${BCF_FILES[$i]}"
      INDEX="${INDEX_FILES[$i]}"
      
      echo "=== Sample: $SAMPLE ===" >> validation_report.txt
      echo "- BCF file: $BCF" >> validation_report.txt
      echo "  Size: $(ls -lh "$BCF" | awk '{print $5}')" >> validation_report.txt
      echo "- BCF index: $INDEX" >> validation_report.txt
      echo "  Size: $(ls -lh "$INDEX" | awk '{print $5}')" >> validation_report.txt
      
      # Validate BCF file format
      if bcftools view -h "$BCF" > /dev/null 2>&1; then
        echo "- BCF format: VALID" >> validation_report.txt
      else
        echo "- BCF format: INVALID" >> validation_report.txt
      fi
      
      # Count variants
      VARIANT_COUNT=$(bcftools view -H "$BCF" | wc -l)
      echo "- Total variants: $VARIANT_COUNT" >> validation_report.txt
      
      echo "" >> validation_report.txt
    done

    # Overall summary
    echo "=== Overall Summary ===" >> validation_report.txt
    TOTAL_VARIANTS=0
    for bcf in "${BCF_FILES[@]}"; do
      COUNT=$(bcftools view -H "$bcf" | wc -l)
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
