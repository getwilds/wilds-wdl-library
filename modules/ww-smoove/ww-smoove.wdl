## WILDS WDL module for structural variant calling using Smoove.
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

workflow smoove_example {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow for structural variant calling via Smoove"
    url: "https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-smoove"
    outputs: {
        smoove_vcf: "Structural variant calls in VCF format from Smoove",
        smoove_vcf_index: "Index file for the VCF output",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    samples: "List of sample objects, each containing name, BAM file, and BAM index"
    reference_genome: "Reference genome object containing name, fasta, and fasta index files"
    exclude_bed: "Optional BED file defining regions to exclude from calling"
    include_bed: "Optional BED file defining regions to include for calling"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    memory_gb: "Memory allocated for each task in the workflow in GB"
  }

  input {
    Array[SampleInfo] samples
    RefGenome reference_genome
    File? exclude_bed
    File? include_bed
    Int cpus = 8
    Int memory_gb = 16
  }

  scatter (sample in samples) {
    call smoove_call {
      input:
        aligned_bam = sample.bam,
        aligned_bam_index = sample.bai,
        sample_name = sample.name,
        reference_fasta = reference_genome.fasta,
        reference_fasta_index = reference_genome.fasta_index,
        exclude_bed = exclude_bed,
        include_bed = include_bed,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  call validate_outputs {
    input:
      vcf_files = smoove_call.vcf,
      vcf_index_files = smoove_call.vcf_index
  }

  output {
    Array[File] smoove_vcfs = smoove_call.vcf
    Array[File] smoove_vcf_indexes = smoove_call.vcf_index
    File validation_report = validate_outputs.report
  }
}

task smoove_call {
  meta {
    description: "Call structural variants using Smoove for a single sample"
  }

  parameter_meta {
    aligned_bam: "Input aligned BAM file containing reads for variant calling"
    aligned_bam_index: "Index file for the aligned BAM above"
    sample_name: "Name of the sample provided for output files"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    exclude_bed: "Optional BED file defining regions to exclude from calling"
    include_bed: "Optional BED file defining regions to include for calling"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
    docker_image: "Docker image containing Smoove and dependencies"
  }

  input {
    File aligned_bam
    File aligned_bam_index
    String sample_name
    File reference_fasta
    File reference_fasta_index
    File? exclude_bed
    File? include_bed
    Int cpu_cores = 8
    Int memory_gb = 16
    String docker_image = "getwilds/smoove:0.2.8"
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
      --name ~{sample_name} \
      --fasta ~{reference_fasta} \
      --processes ~{cpu_cores} \
      ~{if defined(exclude_bed) then "--exclude " + exclude_bed else ""} \
      ~{if defined(include_bed) then "--include " + include_bed else ""} \
      ~{aligned_bam}

    # Move and rename output files for consistency
    mv results/~{sample_name}-smoove.genotyped.vcf.gz ~{vcf_filename}
    mv results/~{sample_name}-smoove.genotyped.vcf.gz.tbi ~{vcf_index_filename}

    # Verify outputs exist
    if [[ ! -f "~{vcf_filename}" ]]; then
      echo "ERROR: Expected VCF output file not found: ~{vcf_filename}" >&2
      exit 1
    fi

    if [[ ! -f "~{vcf_index_filename}" ]]; then
      echo "ERROR: Expected VCF index file not found: ~{vcf_index_filename}" >&2
      exit 1
    fi

    echo "Smoove completed successfully for sample: ~{sample_name}"
  >>>

  output {
    File vcf = vcf_filename
    File vcf_index = vcf_index_filename
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "${memory_gb}GB"
  }
}

task validate_outputs {
  meta {
    description: "Validate Smoove outputs and generate a comprehensive report"
  }

  parameter_meta {
    vcf_files: "Array of VCF files to validate"
    vcf_index_files: "Array of VCF index files to validate"
  }

  input {
    Array[File] vcf_files
    Array[File] vcf_index_files
  }

  command <<<
    set -euo pipefail

    echo "=== Smoove Output Validation Report ===" > validation_report.txt
    echo "Generated: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    # Validate VCF files
    echo "VCF File Validation:" >> validation_report.txt
    vcf_count=0
    for vcf in ~{sep=' ' vcf_files}; do
      vcf_count=$((vcf_count + 1))
      echo "  File $vcf_count: $(basename $vcf)" >> validation_report.txt
      
      # Check file exists and is not empty
      if [[ -f "$vcf" && -s "$vcf" ]]; then
        echo "    Status: PASS - File exists and is not empty" >> validation_report.txt
        
        # Count variants
        variant_count=$(zcat "$vcf" | grep -v "^#" | wc -l || echo "0")
        echo "    Variants: $variant_count" >> validation_report.txt
        
        # Check VCF format validity
        if zcat "$vcf" | head -n 50 | grep -q "^#CHROM"; then
          echo "    Format: PASS - Valid VCF header detected" >> validation_report.txt
        else
          echo "    Format: FAIL - Invalid VCF header" >> validation_report.txt
        fi
      else
        echo "    Status: FAIL - File missing or empty" >> validation_report.txt
      fi
      echo "" >> validation_report.txt
    done

    # Validate index files
    echo "VCF Index File Validation:" >> validation_report.txt
    index_count=0
    for index in ~{sep=' ' vcf_index_files}; do
      index_count=$((index_count + 1))
      echo "  Index $index_count: $(basename $index)" >> validation_report.txt
      
      if [[ -f "$index" && -s "$index" ]]; then
        echo "    Status: PASS - Index file exists and is not empty" >> validation_report.txt
      else
        echo "    Status: FAIL - Index file missing or empty" >> validation_report.txt
      fi
      echo "" >> validation_report.txt
    done

    # Summary
    echo "Summary:" >> validation_report.txt
    echo "  Total VCF files processed: $vcf_count" >> validation_report.txt
    echo "  Total index files processed: $index_count" >> validation_report.txt
    echo "  Validation completed: $(date)" >> validation_report.txt

    echo "Validation completed successfully"
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "ubuntu:20.04"
    cpu: 1
    memory: "2GB"
  }
}
