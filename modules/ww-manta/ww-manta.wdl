## WILDS WDL module for structural variant calling using Manta.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/switch-test-data/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct MantaSample {
    String name
    File bam
    File bai
}

workflow manta_example {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow for structural variant calling via Manta"
    url: "https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-manta"
    outputs: {
        manta_vcf: "Structural variant calls in VCF format",
        manta_vcf_index: "Index file for the VCF output",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    samples: "List of sample objects, each containing name, BAM file, and BAM index"
    ref_fasta: "Optional reference genome FASTA file. If not provided, test data will be used."
    ref_fasta_index: "Optional reference genome FASTA index file. If not provided, test data will be used."
    call_regions_bed: "Optional BED file to restrict calling to specific regions"
    call_regions_index: "Index file for the optional BED file"
    is_rna: "Boolean flag indicating if input data is RNA-seq (default: false for DNA)"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    memory_gb: "Memory allocated for each task in the workflow in GB"
  }

  input {
    Array[MantaSample]? samples
    File? ref_fasta
    File? ref_fasta_index
    File? call_regions_bed
    File? call_regions_index
    Boolean is_rna = false
    Int cpus = 2
    Int memory_gb = 8
  }

  # Determine which genome files to use
  if (!defined(ref_fasta) || !defined(ref_fasta_index)) {
    call ww_testdata.download_ref_data { }
  }
  File genome_fasta = select_first([ref_fasta, download_ref_data.fasta])
  File genome_fasta_index = select_first([ref_fasta_index, download_ref_data.fasta_index])

  # If no samples provided, download demonstration data
  if (!defined(samples)) {
    call ww_testdata.download_bam_data { }
  }

  # Create samples array - either from input or from test data download
  Array[MantaSample] final_samples = if defined(samples) then select_first([samples]) else [
    {
      "name": "demo_sample",
      "bam": select_first([download_bam_data.bam]),
      "bai": select_first([download_bam_data.bai])
    }
  ]

  scatter (sample in final_samples) {
    call manta_call { input:
        aligned_bam = sample.bam,
        aligned_bam_index = sample.bai,
        sample_name = sample.name,
        reference_fasta = genome_fasta,
        reference_fasta_index = genome_fasta_index,
        call_regions_bed = call_regions_bed,
        call_regions_index = call_regions_index,
        is_rna = is_rna,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  call validate_outputs { input:
      vcf_files = manta_call.vcf,
      vcf_index_files = manta_call.vcf_index
  }

  output {
    Array[File] manta_vcf = manta_call.vcf
    Array[File] manta_vcf_index = manta_call.vcf_index
    File validation_report = validate_outputs.report
  }
}

task manta_call {
  meta {
    description: "Call structural variants using Manta on a single sample"
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

task validate_outputs {
  meta {
    description: "Validate Manta outputs and generate a comprehensive report"
    outputs: {
        report: "Validation summary with structural variant calling statistics"
    }
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
    set -eo pipefail
    
    echo "Manta Structural Variant Calling Validation Report" > validation_report.txt
    echo "=================================================" >> validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt
    
    echo "Sample Summary:" >> validation_report.txt
    echo "Total samples processed: ~{length(vcf_files)}" >> validation_report.txt
    echo "" >> validation_report.txt
    
    # Validate each sample's outputs
    vcf_files=(~{sep=" " vcf_files})
    vcf_index_files=(~{sep=" " vcf_index_files})
    
    for i in "${!vcf_files[@]}"; do
        vcf="${vcf_files[$i]}"
        vcf_index="${vcf_index_files[$i]}"
        
        echo "Sample: $vcf" >> validation_report.txt
        
        # Check if VCF file exists and is not empty
        if [[ -f "$vcf" && -s "$vcf" ]]; then
            echo "VCF file present and non-empty" >> validation_report.txt
            variant_count=$(zcat "$vcf" | grep -v '^#' | wc -l)
            echo "Variants called: $variant_count" >> validation_report.txt
        else
            echo "VCF file missing or empty" >> validation_report.txt
        fi
        
        # Check if VCF index exists
        if [[ -f "$vcf_index" ]]; then
            echo "VCF index file present" >> validation_report.txt
        else
            echo "VCF index file missing" >> validation_report.txt
        fi
        
        echo "" >> validation_report.txt
    done
    
    echo "Validation completed successfully" >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/manta:1.6.0"
    memory: "4GB"
    cpu: 1
  }
}
