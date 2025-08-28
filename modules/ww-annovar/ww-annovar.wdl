## WILDS WDL module for variant annotation using Annovar.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow annovar_example {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow for variant annotation via Annovar"
    url: "https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-annovar"
    outputs: {
        annotated_vcf: "Annotated VCF file with variant annotations",
        annotated_table: "Tab-delimited annotation file with detailed variant information",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    vcfs: "Optional list of VCF files to annotate. If not provided, test data will be used."
    ref_name: "Reference genome name (hg19 or hg38)"
    protocols: "Comma-separated list of annotation protocols (e.g., 'refGene,gnomad211_exome,clinvar_20210123')"
    operation: "Comma-separated list of operations corresponding to protocols (e.g. 'g,f,f')"
  }

  input {
    Array[File]? vcfs
    String ref_name = "hg38"
    String protocols = "refGene,knownGene,cosmic70,esp6500siv2_all,clinvar_20180603,gnomad211_exome"
    String operation = "g,f,f,f,f,f"
  }

  # If no VCFs provided, download test data
  if (!defined(vcfs)) {
    call ww_testdata.download_gnomad_vcf { input:
        region = "chr1:1-10000000",
        filter_name = "chr1"
    }
  }

  # Determine which VCFs to use
  Array[File] vcfs_to_process = if defined(vcfs) then select_first([vcfs]) else [select_first([download_gnomad_vcf.gnomad_vcf])]

  scatter (vcf in vcfs_to_process) {
    call annovar_annotate { input:
        vcf_to_annotate = vcf,
        ref_name = ref_name,
        annovar_protocols = protocols,
        annovar_operation = operation
    }
  }

  call validate_outputs { input:
      annotated_vcf_files = annovar_annotate.annotated_vcf,
      annotated_table_files = annovar_annotate.annotated_table
  }

  output {
    Array[File] annotated_vcf = annovar_annotate.annotated_vcf
    Array[File] annotated_table = annovar_annotate.annotated_table
    File validation_report = validate_outputs.report
  }
}

task annovar_annotate {
  meta {
    description: "Annotate variants using Annovar with customizable protocols and operations"
    outputs: {
        annotated_vcf: "VCF file with Annovar annotations added",
        annotated_table: "Tab-delimited table with variant annotations"
    }
  }

  parameter_meta {
    vcf_to_annotate: "Input VCF file to be annotated"
    ref_name: "Reference genome build name for Annovar annotation (e.g., 'hg38', 'hg19')"
    annovar_protocols: "Comma-separated list of annotation protocols to apply"
    annovar_operation: "Comma-separated list of operations corresponding to the protocols"
    cpu_cores: "Number of CPU cores to allocate for the annotation task"
    memory_gb: "Memory in GB to allocate for the annotation task"
  }

  input {
    File vcf_to_annotate
    String ref_name
    String annovar_protocols
    String annovar_operation
    Int cpu_cores = 2
    Int memory_gb = 8
  }

  String base_vcf_name = basename(vcf_to_annotate, ".vcf.gz")

  command <<<
    set -eo pipefail
    perl /annovar/table_annovar.pl "~{vcf_to_annotate}" /annovar/humandb/ \
      -buildver "~{ref_name}" \
      -outfile "~{base_vcf_name}" \
      -remove \
      -protocol "~{annovar_protocols}" \
      -operation "~{annovar_operation}" \
      -nastring . -vcfinput
    sed -i "s/Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\tOtherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10\tOtherinfo11\tOtherinfo12\tOtherinfo13/Otherinfo/g" "~{base_vcf_name}.~{ref_name}_multianno.txt"
  >>>

  output {
    File annotated_vcf = "~{base_vcf_name}.~{ref_name}_multianno.vcf"
    File annotated_table = "~{base_vcf_name}.~{ref_name}_multianno.txt"
  }

  runtime {
    docker: "getwilds/annovar:~{ref_name}"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task validate_outputs {
  meta {
    description: "Validate Annovar outputs and generate summary statistics"
    outputs: {
        report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    annotated_vcf_files: "Array of annotated VCF files to validate"
    annotated_table_files: "Array of annotated table files to validate"
  }

  input {
    Array[File] annotated_vcf_files
    Array[File] annotated_table_files
  }

  command <<<
    set -euo pipefail
    
    echo "=== ANNOVAR OUTPUT VALIDATION REPORT ===" > validation_report.txt
    echo "Generated: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt
    
    # Validate VCF files
    echo "VCF Files Validation:" >> validation_report.txt
    vcf_count=0
    for vcf in ~{sep=" " annotated_vcf_files}; do
      vcf_count=$((vcf_count + 1))
      echo "  File $vcf_count: $(basename $vcf)" >> validation_report.txt
      
      # Check if file exists and is not empty
      if [ -f "$vcf" ] && [ -s "$vcf" ]; then
        echo "    Status: File exists and is not empty" >> validation_report.txt
        
        # Count variants
        variant_count=$(zcat "$vcf" | grep -v "^#" | wc -l || echo "0")
        echo "    Variants: $variant_count" >> validation_report.txt
        
        # Check VCF format
        if zcat "$vcf" | head -1 | grep -q "^##fileformat=VCF"; then
          echo "    Format: Valid VCF header" >> validation_report.txt
        else
          echo "    Format: Invalid VCF header" >> validation_report.txt
        fi
      else
        echo "    Status: File missing or empty" >> validation_report.txt
      fi
      echo "" >> validation_report.txt
    done
    
    # Validate annotation table files
    echo "Annotation Table Files Validation:" >> validation_report.txt
    table_count=0
    for table in ~{sep=" " annotated_table_files}; do
      table_count=$((table_count + 1))
      echo "  File $table_count: $(basename $table)" >> validation_report.txt
      
      if [ -f "$table" ] && [ -s "$table" ]; then
        echo "    Status: File exists and is not empty" >> validation_report.txt
        
        # Count lines (excluding header)
        line_count=$(($(wc -l < "$table") - 1))
        echo "    Annotated variants: $line_count" >> validation_report.txt
        
        # Check for key annotation columns
        if head -1 "$table" | grep -q "Chr.*Start.*End.*Ref.*Alt"; then
          echo "    Format: Contains expected annotation columns" >> validation_report.txt
        else
          echo "    Format: Missing expected annotation columns" >> validation_report.txt
        fi
      else
        echo "    Status: File missing or empty" >> validation_report.txt
      fi
      echo "" >> validation_report.txt
    done
    
    echo "=== VALIDATION SUMMARY ===" >> validation_report.txt
    echo "Total VCF files processed: $vcf_count" >> validation_report.txt
    echo "Total annotation tables generated: $table_count" >> validation_report.txt
    echo "Validation completed: $(date)" >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/bcftools:1.19"
    cpu: 1
    memory: "2 GB"
  }
}
