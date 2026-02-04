version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-annovar/ww-annovar.wdl" as ww_annovar
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-consensus/ww-consensus.wdl" as ww_consensus

workflow consensus_example {
  # Download a small gnomAD VCF subset to use as test data
  call ww_testdata.download_gnomad_vcf { input:
      region = "chr1:1-1000000",
      filter_name = "chr1_1Mb"
  }

  # Annotate the VCF three times (simulating output from three different callers)
  # In a real workflow, these would be VCFs from GATK, bcftools, and Mutect2
  call ww_annovar.annovar_annotate as annotate_gatk { input:
      vcf_to_annotate = download_gnomad_vcf.gnomad_vcf,
      ref_name = "hg38",
      annovar_protocols = "refGene",
      annovar_operation = "g"
  }

  call ww_annovar.annovar_annotate as annotate_bcftools { input:
      vcf_to_annotate = download_gnomad_vcf.gnomad_vcf,
      ref_name = "hg38",
      annovar_protocols = "refGene",
      annovar_operation = "g"
  }

  call ww_annovar.annovar_annotate as annotate_mutect { input:
      vcf_to_annotate = download_gnomad_vcf.gnomad_vcf,
      ref_name = "hg38",
      annovar_protocols = "refGene",
      annovar_operation = "g"
  }

  # Run consensus processing on the three annotated tables
  call ww_consensus.consensus_processing { input:
      haplo_vars = annotate_gatk.annotated_table,
      mpileup_vars = annotate_bcftools.annotated_table,
      mutect_vars = annotate_mutect.annotated_table,
      base_file_name = "test_sample"
  }

  call validate_outputs { input:
      consensus_file = consensus_processing.consensus_tsv
  }

  output {
    File gatk_annotated = annotate_gatk.annotated_table
    File bcftools_annotated = annotate_bcftools.annotated_table
    File mutect_annotated = annotate_mutect.annotated_table
    File consensus_tsv = consensus_processing.consensus_tsv
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  meta {
    description: "Validate consensus processing outputs and generate summary statistics"
    outputs: {
        report: "Validation summary with consensus variant statistics"
    }
  }

  parameter_meta {
    consensus_file: "Consensus TSV file to validate"
  }

  input {
    File consensus_file
  }

  command <<<
    set -eo pipefail

    echo "Consensus Variant Validation Report" > validation_report.txt
    echo "====================================" >> validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    # Validate consensus file
    echo "Consensus File Validation:" >> validation_report.txt
    echo "--------------------------" >> validation_report.txt
    echo "File: $(basename ~{consensus_file})" >> validation_report.txt

    # Check if file exists and is not empty
    if [ -f "~{consensus_file}" ] && [ -s "~{consensus_file}" ]; then
      echo "  File exists and is not empty" >> validation_report.txt

      # Count variants
      VARIANT_COUNT=$(wc -l < "~{consensus_file}")
      # Subtract 1 for header if present
      VARIANT_COUNT=$((VARIANT_COUNT - 1))
      echo "  Total consensus variants: $VARIANT_COUNT" >> validation_report.txt

      # Count by confidence tier
      echo "" >> validation_report.txt
      echo "Confidence Tier Breakdown:" >> validation_report.txt
      HIGH=$(grep -c "HighConfidence" "~{consensus_file}" || echo "0")
      MED=$(grep -c "MediumConfidence" "~{consensus_file}" || echo "0")
      LOW=$(grep -c "LowConfidence" "~{consensus_file}" || echo "0")
      echo "  HighConfidence (3 callers): $HIGH" >> validation_report.txt
      echo "  MediumConfidence (2 callers): $MED" >> validation_report.txt
      echo "  LowConfidence (1 caller): $LOW" >> validation_report.txt
    else
      echo "  WARNING: File missing or empty" >> validation_report.txt
    fi

    echo "" >> validation_report.txt
    echo "Validation completed successfully." >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "ubuntu:22.04"
    memory: "4GB"
    cpu: 1
  }
}
