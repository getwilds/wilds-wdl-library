version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/move-consensus/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/move-consensus/modules/ww-annovar/ww-annovar.wdl" as ww_annovar
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/move-consensus/modules/ww-consensus/ww-consensus.wdl" as ww_consensus

workflow consensus_example {
  # Download Annovar example VCF three times to simulate different callers
  # The same VCF is used for all three to ensure identical column structure
  # which is required by the consensus R script
  call ww_testdata.download_annovar_test_vcf as download_gatk_vcf { }
  call ww_testdata.download_annovar_test_vcf as download_bcftools_vcf { }
  call ww_testdata.download_annovar_test_vcf as download_mutect_vcf { }

  # Annotate each VCF with Annovar using hg19 (matches the example VCF coordinates)
  call ww_annovar.annovar_annotate as annotate_gatk { input:
      vcf_to_annotate = download_gatk_vcf.test_vcf,
      ref_name = "hg19",
      annovar_protocols = "refGene",
      annovar_operation = "g"
  }

  call ww_annovar.annovar_annotate as annotate_bcftools { input:
      vcf_to_annotate = download_bcftools_vcf.test_vcf,
      ref_name = "hg19",
      annovar_protocols = "refGene",
      annovar_operation = "g"
  }

  call ww_annovar.annovar_annotate as annotate_mutect { input:
      vcf_to_annotate = download_mutect_vcf.test_vcf,
      ref_name = "hg19",
      annovar_protocols = "refGene",
      annovar_operation = "g"
  }

  # Run consensus processing on the three annotated tables
  call ww_consensus.consensus_processing { input:
      gatk_vars = annotate_gatk.annotated_table,
      sam_vars = annotate_bcftools.annotated_table,
      mutect_vars = annotate_mutect.annotated_table,
      base_file_name = "test_sample",
      cpu_cores = 1,
      memory_gb = 8
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
    docker: "getwilds/consensus:0.1.1"
    memory: "4GB"
    cpu: 1
  }
}
