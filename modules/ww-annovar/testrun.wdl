version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-annovar/ww-annovar.wdl" as annovar

workflow annovar_example {
  # Download test data for annotation
  call ww_testdata.download_gnomad_vcf { input:
      region = "chr1:1-10000000",
      filter_name = "chr1"
  }

  # Use test VCF data
  Array[File] vcfs_to_process = [download_gnomad_vcf.gnomad_vcf]

  scatter (vcf in vcfs_to_process) {
    call annovar.annovar_annotate { input:
        vcf_to_annotate = vcf,
        ref_name = "hg38",
        annovar_protocols = "refGene,knownGene,cosmic70,esp6500siv2_all,clinvar_20180603,gnomad211_exome",
        annovar_operation = "g,f,f,f,f,f"
    }
  }

  call annovar.validate_outputs { input:
      annotated_vcf_files = annovar_annotate.annotated_vcf,
      annotated_table_files = annovar_annotate.annotated_table
  }

  output {
    Array[File] annotated_vcf = annovar_annotate.annotated_vcf
    Array[File] annotated_table = annovar_annotate.annotated_table
    File validation_report = validate_outputs.report
  }
}
