version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-annotsv/ww-annotsv.wdl" as ww_annotsv

workflow annotsv_example {
  # Download test data for annotation
  call ww_testdata.download_annotsv_vcf { }

  # Use test VCF data
  Array[File] vcfs_to_process = [download_annotsv_vcf.test_vcf]

  scatter (vcf in vcfs_to_process) {
    call ww_annotsv.annotsv_annotate { input:
        raw_vcf = vcf,
        genome_build = "GRCh38",
        sv_min_size = 50,
        annotation_mode = "full",
        include_ci = true,
        overlap_threshold = 70,
        exclude_benign = false,
        cpu_cores = 4,
        memory_gb = 8
    }
  }

  call ww_annotsv.validate_outputs { input:
      annotated_tsv_files = annotsv_annotate.annotated_tsv
  }

  output {
    Array[File] annotated_tsv = annotsv_annotate.annotated_tsv
    File validation_report = validate_outputs.report
  }
}
