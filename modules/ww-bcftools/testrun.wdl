version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bcftools/ww-bcftools.wdl" as ww_bcftools

workflow bcftools_example {
  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_bam_data { }

  # Call variants on test sample
  call ww_bcftools.mpileup_call { input:
      bam_file = download_bam_data.bam,
      bam_index = download_bam_data.bai,
      reference_fasta = download_ref_data.fasta,
      reference_fasta_index = download_ref_data.fasta_index,
      max_depth = 10000,
      max_idepth = 10000,
      memory_gb = 8,
      cpu_cores = 2
  }

  call ww_bcftools.validate_outputs { input:
      vcf_files = [mpileup_call.mpileup_vcf]
  }

  output {
    File variant_vcf = mpileup_call.mpileup_vcf
    File validation_report = validate_outputs.report
  }
}
