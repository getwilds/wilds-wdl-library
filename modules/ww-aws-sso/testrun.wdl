version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-aws-sso/ww-aws-sso.wdl" as ww_aws_sso

workflow aws_sso_example {
  # Test listing bucket contents
  call ww_aws_sso.s3_list_bucket { input:
      s3_uri = "s3://gatk-test-data/wgs_fastq/",
      cpu_cores = 2,
      memory_gb = 4
  }

  # Test file download
  call ww_aws_sso.s3_download_file { input:
      s3_uri = "s3://gatk-test-data/wgs_fastq/NA12878_20k/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq",
      cpu_cores = 2,
      memory_gb = 4
  }

  # Validate all operations
  call ww_aws_sso.validate_outputs { input:
      downloaded_files = [s3_download_file.downloaded_file],
      bucket_listing = s3_list_bucket.file_list,
      object_count = s3_list_bucket.object_count
  }

  output {
    File bucket_listing = s3_list_bucket.file_list
    File downloaded_test_file = s3_download_file.downloaded_file
    File validation_report = validate_outputs.report
  }
}
