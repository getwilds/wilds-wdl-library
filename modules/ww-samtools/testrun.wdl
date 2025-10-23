version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-samtools/ww-samtools.wdl" as ww_samtools

workflow samtools_example {
  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_cram_data { input:
      ref_fasta = download_ref_data.fasta
  }

  # Convert CRAM to FASTQ
  call ww_samtools.crams_to_fastq { input:
      cram_files = [download_cram_data.cram],
      ref = download_ref_data.fasta,
      name = "test_sample",
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    File r1_fastqs = crams_to_fastq.r1_fastq
    File r2_fastqs = crams_to_fastq.r2_fastq
  }
}
