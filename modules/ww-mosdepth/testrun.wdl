version 1.0

import "ww-mosdepth.wdl" as ww_mosdepth
import "../ww-testdata/ww-testdata.wdl" as ww_testdata

struct MosdepthSample {
    String name
    File bam
    File bai
}

workflow mosdepth_example {
  # Auto-download reference and BAM data
  call ww_testdata.download_ref_data {
    input:
      chromo = "chr1",
      version = "hg38"
  }
  
  call ww_testdata.download_bam_data {
    input:
      filename = "test_sample.bam"
  }

  Array[MosdepthSample] final_samples = [
    {
      "name": "demo_sample_1",
      "bam": download_bam_data.bam,
      "bai": download_bam_data.bai
    }
  ]

  scatter (sample in final_samples) {
    call ww_mosdepth.calculate_depth {
        input:
          sample_name = sample.name,
          input_bam = sample.bam,
          input_bam_index = sample.bai,
          ref_fasta = download_ref_data.fasta,
          cpu_cores = 1,
          memory_gb = 4
    }
  }

  output {
    Array[File] depth_files = calculate_depth.depth_per_base
  }
}
