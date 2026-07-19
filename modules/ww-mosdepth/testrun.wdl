version 1.0

import "./ww-mosdepth.wdl" as ww_mosdepth
import "../ww-testdata/ww-testdata.wdl" as ww_testdata

workflow mosdepth_example {
  call ww_testdata.download_ref_data {
    input:
      chromo = "chr1",
      version = "hg38"
  }

  call ww_testdata.download_bam_data

  call ww_mosdepth.calculate_depth {
    input:
      sample_name = "demo_sample_1",
      input_bam = download_bam_data.bam,
      input_bam_index = download_bam_data.bai,
      ref_fasta = download_ref_data.fasta,
      cpu_cores = 1,
      memory_gb = 4
  }

  output {
    File depth_per_base = calculate_depth.depth_per_base
    File depth_summary = calculate_depth.depth_summary
    File region_depth = calculate_depth.region_depth
  }
}
