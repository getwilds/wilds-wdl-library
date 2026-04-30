version 1.0

# Local relative imports — flip to GitHub raw URLs at merge time.
import "ww-umi-tools.wdl" as ww_umi_tools
import "../ww-testdata/ww-testdata.wdl" as ww_testdata

struct DedupSample {
  String name
  File bam
  File bai
}

#### TEST WORKFLOW DEFINITION ####
# Demonstrates umi_tools dedup on a small paired-end BAM. The test data BAM from
# ww-testdata has no UMIs in its read names, so we inject synthetic UMIs via the
# ww-testdata helper task before running dedup.

workflow umi_tools_example {
  call ww_testdata.download_bam_data as download_demo_bam { }

  call ww_testdata.inject_synthetic_umis as add_umis { input:
      input_bam = download_demo_bam.bam,
      input_bai = download_demo_bam.bai,
      sample_name = "demo_sample"
  }

  Array[DedupSample] samples = [
    {
      "name": "demo_sample",
      "bam": add_umis.umi_bam,
      "bai": add_umis.umi_bai
    }
  ]

  scatter (sample in samples) {
    call ww_umi_tools.dedup { input:
        input_bam = sample.bam,
        input_bai = sample.bai,
        sample_name = sample.name,
        paired = true,
        umi_separator = ":",
        cpu_cores = 2,
        memory_gb = 4
    }
  }

  output {
    Array[File] deduped_bams = dedup.deduped_bam
    Array[File] deduped_bais = dedup.deduped_bai
    Array[File] dedup_logs = dedup.log
  }
}
