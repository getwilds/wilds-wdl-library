version 1.0

import "ww-deepvariant.wdl" as ww_deepvariant
import "../ww-testdata/ww-testdata.wdl" as ww_testdata

struct DeepVariantSample {
    String name
    File bam
    File bai
}

workflow deepvariant_example {
  # Download reference genome and BAM test data
  call ww_testdata.download_ref_data as download_reference { input:
    chromo = "chr1"
  }

  call ww_testdata.download_bam_data as download_bam_1 { input:
    filename = "sample1.bam"
  }

  call ww_testdata.download_bam_data as download_bam_2 { input:
    filename = "sample2.bam"
  }

  # Create samples array using test data
  Array[DeepVariantSample] samples = [
    {
      "name": "sample1",
      "bam": download_bam_1.bam,
      "bai": download_bam_1.bai
    },
    {
      "name": "sample2",
      "bam": download_bam_2.bam,
      "bai": download_bam_2.bai
    }
  ]

  # Call DeepVariant on each sample
  scatter (sample in samples) {
    call ww_deepvariant.run_deepvariant { input:
      sample_name = sample.name,
      input_bam = sample.bam,
      input_bam_index = sample.bai,
      ref_fasta = download_reference.fasta,
      ref_fasta_index = download_reference.fasta_index,
      model_type = "WGS",
      regions = "chr1",
      cpu_cores = 2,
      memory_gb = 8
    }
  }

  output {
    Array[File] output_vcfs = run_deepvariant.output_vcf
    Array[File] output_vcf_indices = run_deepvariant.output_vcf_index
  }
}
