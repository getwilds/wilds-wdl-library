version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-bowtie/modules/ww-bowtie/ww-bowtie.wdl" as ww_bowtie
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-bowtie/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct BowtieSample {
    String name
    File reads
    File? mates
}

workflow bowtie_example {
  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_fastq_data { }

  # Build Bowtie index from reference
  call ww_bowtie.bowtie_build { input:
      reference_fasta = download_ref_data.fasta,
      cpu_cores = 2,
      memory_gb = 8
  }

  # Create test samples array
  Array[BowtieSample] final_samples = [
    object {
      name: "demo_sample_paired",
      reads: download_fastq_data.r1_fastq,
      mates: download_fastq_data.r2_fastq
    },
    object {
      name: "demo_sample_single",
      reads: download_fastq_data.r1_fastq
    }
  ]

  # Align each sample
  scatter (sample in final_samples) {
    call ww_bowtie.bowtie_align { input:
        bowtie_index_tar = bowtie_build.bowtie_index_tar,
        reads = sample.reads,
        mates = sample.mates,
        name = sample.name,
        cpu_cores = 2,
        memory_gb = 8
    }
  }

  output {
    Array[File] aligned_bams = bowtie_align.sorted_bam
    Array[File] aligned_bais = bowtie_align.sorted_bai
  }
}
