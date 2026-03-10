version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bowtie2/ww-bowtie2.wdl" as ww_bowtie2
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct Bowtie2Sample {
    String name
    File reads
    File? mates
}

workflow bowtie2_example {
  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_fastq_data { }

  # Build Bowtie 2 index from reference
  call ww_bowtie2.bowtie2_build { input:
      reference_fasta = download_ref_data.fasta,
      cpu_cores = 2,
      memory_gb = 8
  }

  # Create test samples array
  Array[Bowtie2Sample] final_samples = [
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
    call ww_bowtie2.bowtie2_align { input:
        bowtie2_index_tar = bowtie2_build.bowtie2_index_tar,
        reads = sample.reads,
        mates = sample.mates,
        name = sample.name,
        cpu_cores = 2,
        memory_gb = 8
    }
  }

  output {
    Array[File] aligned_bams = bowtie2_align.sorted_bam
    Array[File] aligned_bais = bowtie2_align.sorted_bai
  }
}
