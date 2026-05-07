version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-proseq/modules/ww-bowtie2/ww-bowtie2.wdl" as ww_bowtie2
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

  # Align each sample with default settings
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

  # Exercise the PRO-seq-style configuration on the paired sample: a sensitivity
  # preset, unaligned-read capture (as upstream rRNA depletion would use),
  # MAPQ filtering, a SAM-flag filter for proper pairs, and extra bowtie2 args.
  call ww_bowtie2.bowtie2_align as bowtie2_align_proseq { input:
      bowtie2_index_tar = bowtie2_build.bowtie2_index_tar,
      reads = download_fastq_data.r1_fastq,
      mates = download_fastq_data.r2_fastq,
      name = "demo_sample_proseq",
      preset = "very-sensitive-local",
      capture_unaligned = true,
      min_mapq = 10,
      samtools_filter_flags = "-f 2",
      extra_bowtie2_args = "--no-mixed --no-discordant",
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    Array[File] aligned_bams = bowtie2_align.sorted_bam
    Array[File] aligned_bais = bowtie2_align.sorted_bai
    File proseq_bam = bowtie2_align_proseq.sorted_bam
    File proseq_bai = bowtie2_align_proseq.sorted_bai
    File? proseq_unaligned_r1 = bowtie2_align_proseq.unaligned_r1
    File? proseq_unaligned_r2 = bowtie2_align_proseq.unaligned_r2
  }
}
