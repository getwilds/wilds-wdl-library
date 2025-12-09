version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-varscan/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-varscan/modules/ww-samtools/ww-samtools.wdl" as ww_samtools

workflow samtools_example {
  # Download test data
  call ww_testdata.download_ref_data { }

  # Download two CRAM files to test merging in crams_to_fastq
  call ww_testdata.download_cram_data as download_cram_1 { input:
      ref_fasta = download_ref_data.fasta
  }
  call ww_testdata.download_cram_data as download_cram_2 { input:
      ref_fasta = download_ref_data.fasta
  }

  # Download two BAM files to test merging in merge_bams_to_cram
  call ww_testdata.download_bam_data as download_bam_1 { }
  call ww_testdata.download_bam_data as download_bam_2 { }

  # Convert multiple CRAMs to FASTQ (tests the merge functionality)
  call ww_samtools.crams_to_fastq { input:
      cram_files = [download_cram_1.cram, download_cram_2.cram],
      ref = download_ref_data.fasta,
      name = "test_sample",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Test merge_bams_to_cram task with multiple BAMs
  call ww_samtools.merge_bams_to_cram { input:
      bams_to_merge = [download_bam_1.bam, download_bam_2.bam],
      base_file_name = "test_merged",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Test mpileup task with a BAM file
  call ww_samtools.mpileup { input:
      bamfile = download_bam_1.bam,
      ref_fasta = download_ref_data.fasta,
      sample_name = "test_sample",
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    File r1_fastqs = crams_to_fastq.r1_fastq
    File r2_fastqs = crams_to_fastq.r2_fastq
    File merged_cram = merge_bams_to_cram.cram
    File pileup_file = mpileup.pileup
  }
}
