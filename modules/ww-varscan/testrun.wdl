version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-varscan/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-varscan/modules/ww-samtools/ww-samtools.wdl" as ww_samtools
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-varscan/modules/ww-varscan/ww-varscan.wdl" as ww_varscan

workflow varscan_example {
  # Download test data
  call ww_testdata.download_ref_data { }

  # Download BAM files to simulate tumor and normal samples
  call ww_testdata.download_bam_data as download_normal_bam { }
  call ww_testdata.download_bam_data as download_tumor_bam { }

  # Generate mpileup for normal sample
  call ww_samtools.mpileup as normal_mpileup {
    input:
      bamfile = download_normal_bam.bam,
      ref_fasta = download_ref_data.fasta,
      sample_name = "test_normal",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Generate mpileup for tumor sample
  call ww_samtools.mpileup as tumor_mpileup {
    input:
      bamfile = download_tumor_bam.bam,
      ref_fasta = download_ref_data.fasta,
      sample_name = "test_tumor",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Run VarScan somatic variant calling
  call ww_varscan.somatic {
    input:
      sample_name = "test_sample",
      normal_pileup = normal_mpileup.pileup,
      tumor_pileup = tumor_mpileup.pileup,
      cpu_cores = 2,
      memory_gb = 8
  }

  output {
    File somatic_snvs = somatic.somatic_snvs_vcf
    File somatic_indels = somatic.somatic_indels_vcf
  }
}
