version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/vignettes/ww-bwa-gatk/ww-bwa-gatk.wdl" as bwa_gatk_workflow

struct BwaSample {
    String name
    File reads
    File? mates
}

workflow bwa_gatk_example {
  # Call ww_testdata tasks to get test data
  call ww_testdata.download_ref_data {
    input:
      chromo = "chr1",
      version = "hg38"
  }

  call ww_testdata.download_fastq_data { }

  call ww_testdata.download_dbsnp_vcf {
    input:
      region = "NC_000001.11:1-10000000",
      filter_name = "chr1"
  }

  call ww_testdata.download_known_indels_vcf {
    input:
      region = "chr1:1-10000000",
      filter_name = "chr1"
  }

  # Construct BwaSample struct from task outputs
  BwaSample test_sample = {
    "name": "Sample1",
    "reads": download_fastq_data.r1_fastq,
    "mates": download_fastq_data.r2_fastq
  }

  # Run actual  BWA GATK workflow
  call bwa_gatk_workflow.bwa_gatk {
    input:
      reference_fasta = download_ref_data.fasta,
      reference_fasta_index = download_ref_data.fasta_index,
      dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf,
      known_indels_vcf = [download_known_indels_vcf.known_indels_vcf],
      samples = [test_sample]
  }

  output {
    Array[File] duplicate_metrics = bwa_gatk.duplicate_metrics
    Array[File] recalibrated_bam = bwa_gatk.recalibrated_bam
    Array[File] recalibrated_bai = bwa_gatk.recalibrated_bai
    Array[File] recalibration_report = bwa_gatk.recalibration_report
  }
}

