version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/bwa-vignette/vignettes/ww-bwa-gatk/ww-bwa-gatk.wdl" as bwagatk

struct BwaSample {
    String name
    File reads
    File? mates
}

workflow bwa_gatk_example {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Use test data to run the ww-bwa-gatk workflow"
    url: "https://github.com/getwilds/wilds-wdl-library/vignettes/ww-bwa-gatk"
    outputs: {
        duplicate_metrics: "Array of duplicate marking statistics for each sample",
        recalibrated_bam: "Array of BAM files with recalibrated base quality scores",
        recalibrated_bai: "Array of corresponding index files for each recalibrated bam file",
        recalibration_report: "Array of base recalibration report tables"
    }
  }

  # Download test data
  call testdata.download_ref_data {
    input:
      chromo = "chr1",
      version = "hg38"
  }

  call testdata.download_fastq_data { }

  call testdata.download_dbsnp_vcf {
    input:
      region = "NC_000001.11:1-10000000",
      filter_name = "chr1"
  }

  call testdata.download_known_indels_vcf {
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

  # Run BWA GATK workflow
  call bwagatk.bwa_gatk {
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
