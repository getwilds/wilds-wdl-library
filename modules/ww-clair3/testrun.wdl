version 1.0

# Import module in question as well as the testdata module for automatic demo functionality
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-clair3/modules/ww-clair3/ww-clair3.wdl" as ww_clair3
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-clair3/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

#### TEST WORKFLOW DEFINITION ####

workflow clair3_example {
  # Download reference genome (restricted to first 5 Mb for faster testing) and BAM test data
  call ww_testdata.download_ref_data as download_reference { input:
    chromo = "chr1",
    region = "1-5000000"
  }

  call ww_testdata.download_bam_data as download_bam { input:
    filename = "sample1.bam"
  }

  # Call Clair3 using Illumina model with pileup-only for speed (test data is Illumina short-read)
  call ww_clair3.run_clair3 { input:
    sample_name = "sample1",
    input_bam = download_bam.bam,
    input_bam_index = download_bam.bai,
    ref_fasta = download_reference.fasta,
    ref_fasta_index = download_reference.fasta_index,
    platform = "ilmn",
    model_path = "/opt/models/ilmn",
    bed_file = download_reference.bed,
    pileup_only = true,
    cpu_cores = 2,
    memory_gb = 8
  }

  # Test gVCF output with full-alignment mode (pileup_only doesn't support gVCF)
  call ww_clair3.run_clair3 as run_clair3_gvcf { input:
    sample_name = "sample1_gvcf",
    input_bam = download_bam.bam,
    input_bam_index = download_bam.bai,
    ref_fasta = download_reference.fasta,
    ref_fasta_index = download_reference.fasta_index,
    platform = "ilmn",
    model_path = "/opt/models/ilmn",
    bed_file = download_reference.bed,
    gvcf_enabled = true,
    cpu_cores = 2,
    memory_gb = 8
  }

  output {
    File output_vcf = run_clair3.output_vcf
    File output_vcf_index = run_clair3.output_vcf_index
    File output_gvcf_vcf = run_clair3_gvcf.output_vcf
    File output_gvcf_vcf_index = run_clair3_gvcf.output_vcf_index
    Array[File] output_gvcf = run_clair3_gvcf.output_gvcf
    Array[File] output_gvcf_index = run_clair3_gvcf.output_gvcf_index
  }
}
