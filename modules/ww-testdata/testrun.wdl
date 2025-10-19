version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/nested-vignette/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow testdata_example {
  # Pull down reference genome and index files for chr1
  call ww_testdata.download_ref_data { input:
      chromo = "chr1",
      version = "hg38"
  }

  call ww_testdata.download_fastq_data { }

  call ww_testdata.interleave_fastq { input:
    r1_fq = download_fastq_data.r1_fastq,
    r2_fq = download_fastq_data.r2_fastq
  }

  call ww_testdata.download_cram_data { input:
    ref_fasta = download_ref_data.fasta
  }

  call ww_testdata.download_bam_data { }

  call ww_testdata.download_ichor_data { }

  call ww_testdata.download_dbsnp_vcf { input:
    region = "NC_000001.11:1-10000000",
    filter_name = "chr1"
  }

  call ww_testdata.download_known_indels_vcf { input:
    region = "chr1:1-10000000",
    filter_name = "chr1"
  }

  call ww_testdata.download_gnomad_vcf { input:
    region = "chr1:1-10000000",
    filter_name = "chr1"
  }

  call ww_testdata.download_annotsv_vcf { }

  call ww_testdata.generate_pasilla_counts { }

  call ww_testdata.validate_outputs { input:
    ref_fasta = download_ref_data.fasta,
    ref_fasta_index = download_ref_data.fasta_index,
    ref_gtf = download_ref_data.gtf,
    ref_bed = download_ref_data.bed,
    r1_fastq = download_fastq_data.r1_fastq,
    r2_fastq = download_fastq_data.r2_fastq,
    inter_fastq = interleave_fastq.inter_fastq,
    cram = download_cram_data.cram,
    crai = download_cram_data.crai,
    bam = download_bam_data.bam,
    bai = download_bam_data.bai,
    ichor_gc_wig = download_ichor_data.wig_gc,
    ichor_map_wig = download_ichor_data.wig_map,
    ichor_centromeres = download_ichor_data.centromeres,
    ichor_panel_of_norm_rds = download_ichor_data.panel_of_norm_rds,
    dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf,
    dbsnp_vcf_index = download_dbsnp_vcf.dbsnp_vcf_index,
    known_indels_vcf = download_known_indels_vcf.known_indels_vcf,
    gnomad_vcf = download_gnomad_vcf.gnomad_vcf,
    gnomad_vcf_index = download_gnomad_vcf.gnomad_vcf_index,
    annotsv_test_vcf = download_annotsv_vcf.test_vcf,
    pasilla_counts = generate_pasilla_counts.individual_count_files,
    pasilla_gene_info = generate_pasilla_counts.gene_info
  }

  output {
    # Outputs from the reference data download
    File ref_fasta = download_ref_data.fasta
    File ref_fasta_index = download_ref_data.fasta_index
    File ref_gtf = download_ref_data.gtf
    File ref_bed = download_ref_data.bed
    # Outputs from the fastq, cram, and bam data downloads
    File r1_fastq = download_fastq_data.r1_fastq
    File r2_fastq = download_fastq_data.r2_fastq
    File cram = download_cram_data.cram
    File crai = download_cram_data.crai
    File bam = download_bam_data.bam
    File bai = download_bam_data.bai
    # Outputs from the ichorCNA data download
    File ichor_gc_wig = download_ichor_data.wig_gc
    File ichor_map_wig = download_ichor_data.wig_map
    File ichor_centromeres = download_ichor_data.centromeres
    File ichor_panel_of_norm_rds = download_ichor_data.panel_of_norm_rds
    # Outputs from VCF downloads
    File dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf
    File dbsnp_vcf_index = download_dbsnp_vcf.dbsnp_vcf_index
    File known_indels_vcf = download_known_indels_vcf.known_indels_vcf
    File gnomad_vcf = download_gnomad_vcf.gnomad_vcf
    File gnomad_vcf_index = download_gnomad_vcf.gnomad_vcf_index
    File annotsv_test_vcf = download_annotsv_vcf.test_vcf
    # Outputs from Pasilla DESeq2 count generation
    Array[File] pasilla_counts = generate_pasilla_counts.individual_count_files
    Array[String] pasilla_sample_names = generate_pasilla_counts.sample_names
    Array[String] pasilla_sample_conditions = generate_pasilla_counts.sample_conditions
    File pasilla_gene_info = generate_pasilla_counts.gene_info
    # Validation report summarizing all outputs
    File validation_report = validate_outputs.report
  }
}

