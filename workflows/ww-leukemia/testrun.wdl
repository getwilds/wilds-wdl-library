version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/nested-vignette/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/nested-vignette/workflows/ww-leukemia/ww-leukemia.wdl" as leukemia

workflow leukemia_test {
  # Download test reference data
  call ww_testdata.download_ref_data { }

  # Download test CRAM data
  call ww_testdata.download_cram_data {
    input:
      ref_fasta = download_ref_data.fasta
  }

  # Download GATK resource files
  call ww_testdata.download_dbsnp_vcf {
    input:
      region = "chr1:1-50000000",
      filter_name = "chr1"
  }

  call ww_testdata.download_known_indels_vcf {
    input:
      region = "chr1:1-50000000",
      filter_name = "chr1"
  }

  call ww_testdata.download_gnomad_vcf {
    input:
      region = "chr1:1-50000000",
      filter_name = "chr1"
  }

  # Download ichorCNA reference data
  call ww_testdata.download_ichor_data { }

  # Call the leukemia workflow with test data
  call leukemia.ww_leukemia {
    input:
      ref_fasta = download_ref_data.fasta,
      ref_fasta_index = download_ref_data.fasta_index,
      ref_dict = download_ref_data.dict,
      dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf,
      af_only_gnomad = download_gnomad_vcf.gnomad_vcf,
      wig_gc = download_ichor_data.wig_gc,
      wig_map = download_ichor_data.wig_map,
      panel_of_norm_rds = download_ichor_data.panel_of_norm_rds,
      centromeres = download_ichor_data.centromeres,
      known_indels_sites_vcfs = [download_known_indels_vcf.known_indels_vcf],
      samples = [{
        "name": "test_sample",
        "cramfiles": [download_cram_data.cram]
      }],
      ref_name = "hg38",
      annovar_protocols = "refGene",
      annovar_operation = "g",
      scatter_count = 2  # Small scatter for testing
  }

  output {
    # Key validation outputs
    Array[File] haplotype_vcf = ww_leukemia.haplotype_vcf
    Array[File] mutect_vcf = ww_leukemia.mutect_vcf
    Array[File] mpileup_vcf = ww_leukemia.mpileup_vcf
    Array[File] consensus_variants = ww_leukemia.consensus_variants
    Array[File] manta_sv_vcf = ww_leukemia.manta_sv_vcf
    Array[File] ichorcna_genomewide_pdf = ww_leukemia.ichorcna_genomewide_pdf
  }
}
