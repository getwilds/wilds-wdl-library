version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/move-consensus/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-leukemia/ww-leukemia.wdl" as leukemia

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
      samples = [
        object {
          name: "test_sample",
          cramfiles: [download_cram_data.cram]
        }
      ],
      ref_name = "hg38",
      annovar_protocols = "refGene",
      annovar_operation = "g",
      # ichorCNA chromosome settings for chr1-only test data
      ichorcna_chromosomes = ["chr1"],
      ichorcna_chrs_string = "c('1')",
      # Testing-friendly resource settings
      scatter_count = 2,  # Small scatter for testing
      high_intensity_cpus = 2,  # Reduced from production default of 8
      high_intensity_memory_gb = 4,  # Reduced from production default of 16
      standard_cpus = 2,  # Reduced from production default of 4
      standard_memory_gb = 4,  # Reduced from production default of 8
      # Skip annotations to reduce disk usage in CI/CD (saves ~13GB Docker images)
      skip_annotations = true
  }

  output {
    # Key validation outputs (annotation outputs excluded since skip_annotations = true)
    Array[File] haplotype_vcf = ww_leukemia.haplotype_vcf
    Array[File] mutect_vcf = ww_leukemia.mutect_vcf
    Array[File] mpileup_vcf = ww_leukemia.mpileup_vcf
    Array[File] manta_sv_vcf = ww_leukemia.manta_sv_vcf
    Array[File] ichorcna_genomewide_pdf = ww_leukemia.ichorcna_genomewide_pdf
  }
}
