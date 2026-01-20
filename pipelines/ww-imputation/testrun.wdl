version 1.0

# Import the imputation pipeline and testdata module for automatic test data provisioning
import "../../modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "ww-imputation.wdl" as ww_imputation

#### TEST WORKFLOW DEFINITION ####
# This workflow demonstrates the ww-imputation pipeline with automatic test data download.
# It downloads a small region of chr1 data from 1000 Genomes and performs imputation
# on a single sample from CRAM input.

workflow imputation_testrun {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Zero-configuration test workflow for ww-imputation. Downloads test data and runs imputation on a small chr1 region."
  }

  # Test region on chr1 - small enough for CI/CD but large enough for meaningful testing
  String test_region = "chr1:1-10000000"
  String test_chromosome = "chr1"

  # Step 1: Download reference genome for chr1
  call ww_testdata.download_ref_data as download_reference {
    input:
      chromo = test_chromosome,
      version = "hg38",
      region = "1-10000000",
      output_name = "chr1_test"
  }

  # Step 2: Download genetic map for chr1
  call ww_testdata.download_glimpse2_genetic_map as download_genetic_map {
    input:
      chromosome = test_chromosome,
      genome_build = "b38"
  }

  # Step 3: Download 1000 Genomes reference panel (excluding NA12878 for validation)
  call ww_testdata.download_glimpse2_reference_panel as download_reference_panel {
    input:
      chromosome = test_chromosome,
      region = test_region,
      exclude_samples = "NA12878"
  }

  # Step 4: Download test CRAM file (NA12878 from existing testdata task)
  call ww_testdata.download_cram_data as download_cram {
    input:
      ref_fasta = download_reference.fasta
  }

  # Step 5: Run the imputation pipeline
  call ww_imputation.imputation {
    input:
      samples = [
        {
          "sample_id": "NA12878",
          "cram": download_cram.cram,
          "cram_index": download_cram.crai
        }
      ],
      chromosomes = [
        {
          "chromosome": test_chromosome,
          "reference_vcf": download_reference_panel.reference_vcf,
          "reference_vcf_index": download_reference_panel.reference_vcf_index,
          "genetic_map": download_genetic_map.genetic_map
        }
      ],
      reference_fasta = download_reference.fasta,
      reference_fasta_index = download_reference.fasta_index,
      output_format = "bcf",
      # Use smaller resources for testing
      chunk_cpu_cores = 2,
      chunk_memory_gb = 4,
      phase_cpu_cores = 2,
      phase_memory_gb = 4,
      ligate_cpu_cores = 2,
      ligate_memory_gb = 4
  }

  output {
    # Test data outputs
    File reference_fasta = download_reference.fasta
    File genetic_map = download_genetic_map.genetic_map
    File reference_panel = download_reference_panel.reference_vcf
    File input_cram = download_cram.cram

    # Final imputed outputs
    Array[File] imputed_vcfs = imputation.imputed_vcfs
    Array[File] imputed_vcf_indices = imputation.imputed_vcf_indices
  }
}
