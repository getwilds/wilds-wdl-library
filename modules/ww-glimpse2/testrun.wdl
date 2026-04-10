version 1.1

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-glimpse2/ww-glimpse2.wdl" as ww_glimpse2
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

task create_cram_directory {
  input {
    Array[File] cram_files
    Array[File] crai_files
  }

  command <<<
    set -eo pipefail
    mkdir cram_dir
    for f in ~{sep=' ' cram_files}; do
      cp "$f" cram_dir/
    done
    for f in ~{sep=' ' crai_files}; do
      cp "$f" cram_dir/
    done
  >>>

  output {
    Directory cram_directory = "cram_dir"
  }

  runtime {
    docker: "ubuntu:24.04"
    cpu: 1
    memory: "2 GB"
  }
}

workflow glimpse2_example {
  # Test region on chr1 - small enough for CI/CD but large enough for meaningful testing
  String test_region = "chr1:1-10000000"
  String test_chromosome = "chr1"
  String output_prefix = "glimpse2_test"

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

  # Step 3: Download and prepare 1000 Genomes reference panel
  call ww_testdata.download_glimpse2_reference_panel as download_reference_panel {
    input:
      chromosome = test_chromosome,
      region = test_region,
      exclude_samples = "NA12878"
  }

  # Step 4a: Download test CRAM file (NA12878) for glimpse2_phase_cram test
  call ww_testdata.download_cram_data as download_cram {
    input:
      ref_fasta = download_reference.fasta
  }

  # Step 4a-ii: Stage CRAM and index into a directory for glimpse2_phase_cram
  call create_cram_directory {
    input:
      cram_files = [download_cram.cram],
      crai_files = [download_cram.crai]
  }

  # Step 4b: Download test VCF with genotype likelihoods (NA12878 from 1000 Genomes low-coverage)
  call ww_testdata.download_glimpse2_test_gl_vcf as download_gl_vcf {
    input:
      chromosome = test_chromosome,
      region = test_region,
      sample_name = "NA12878"
  }

  # Step 4c: Download truth VCF for concordance evaluation (NA12878 high-coverage genotypes)
  call ww_testdata.download_glimpse2_truth_vcf as download_truth_vcf {
    input:
      chromosome = test_chromosome,
      region = test_region,
      sample_name = "NA12878"
  }

  # Step 5: Create chunks for the test region
  call ww_glimpse2.glimpse2_chunk {
    input:
      reference_vcf = download_reference_panel.reference_vcf,
      reference_vcf_index = download_reference_panel.reference_vcf_index,
      genetic_map = download_genetic_map.genetic_map,
      region = test_region,
      output_prefix = output_prefix,
      window_size_cm = 2.0,
      buffer_size_cm = 0.2
  }

  # Step 6: Parse chunks file to get regions for parallel processing
  call ww_glimpse2.parse_chunks_file {
    input:
      chunks_file = glimpse2_chunk.chunks_file
  }

  # Step 7: Create binary reference chunks and perform imputation for each chunk (parallel)
  scatter (idx in range(length(parse_chunks_file.input_regions))) {
    call ww_glimpse2.glimpse2_split_reference {
      input:
        reference_vcf = download_reference_panel.reference_vcf,
        reference_vcf_index = download_reference_panel.reference_vcf_index,
        genetic_map = download_genetic_map.genetic_map,
        input_region = parse_chunks_file.input_regions[idx],
        output_region = parse_chunks_file.output_regions[idx],
        output_prefix = "~{output_prefix}_chunk_~{idx}"
    }

    call ww_glimpse2.glimpse2_phase {
      input:
        input_vcf = download_gl_vcf.gl_vcf,
        input_vcf_index = download_gl_vcf.gl_vcf_index,
        reference_chunk = glimpse2_split_reference.reference_chunk,
        output_prefix = "~{output_prefix}_imputed_~{idx}"
    }

    call ww_glimpse2.glimpse2_phase_cram {
      input:
        input_cram_dir = create_cram_directory.cram_directory,
        reference_fasta = download_reference.fasta,
        reference_fasta_index = download_reference.fasta_index,
        reference_chunk = glimpse2_split_reference.reference_chunk,
        output_prefix = "~{output_prefix}_cram_imputed_~{idx}"
    }
  }

  # Step 8: Ligate all imputed chunks
  call ww_glimpse2.glimpse2_ligate {
    input:
      imputed_chunks = glimpse2_phase.imputed_chunk,
      imputed_chunks_indices = glimpse2_phase.imputed_chunk_index,
      output_prefix = "~{output_prefix}_final"
  }

  # Step 9: Evaluate imputation accuracy with concordance
  call ww_glimpse2.glimpse2_concordance {
    input:
      imputed_vcf = glimpse2_ligate.ligated_vcf,
      imputed_vcf_index = glimpse2_ligate.ligated_vcf_index,
      truth_vcf = download_truth_vcf.truth_vcf,
      truth_vcf_index = download_truth_vcf.truth_vcf_index,
      output_prefix = "~{output_prefix}_concordance"
  }

  output {
    # Test data outputs
    File reference_fasta = download_reference.fasta
    File genetic_map = download_genetic_map.genetic_map
    File reference_panel = download_reference_panel.reference_vcf
    File input_gl_vcf = download_gl_vcf.gl_vcf
    File truth_vcf = download_truth_vcf.truth_vcf

    # GLIMPSE2 intermediate outputs
    File chunks_file = glimpse2_chunk.chunks_file
    Array[File] reference_chunks = glimpse2_split_reference.reference_chunk
    Array[File] imputed_chunks = glimpse2_phase.imputed_chunk
    Array[File] cram_imputed_chunks = glimpse2_phase_cram.imputed_chunk

    # Final imputed output
    File final_imputed_vcf = glimpse2_ligate.ligated_vcf
    File final_imputed_vcf_index = glimpse2_ligate.ligated_vcf_index

    # Concordance output
    Array[File] concordance_results = glimpse2_concordance.concordance_output
  }
}
