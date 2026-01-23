## WILDS WDL pipeline for genotype imputation using GLIMPSE2.
## This pipeline performs low-coverage whole genome sequencing imputation
## from CRAM/BAM files using a reference panel.
##
## Inspired by the Broad Institute's GLIMPSE Imputation Pipeline:
## https://github.com/broadinstitute/palantir-workflows/tree/main/GlimpseImputationPipeline

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-glimpse2/modules/ww-glimpse2/ww-glimpse2.wdl" as glimpse2
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-glimpse2/modules/ww-bcftools/ww-bcftools.wdl" as bcftools

struct ImputationSample {
    String sample_id
    File cram
    File cram_index
}

struct ChromosomeData {
    String chromosome
    File reference_vcf
    File reference_vcf_index
    File genetic_map
}

workflow imputation {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow for genotype imputation from low-coverage WGS data using GLIMPSE2. Processes CRAM/BAM files against a reference panel to produce imputed VCF files."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/pipelines/ww-imputation/ww-imputation.wdl"
    outputs: {
        imputed_vcfs: "Array of imputed VCF/BCF files, one per sample (all chromosomes ligated together)",
        imputed_vcf_indices: "Array of index files for imputed VCFs"
    }
  }

  parameter_meta {
    samples: "Array of ImputationSample objects containing sample_id, CRAM file, and CRAM index"
    chromosomes: "Array of ChromosomeData objects containing chromosome name, reference panel VCF, index, and genetic map"
    reference_fasta: "Reference genome FASTA file (must match CRAM reference)"
    reference_fasta_index: "Reference genome FASTA index file (.fai)"
    output_format: "Output format for imputed files: bcf, vcf, or vcf.gz (default: bcf)"
    impute_reference_only_variants: "Only impute variants present in reference panel (default: false)"
    window_size_cm: "Chunk window size in centiMorgans (default: 2.0)"
    buffer_size_cm: "Chunk buffer size in centiMorgans (default: 0.2)"
    n_burnin: "Number of burn-in Markov chain Monte Carlo iterations (default: 5)"
    n_main: "Number of main Markov chain Monte Carlo iterations (default: 15)"
    effective_population_size: "Effective population size for hidden Markov model (default: 15000)"
    chunk_cpu_cores: "CPU cores for chunking tasks"
    chunk_memory_gb: "Memory in GB for chunking tasks"
    phase_cpu_cores: "CPU cores for phasing/imputation tasks"
    phase_memory_gb: "Memory in GB for phasing/imputation tasks"
    ligate_cpu_cores: "CPU cores for ligation tasks"
    ligate_memory_gb: "Memory in GB for ligation tasks"
  }

  input {
    # Required inputs
    Array[ImputationSample] samples
    Array[ChromosomeData] chromosomes
    File reference_fasta
    File reference_fasta_index

    # Output options
    String output_format = "bcf"

    # Imputation parameters
    Boolean impute_reference_only_variants = false
    Float window_size_cm = 2.0
    Float buffer_size_cm = 0.2
    Int n_burnin = 5
    Int n_main = 15
    Int effective_population_size = 15000

    # Resource allocation
    Int chunk_cpu_cores = 4
    Int chunk_memory_gb = 8
    Int phase_cpu_cores = 4
    Int phase_memory_gb = 8
    Int ligate_cpu_cores = 4
    Int ligate_memory_gb = 16
  }

  # Step 1: Process each chromosome - chunk and split reference
  scatter (chrom_data in chromosomes) {
    # Generate chunks for this chromosome
    call glimpse2.glimpse2_chunk {
      input:
        reference_vcf = chrom_data.reference_vcf,
        reference_vcf_index = chrom_data.reference_vcf_index,
        genetic_map = chrom_data.genetic_map,
        region = chrom_data.chromosome,
        output_prefix = chrom_data.chromosome,
        window_size_cm = window_size_cm,
        buffer_size_cm = buffer_size_cm,
        cpu_cores = chunk_cpu_cores,
        memory_gb = chunk_memory_gb
    }

    # Parse the chunks file to get regions
    call glimpse2.parse_chunks_file {
      input:
        chunks_file = glimpse2_chunk.chunks_file
    }

    # Split reference into binary chunks for this chromosome
    scatter (idx in range(length(parse_chunks_file.input_regions))) {
      call glimpse2.glimpse2_split_reference {
        input:
          reference_vcf = chrom_data.reference_vcf,
          reference_vcf_index = chrom_data.reference_vcf_index,
          genetic_map = chrom_data.genetic_map,
          input_region = parse_chunks_file.input_regions[idx],
          output_region = parse_chunks_file.output_regions[idx],
          output_prefix = chrom_data.chromosome + "_chunk_" + parse_chunks_file.chunk_ids[idx],
          cpu_cores = chunk_cpu_cores,
          memory_gb = chunk_memory_gb
      }
    }

    # Collect all reference chunks for this chromosome
    Array[File] chrom_reference_chunks = glimpse2_split_reference.reference_chunk
  }

  # Step 2: Process each sample across all chromosomes
  scatter (sample in samples) {
    # Process each chromosome for this sample
    scatter (chrom_idx in range(length(chromosomes))) {
      # Phase each chunk for this chromosome
      scatter (chunk_idx in range(length(chrom_reference_chunks[chrom_idx]))) {
        call glimpse2.glimpse2_phase_cram as phase_chunk {
          input:
            input_cram = sample.cram,
            input_cram_index = sample.cram_index,
            reference_fasta = reference_fasta,
            reference_fasta_index = reference_fasta_index,
            reference_chunk = chrom_reference_chunks[chrom_idx][chunk_idx],
            output_prefix = "~{sample.sample_id}_~{chromosomes[chrom_idx].chromosome}_chunk_~{chunk_idx}",
            impute_reference_only_variants = impute_reference_only_variants,
            n_burnin = n_burnin,
            n_main = n_main,
            effective_population_size = effective_population_size,
            cpu_cores = phase_cpu_cores,
            memory_gb = phase_memory_gb
        }
      }

      # Ligate chunks for this chromosome
      call glimpse2.glimpse2_ligate as ligate_chromosome {
        input:
          imputed_chunks = phase_chunk.imputed_chunk,
          imputed_chunks_indices = phase_chunk.imputed_chunk_index,
          output_prefix = sample.sample_id + "_" + chromosomes[chrom_idx].chromosome + "_imputed",
          output_format = output_format,
          cpu_cores = ligate_cpu_cores,
          memory_gb = ligate_memory_gb
      }
    }

    # Concatenate all chromosomes into a single file per sample
    call bcftools.concat as concatenate_chromosomes {
      input:
        vcf_files = ligate_chromosome.ligated_vcf,
        vcf_indices = ligate_chromosome.ligated_vcf_index,
        output_prefix = sample.sample_id + "_imputed_all_chromosomes",
        output_format = output_format
    }
  }

  output {
    Array[File] imputed_vcfs = concatenate_chromosomes.concatenated_vcf
    Array[File] imputed_vcf_indices = concatenate_chromosomes.concatenated_vcf_index
  }
}
