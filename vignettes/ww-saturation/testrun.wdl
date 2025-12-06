version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/vignettes/ww-saturation/ww-saturation.wdl" as saturation_workflow

struct SaturationSample {
    String name
    File reads
    File? mates
}

workflow saturation_mutagenesis_example {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Example workflow demonstrating saturation mutagenesis analysis using test data"
    outputs: {
        variant_counts: "Array of variant count tables for each sample",
        aa_counts: "Array of amino acid count tables for each sample",
        aa_fractions: "Array of amino acid fraction tables for each sample",
        codon_counts: "Array of codon count tables for each sample",
        codon_fractions: "Array of codon fraction tables for each sample",
        cov_length_counts: "Array of coverage length count tables for each sample",
        read_counts: "Array of read count tables for each sample",
        ref_coverage: "Array of reference coverage tables for each sample"
    }
  }

  parameter_meta {
    orf_range: "Open reading frame range to analyze (e.g., '1-100')"
    cpu_cores: "Number of CPU cores to use for processing"
    memory_gb: "Memory allocation in GB"
  }

  input {
    String orf_range = "1-99"
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  # Download test reference data (using a small region of chr1 for CI/CD memory constraints)
  # Using a 5KB region to keep memory usage under GitHub runner limits
  # Region 1000000-1005000 chosen to avoid stop codons in the test ORF
  call ww_testdata.download_ref_data {
    input:
      chromo = "chr1",
      version = "hg38",
      region = "1000000-1005000"
  }

  # Clean the reference to remove any N bases (required by GATK AnalyzeSaturationMutagenesis)
  call ww_testdata.create_clean_amplicon_reference {
    input:
      input_fasta = download_ref_data.fasta,
      output_name = "chr1_saturation_amplicon",
      replace_n_with = "A"
  }

  # Download test FASTQ data
  call ww_testdata.download_fastq_data { }

  # Construct SaturationSample struct from test data
  SaturationSample test_sample = {
    "name": "saturation_test_sample",
    "reads": download_fastq_data.r1_fastq,
    "mates": download_fastq_data.r2_fastq
  }

  # Run saturation mutagenesis workflow with cleaned reference
  call saturation_workflow.saturation_mutagenesis {
    input:
      samples = [test_sample],
      reference_fasta = create_clean_amplicon_reference.clean_fasta,
      reference_fasta_index = create_clean_amplicon_reference.clean_fasta_index,
      reference_dict = create_clean_amplicon_reference.clean_dict,
      orf_range = orf_range,
      cpu_cores = cpu_cores,
      memory_gb = memory_gb
  }

  output {
    Array[File] variant_counts = saturation_mutagenesis.variant_counts
    Array[File] aa_counts = saturation_mutagenesis.aa_counts
    Array[File] aa_fractions = saturation_mutagenesis.aa_fractions
    Array[File] codon_counts = saturation_mutagenesis.codon_counts
    Array[File] codon_fractions = saturation_mutagenesis.codon_fractions
    Array[File] cov_length_counts = saturation_mutagenesis.cov_length_counts
    Array[File] read_counts = saturation_mutagenesis.read_counts
    Array[File] ref_coverage = saturation_mutagenesis.ref_coverage
  }
}
