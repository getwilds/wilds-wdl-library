## WILDS WDL module for GLIMPSE2 genotype imputation.
## GLIMPSE2 is a set of tools for low-coverage whole genome sequencing imputation.
## This module provides tasks for reference panel preparation, phasing, and imputation.
## Inspired by the Broad Institute's GLIMPSE Imputation Pipeline.

version 1.0

task glimpse2_chunk {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Split a genomic region into chunks for parallel imputation using GLIMPSE2_chunk"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-glimpse2/ww-glimpse2.wdl"
    outputs: {
      chunks_file: "Text file containing chunk definitions with input/output regions"
    }
  }

  parameter_meta {
    reference_vcf: "Reference panel VCF/BCF file"
    reference_vcf_index: "Index file for reference panel"
    genetic_map: "Genetic map file for the chromosome"
    region: "Genomic region to process (e.g., chr22 or chr22:16000000-20000000)"
    output_prefix: "Prefix for output files"
    window_size_cm: "Minimal window size in centiMorgans (default: 2.0)"
    buffer_size_cm: "Buffer size in centiMorgans (default: 0.2)"
    uniform_number_variants: "Use uniform number of variants per chunk instead of cM-based"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File reference_vcf
    File reference_vcf_index
    File genetic_map
    String region
    String output_prefix
    Float window_size_cm = 2.0
    Float buffer_size_cm = 0.2
    Boolean uniform_number_variants = false
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    GLIMPSE2_chunk \
      --input "~{reference_vcf}" \
      --region "~{region}" \
      --map "~{genetic_map}" \
      --window-cm ~{window_size_cm} \
      --buffer-cm ~{buffer_size_cm} \
      ~{if uniform_number_variants then "--uniform-number-variants" else ""} \
      --sequential \
      --output "~{output_prefix}.chunks.txt" \
      --threads ~{cpu_cores}
  >>>

  output {
    File chunks_file = "~{output_prefix}.chunks.txt"
  }

  runtime {
    docker: "getwilds/glimpse2:2.0.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task glimpse2_split_reference {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Convert reference panel VCF to binary format for a specific chunk using GLIMPSE2_split_reference"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-glimpse2/ww-glimpse2.wdl"
    outputs: {
      reference_chunk: "Binary reference chunk file for imputation"
    }
  }

  parameter_meta {
    reference_vcf: "Reference panel VCF/BCF file"
    reference_vcf_index: "Index file for reference panel"
    genetic_map: "Genetic map file for the chromosome"
    input_region: "Input region from chunks file (column 3)"
    output_region: "Output region from chunks file (column 4)"
    output_prefix: "Prefix for output files"
    keep_monomorphic_ref_sites: "Keep monomorphic reference sites in output"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File reference_vcf
    File reference_vcf_index
    File genetic_map
    String input_region
    String output_region
    String output_prefix
    Boolean keep_monomorphic_ref_sites = false
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    GLIMPSE2_split_reference \
      --reference "~{reference_vcf}" \
      --map "~{genetic_map}" \
      --input-region "~{input_region}" \
      --output-region "~{output_region}" \
      ~{if keep_monomorphic_ref_sites then "--keep-monomorphic-ref-sites" else ""} \
      --output "~{output_prefix}" \
      --threads ~{cpu_cores}
  >>>

  output {
    File reference_chunk = "~{output_prefix}.bin"
  }

  runtime {
    docker: "getwilds/glimpse2:2.0.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task glimpse2_phase {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Perform imputation and phasing for a single chunk using GLIMPSE2_phase"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-glimpse2/ww-glimpse2.wdl"
    outputs: {
      imputed_chunk: "Imputed and phased BCF file for the chunk",
      imputed_chunk_index: "Index file for imputed BCF"
    }
  }

  parameter_meta {
    input_vcf: "Input VCF/BCF file with genotype likelihoods (GL/PL field required)"
    input_vcf_index: "Index file for input VCF"
    reference_chunk: "Binary reference chunk from glimpse2_split_reference"
    output_prefix: "Prefix for output files"
    impute_reference_only_variants: "Impute variants only present in reference panel"
    n_burnin: "Number of burn-in iterations (default: 5)"
    n_main: "Number of main iterations (default: 15)"
    effective_population_size: "Effective population size (default: 15000)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File input_vcf
    File input_vcf_index
    File reference_chunk
    String output_prefix
    Boolean impute_reference_only_variants = false
    Int n_burnin = 5
    Int n_main = 15
    Int effective_population_size = 15000
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    GLIMPSE2_phase \
      --input "~{input_vcf}" \
      --reference "~{reference_chunk}" \
      --burnin ~{n_burnin} \
      --main ~{n_main} \
      --ne ~{effective_population_size} \
      ~{if impute_reference_only_variants then "--impute-reference-only-variants" else ""} \
      --output "~{output_prefix}.bcf" \
      --threads ~{cpu_cores}

    bcftools index "~{output_prefix}.bcf"
  >>>

  output {
    File imputed_chunk = "~{output_prefix}.bcf"
    File imputed_chunk_index = "~{output_prefix}.bcf.csi"
  }

  runtime {
    docker: "getwilds/glimpse2:2.0.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task glimpse2_phase_cram {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Perform imputation directly from CRAM/BAM files using GLIMPSE2_phase"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-glimpse2/ww-glimpse2.wdl"
    outputs: {
      imputed_chunk: "Imputed and phased BCF file for the chunk",
      imputed_chunk_index: "Index file for imputed BCF"
    }
  }

  parameter_meta {
    input_cram: "Input CRAM or BAM file"
    input_cram_index: "Index file for input CRAM/BAM"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Reference genome FASTA index file"
    reference_chunk: "Binary reference chunk from glimpse2_split_reference"
    output_prefix: "Prefix for output files"
    impute_reference_only_variants: "Impute variants only present in reference panel"
    n_burnin: "Number of burn-in iterations (default: 5)"
    n_main: "Number of main iterations (default: 15)"
    effective_population_size: "Effective population size (default: 15000)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File input_cram
    File input_cram_index
    File reference_fasta
    File reference_fasta_index
    File reference_chunk
    String output_prefix
    Boolean impute_reference_only_variants = false
    Int n_burnin = 5
    Int n_main = 15
    Int effective_population_size = 15000
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    GLIMPSE2_phase \
      --bam-file "~{input_cram}" \
      --fasta "~{reference_fasta}" \
      --reference "~{reference_chunk}" \
      --burnin ~{n_burnin} \
      --main ~{n_main} \
      --ne ~{effective_population_size} \
      ~{if impute_reference_only_variants then "--impute-reference-only-variants" else ""} \
      --output "~{output_prefix}.bcf" \
      --threads ~{cpu_cores}

    bcftools index "~{output_prefix}.bcf"
  >>>

  output {
    File imputed_chunk = "~{output_prefix}.bcf"
    File imputed_chunk_index = "~{output_prefix}.bcf.csi"
  }

  runtime {
    docker: "getwilds/glimpse2:2.0.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task glimpse2_ligate {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Ligate multiple imputed chunks into a single chromosome file using GLIMPSE2_ligate"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-glimpse2/ww-glimpse2.wdl"
    outputs: {
      ligated_vcf: "Ligated VCF/BCF file containing all imputed variants",
      ligated_vcf_index: "Index file for ligated VCF"
    }
  }

  parameter_meta {
    imputed_chunks: "Array of imputed chunk BCF files from glimpse2_phase"
    imputed_chunks_indices: "Array of index files for imputed chunks"
    output_prefix: "Prefix for output files"
    output_format: "Output format: bcf, vcf, or vcf.gz (default: bcf)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    Array[File] imputed_chunks
    Array[File] imputed_chunks_indices
    String output_prefix
    String output_format = "bcf"
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Debug: show what files we received
    echo "=== Number of imputed chunks: ~{length(imputed_chunks)} ===" >&2
    echo "=== Contents of write_lines file ===" >&2
    cat "~{write_lines(imputed_chunks)}" >&2
    echo "=== End write_lines file ===" >&2

    # Determine output extension and bcftools index type
    output_ext="~{output_format}"
    if [ "$output_ext" == "bcf" ]; then
      index_type="csi"
    else
      index_type="tbi"
    fi

    GLIMPSE2_ligate \
      --input "~{write_lines(imputed_chunks)}" \
      --output "~{output_prefix}.${output_ext}" \
      --threads ~{cpu_cores}

    # Index the output
    if [ "$output_ext" == "bcf" ]; then
      bcftools index "~{output_prefix}.${output_ext}"
    else
      bcftools index -t "~{output_prefix}.${output_ext}"
    fi
  >>>

  output {
    File ligated_vcf = "~{output_prefix}.~{output_format}"
    File ligated_vcf_index = if output_format == "bcf" then "~{output_prefix}.~{output_format}.csi" else "~{output_prefix}.~{output_format}.tbi"
  }

  runtime {
    docker: "getwilds/glimpse2:2.0.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task glimpse2_concordance {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Compute concordance metrics between imputed and truth genotypes using GLIMPSE2_concordance"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-glimpse2/ww-glimpse2.wdl"
    outputs: {
      concordance_output: "Concordance metrics output file"
    }
  }

  parameter_meta {
    imputed_vcf: "Imputed VCF/BCF file to evaluate"
    imputed_vcf_index: "Index file for imputed VCF"
    truth_vcf: "Truth/validation VCF/BCF file"
    truth_vcf_index: "Index file for truth VCF"
    allele_frequencies: "File with allele frequencies for binning"
    output_prefix: "Prefix for output files"
    region: "Genomic region to evaluate"
    min_val_dp: "Minimum depth in validation data (default: 0)"
    min_val_gq: "Minimum genotype quality in validation data (default: 0)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File imputed_vcf
    File imputed_vcf_index
    File truth_vcf
    File truth_vcf_index
    File? allele_frequencies
    String output_prefix
    String? region
    Int min_val_dp = 0
    Int min_val_gq = 0
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    GLIMPSE2_concordance \
      --input "~{imputed_vcf}" \
      --truth "~{truth_vcf}" \
      ~{if defined(allele_frequencies) then "--allele-frequencies " + allele_frequencies else ""} \
      ~{if defined(region) then "--region " + region else ""} \
      --min-val-dp ~{min_val_dp} \
      --min-val-gq ~{min_val_gq} \
      --output "~{output_prefix}" \
      --threads ~{cpu_cores}
  >>>

  output {
    Array[File] concordance_output = glob("~{output_prefix}*")
  }

  runtime {
    docker: "getwilds/glimpse2:2.0.0"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task parse_chunks_file {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Parse the chunks file to extract input and output regions for parallel processing"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-glimpse2/ww-glimpse2.wdl"
    outputs: {
      input_regions: "Array of input regions from chunks file",
      output_regions: "Array of output regions from chunks file",
      chunk_ids: "Array of chunk identifiers"
    }
  }

  parameter_meta {
    chunks_file: "Chunks file from glimpse2_chunk task"
  }

  input {
    File chunks_file
  }

  command <<<
    set -eo pipefail

    # Debug: show the chunks file content
    echo "=== Chunks file content ===" >&2
    cat "~{chunks_file}" >&2
    echo "=== End chunks file ===" >&2

    # GLIMPSE2_chunk with --sequential outputs a file with format:
    # chunk_id index input_region output_region [additional_columns...]
    # Extract input regions (column 3), output regions (column 4), and chunk IDs (columns 1-2)
    awk '{print $3}' "~{chunks_file}" > input_regions.txt
    awk '{print $4}' "~{chunks_file}" > output_regions.txt
    awk '{print $1"_"$2}' "~{chunks_file}" > chunk_ids.txt

    # Debug: show parsed content
    echo "=== Parsed input_regions ===" >&2
    cat input_regions.txt >&2
    echo "=== Parsed output_regions ===" >&2
    cat output_regions.txt >&2
    echo "=== Parsed chunk_ids ===" >&2
    cat chunk_ids.txt >&2
  >>>

  output {
    Array[String] input_regions = read_lines("input_regions.txt")
    Array[String] output_regions = read_lines("output_regions.txt")
    Array[String] chunk_ids = read_lines("chunk_ids.txt")
  }

  runtime {
    docker: "getwilds/glimpse2:2.0.0"
    cpu: 1
    memory: "2 GB"
  }
}
