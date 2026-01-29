## WILDS WDL module for JCAST (Junction Centric Alternative Splicing Translator).
## Creates custom protein sequence databases from RNA-seq alternative splicing data
## for mass spectrometry proteomics analysis.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task jcast {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Translates alternative splicing events from rMATS output into protein sequences for proteomics analysis"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-jcast/ww-jcast.wdl"
    outputs: {
        output_fasta: "FASTA file containing translated protein sequences from alternative splicing events",
        output_directory: "Directory containing all JCAST output files"
    }
  }

  parameter_meta {
    rmats_directory: "Directory containing rMATS output files (e.g., SE.MATS.JC.txt, MXE.MATS.JC.txt, etc.)"
    gtf_file: "Ensembl GTF annotation file for the reference genome"
    genome_fasta: "Reference genome FASTA file (should be unmasked - no N characters)"
    output_name: "Prefix for output file names (default: jcast_output)"
    min_read_count: "Minimum skipped junction read count for translation (default: 1)"
    use_gmm: "Use Gaussian mixture model for read count cutoff determination (default: false)"
    write_canonical: "Write canonical protein sequences even if splice variants are untranslatable (default: false)"
    qvalue_min: "Minimum rMATS FDR q-value threshold (default: 0)"
    qvalue_max: "Maximum rMATS FDR q-value threshold (default: 1)"
    splice_types: "Comma-separated list of splice types to process (MXE,RI,SE,A3SS,A5SS). If empty, process all types."
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File rmats_directory
    File gtf_file
    File genome_fasta
    String output_name = "jcast_output"
    Int min_read_count = 1
    Boolean use_gmm = false
    Boolean write_canonical = false
    Float qvalue_min = 0
    Float qvalue_max = 1
    String splice_types = ""
    Int cpu_cores = 2
    Int memory_gb = 8
  }

  String gmm_flag = if use_gmm then "-m" else ""
  String canonical_flag = if write_canonical then "-c" else ""
  String splice_flag = if splice_types != "" then "-s " + splice_types else ""

  command <<<
    set -eo pipefail

    # Extract rMATS directory if it's a tarball
    RMATS_DIR="rmats_input"
    mkdir -p "${RMATS_DIR}"

    if [[ "~{rmats_directory}" == *.tar.gz ]] || [[ "~{rmats_directory}" == *.tgz ]]; then
      tar -xzf "~{rmats_directory}" -C "${RMATS_DIR}" --strip-components=1
    elif [[ "~{rmats_directory}" == *.tar ]]; then
      tar -xf "~{rmats_directory}" -C "${RMATS_DIR}" --strip-components=1
    elif [[ "~{rmats_directory}" == *.zip ]]; then
      unzip -q "~{rmats_directory}" -d "${RMATS_DIR}"
    else
      # Assume it's already a directory or single file
      cp -r "~{rmats_directory}" "${RMATS_DIR}/"
    fi

    echo "rMATS input directory contents:"
    ls -la "${RMATS_DIR}"

    # Run JCAST
    echo "Running JCAST..."
    jcast \
      "${RMATS_DIR}" \
      "~{gtf_file}" \
      "~{genome_fasta}" \
      -o "~{output_name}" \
      -r ~{min_read_count} \
      -q ~{qvalue_min} ~{qvalue_max} \
      ~{gmm_flag} \
      ~{canonical_flag} \
      ~{splice_flag}

    echo "JCAST completed successfully"

    # JCAST creates a nested directory structure: {output_name}/jcast_{timestamp}/*.fasta
    # Combine all protein sequence FASTA files into a single output
    find "~{output_name}" -name "*.fasta" -exec cat {} + > "~{output_name}_combined.fasta" 2>/dev/null || \
    touch "~{output_name}_combined.fasta"

    echo "Combined FASTA sequences: $(grep -c '^>' "~{output_name}_combined.fasta" || echo 0)"

    # Create output tarball of all results
    tar -czf "~{output_name}_results.tar.gz" "~{output_name}"
  >>>

  output {
    File output_fasta = "~{output_name}_combined.fasta"
    File output_directory = "~{output_name}_results.tar.gz"
  }

  runtime {
    docker: "getwilds/jcast:0.3.5"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
