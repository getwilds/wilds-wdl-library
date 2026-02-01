## WILDS WDL module for rMATS-turbo (Replicate Multivariate Analysis of Transcript Splicing).
## Detects differential alternative splicing events from RNA-Seq data.
## Supports skipped exon (SE), alternative 5'/3' splice sites (A5SS/A3SS),
## mutually exclusive exons (MXE), and retained intron (RI) events.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task rmats {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Detects and quantifies differential alternative splicing events between two sample groups from RNA-seq BAM files"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-rmats-turbo/ww-rmats-turbo.wdl"
    outputs: {
        output_directory: "Tarball containing all rMATS output files including splice event tables and summary statistics"
    }
  }

  parameter_meta {
    gtf_file: "Gene annotation file in GTF format (required)"
    sample1_bams: "Array of BAM files for sample group 1 (condition 1)"
    sample2_bams: "Array of BAM files for sample group 2 (condition 2)"
    read_length: "The length of each read in the RNA-seq data (required)"
    read_type: "Type of reads: 'paired' for paired-end or 'single' for single-end (default: paired)"
    library_type: "Library type: 'fr-unstranded', 'fr-firststrand', or 'fr-secondstrand' (default: fr-unstranded)"
    output_name: "Prefix for output directory name (default: rmats_output)"
    variable_read_length: "Allow reads with lengths different from read_length to be processed (default: false)"
    anchor_length: "Minimum number of nucleotides mapped to each end of a splice junction (default: 1)"
    novel_splice_sites: "Enable detection of unannotated splice sites (default: false)"
    stat_off: "Skip the statistical analysis (default: false)"
    paired_stats: "Use the paired statistical model (default: false)"
    cstat: "Cutoff splicing difference for null hypothesis test (default: 0.0001)"
    individual_counts: "Output individual count files for each sample (default: false)"
    allow_clipping: "Allow alignments with soft or hard clipping (default: false)"
    min_intron_length: "Minimum intron length for novel splice site detection (default: 50)"
    max_exon_length: "Maximum exon length for novel splice site detection (default: 500)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File gtf_file
    Array[File] sample1_bams
    Array[File] sample2_bams
    Int read_length
    String read_type = "paired"
    String library_type = "fr-unstranded"
    String output_name = "rmats_output"
    Boolean variable_read_length = false
    Int anchor_length = 1
    Boolean novel_splice_sites = false
    Boolean stat_off = false
    Boolean paired_stats = false
    Float cstat = 0.0001
    Boolean individual_counts = false
    Boolean allow_clipping = false
    Int min_intron_length = 50
    Int max_exon_length = 500
    Int cpu_cores = 4
    Int memory_gb = 16
  }

  String variable_read_length_flag = if variable_read_length then "--variable-read-length" else ""
  String novel_ss_flag = if novel_splice_sites then "--novelSS" else ""
  String stat_off_flag = if stat_off then "--statoff" else ""
  String paired_stats_flag = if paired_stats then "--paired-stats" else ""
  String individual_counts_flag = if individual_counts then "--individual-counts" else ""
  String allow_clipping_flag = if allow_clipping then "--allow-clipping" else ""

  command <<<
    set -eo pipefail

    # Create output and temp directories
    mkdir -p "~{output_name}"
    mkdir -p "rmats_tmp"

    # Create BAM list files for sample groups
    echo "~{sep=',' sample1_bams}" > b1.txt
    echo "~{sep=',' sample2_bams}" > b2.txt

    echo "Sample 1 BAM files:"
    cat b1.txt
    echo ""
    echo "Sample 2 BAM files:"
    cat b2.txt
    echo ""

    # Run rMATS-turbo
    echo "Running rMATS-turbo..."
    python /rmats/rmats.py \
      --gtf "~{gtf_file}" \
      --b1 b1.txt \
      --b2 b2.txt \
      --od "~{output_name}" \
      --tmp rmats_tmp \
      -t ~{read_type} \
      --libType ~{library_type} \
      --readLength ~{read_length} \
      --anchorLength ~{anchor_length} \
      --nthread ~{cpu_cores} \
      --cstat ~{cstat} \
      --mil ~{min_intron_length} \
      --mel ~{max_exon_length} \
      ~{variable_read_length_flag} \
      ~{novel_ss_flag} \
      ~{stat_off_flag} \
      ~{paired_stats_flag} \
      ~{individual_counts_flag} \
      ~{allow_clipping_flag}

    echo "rMATS-turbo completed successfully"

    # List outputs
    echo "Output directory contents:"
    ls -la "~{output_name}"

    # Create output tarball of all results
    tar -czf "~{output_name}.tar.gz" "~{output_name}"
  >>>

  output {
    File output_directory = "~{output_name}.tar.gz"
  }

  runtime {
    docker: "getwilds/rmats-turbo:latest"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task rmats_prep {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Preprocesses BAM files and generates .rmats intermediate files (prep step only)"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-rmats-turbo/ww-rmats-turbo.wdl"
    outputs: {
        prep_output: "Tarball containing .rmats intermediate files from preprocessing"
    }
  }

  parameter_meta {
    gtf_file: "Gene annotation file in GTF format (required)"
    sample_bams: "Array of BAM files to preprocess"
    read_length: "The length of each read in the RNA-seq data (required)"
    read_type: "Type of reads: 'paired' for paired-end or 'single' for single-end (default: paired)"
    library_type: "Library type: 'fr-unstranded', 'fr-firststrand', or 'fr-secondstrand' (default: fr-unstranded)"
    output_name: "Prefix for output directory name (default: rmats_prep)"
    variable_read_length: "Allow reads with lengths different from read_length to be processed (default: false)"
    anchor_length: "Minimum number of nucleotides mapped to each end of a splice junction (default: 1)"
    novel_splice_sites: "Enable detection of unannotated splice sites (default: false)"
    allow_clipping: "Allow alignments with soft or hard clipping (default: false)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File gtf_file
    Array[File] sample_bams
    Int read_length
    String read_type = "paired"
    String library_type = "fr-unstranded"
    String output_name = "rmats_prep"
    Boolean variable_read_length = false
    Int anchor_length = 1
    Boolean novel_splice_sites = false
    Boolean allow_clipping = false
    Int cpu_cores = 4
    Int memory_gb = 16
  }

  String variable_read_length_flag = if variable_read_length then "--variable-read-length" else ""
  String novel_ss_flag = if novel_splice_sites then "--novelSS" else ""
  String allow_clipping_flag = if allow_clipping then "--allow-clipping" else ""

  command <<<
    set -eo pipefail

    # Create output and temp directories
    mkdir -p "~{output_name}"
    mkdir -p "rmats_tmp"

    # Create BAM list file
    echo "~{sep=',' sample_bams}" > bams.txt

    echo "BAM files to preprocess:"
    cat bams.txt
    echo ""

    # Run rMATS-turbo prep step
    echo "Running rMATS-turbo prep step..."
    python /rmats/rmats.py \
      --gtf "~{gtf_file}" \
      --b1 bams.txt \
      --od "~{output_name}" \
      --tmp rmats_tmp \
      -t ~{read_type} \
      --libType ~{library_type} \
      --readLength ~{read_length} \
      --anchorLength ~{anchor_length} \
      --nthread ~{cpu_cores} \
      --task prep \
      ~{variable_read_length_flag} \
      ~{novel_ss_flag} \
      ~{allow_clipping_flag}

    echo "rMATS-turbo prep step completed successfully"

    # List outputs
    echo "Prep output directory contents:"
    ls -la rmats_tmp

    # Create output tarball of prep results
    tar -czf "~{output_name}.tar.gz" rmats_tmp
  >>>

  output {
    File prep_output = "~{output_name}.tar.gz"
  }

  runtime {
    docker: "getwilds/rmats-turbo:latest"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task rmats_post {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Loads .rmats files and performs alternative splicing event detection and statistical analysis (post step)"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-rmats-turbo/ww-rmats-turbo.wdl"
    outputs: {
        output_directory: "Tarball containing all rMATS output files including splice event tables and summary statistics"
    }
  }

  parameter_meta {
    gtf_file: "Gene annotation file in GTF format (required)"
    prep_outputs: "Array of tarballs containing .rmats files from prep step"
    read_length: "The length of each read in the RNA-seq data (required)"
    read_type: "Type of reads: 'paired' for paired-end or 'single' for single-end (default: paired)"
    output_name: "Prefix for output directory name (default: rmats_output)"
    stat_off: "Skip the statistical analysis (default: false)"
    paired_stats: "Use the paired statistical model (default: false)"
    cstat: "Cutoff splicing difference for null hypothesis test (default: 0.0001)"
    individual_counts: "Output individual count files for each sample (default: false)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File gtf_file
    Array[File] prep_outputs
    Int read_length
    String read_type = "paired"
    String output_name = "rmats_output"
    Boolean stat_off = false
    Boolean paired_stats = false
    Float cstat = 0.0001
    Boolean individual_counts = false
    Int cpu_cores = 4
    Int memory_gb = 16
  }

  String stat_off_flag = if stat_off then "--statoff" else ""
  String paired_stats_flag = if paired_stats then "--paired-stats" else ""
  String individual_counts_flag = if individual_counts then "--individual-counts" else ""

  command <<<
    set -eo pipefail

    # Create output and temp directories
    mkdir -p "~{output_name}"
    mkdir -p "rmats_tmp"

    # Extract all prep outputs into the temp directory
    echo "Extracting prep outputs..."
    for prep_file in ~{sep=' ' prep_outputs}; do
      tar -xzf "${prep_file}" -C . --strip-components=0
    done

    echo "Extracted prep files:"
    ls -la rmats_tmp

    # Run rMATS-turbo post step
    echo "Running rMATS-turbo post step..."
    python /rmats/rmats.py \
      --gtf "~{gtf_file}" \
      --od "~{output_name}" \
      --tmp rmats_tmp \
      -t ~{read_type} \
      --readLength ~{read_length} \
      --nthread ~{cpu_cores} \
      --cstat ~{cstat} \
      --task post \
      ~{stat_off_flag} \
      ~{paired_stats_flag} \
      ~{individual_counts_flag}

    echo "rMATS-turbo post step completed successfully"

    # List outputs
    echo "Output directory contents:"
    ls -la "~{output_name}"

    # Create output tarball of all results
    tar -czf "~{output_name}.tar.gz" "~{output_name}"
  >>>

  output {
    File output_directory = "~{output_name}.tar.gz"
  }

  runtime {
    docker: "getwilds/rmats-turbo:latest"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task rmats_stat {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Runs statistical analysis on existing rMATS output files"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-rmats-turbo/ww-rmats-turbo.wdl"
    outputs: {
        output_directory: "Tarball containing rMATS output files with updated statistical results"
    }
  }

  parameter_meta {
    gtf_file: "Gene annotation file in GTF format (required)"
    existing_output: "Tarball containing existing rMATS output files to run statistics on"
    read_length: "The length of each read in the RNA-seq data (required)"
    read_type: "Type of reads: 'paired' for paired-end or 'single' for single-end (default: paired)"
    output_name: "Prefix for output directory name (default: rmats_output)"
    paired_stats: "Use the paired statistical model (default: false)"
    cstat: "Cutoff splicing difference for null hypothesis test (default: 0.0001)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File gtf_file
    File existing_output
    Int read_length
    String read_type = "paired"
    String output_name = "rmats_output"
    Boolean paired_stats = false
    Float cstat = 0.0001
    Int cpu_cores = 4
    Int memory_gb = 16
  }

  String paired_stats_flag = if paired_stats then "--paired-stats" else ""

  command <<<
    set -eo pipefail

    # Create temp directory
    mkdir -p "rmats_tmp"

    # Extract existing output
    echo "Extracting existing rMATS output..."
    tar -xzf "~{existing_output}" -C . --strip-components=0

    # Find the extracted directory
    EXISTING_DIR=$(find . -maxdepth 1 -type d -name "*rmats*" -o -name "*output*" | head -1)
    if [ -z "$EXISTING_DIR" ]; then
      EXISTING_DIR="."
    fi

    echo "Using existing output directory: ${EXISTING_DIR}"
    ls -la "${EXISTING_DIR}"

    # Run rMATS-turbo stat step
    echo "Running rMATS-turbo stat step..."
    python /rmats/rmats.py \
      --gtf "~{gtf_file}" \
      --od "${EXISTING_DIR}" \
      --tmp rmats_tmp \
      -t ~{read_type} \
      --readLength ~{read_length} \
      --nthread ~{cpu_cores} \
      --cstat ~{cstat} \
      --task stat \
      ~{paired_stats_flag}

    echo "rMATS-turbo stat step completed successfully"

    # List outputs
    echo "Output directory contents:"
    ls -la "${EXISTING_DIR}"

    # Create output tarball of all results
    tar -czf "~{output_name}.tar.gz" "${EXISTING_DIR}"
  >>>

  output {
    File output_directory = "~{output_name}.tar.gz"
  }

  runtime {
    docker: "getwilds/rmats-turbo:latest"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
