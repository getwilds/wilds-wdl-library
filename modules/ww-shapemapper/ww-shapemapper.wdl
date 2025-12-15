## WILDS WDL ShapeMapper Module
## ShapeMapper analyzes RNA structure probing data to determine nucleotide flexibility
## and secondary structure. This module wraps ShapeMapper 2 for processing SHAPE-MaP
## and related chemical probing datasets.

version 1.0

task run_shapemapper {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Run ShapeMapper to analyze RNA structure probing data and generate reactivity profiles"
    outputs: {
        shape_file: "SHAPE reactivity profile with per-nucleotide reactivity scores",
        log_file: "ShapeMapper log file with processing details and quality metrics"
    }
  }

  parameter_meta {
    sample_name: "Name identifier for the sample"
    target_fa: "FASTA file containing the target RNA sequence(s)"
    modified_r1: "R1 FASTQ file from modified (chemically treated) sample"
    modified_r2: "R2 FASTQ file from modified (chemically treated) sample"
    untreated_r1: "R1 FASTQ file from untreated control sample"
    untreated_r2: "R2 FASTQ file from untreated control sample"
    primers_fa: "Optional FASTA file containing primer sequences for trimming"
    min_depth: "Minimum read depth required for calculating reactivity (default: 5000)"
    is_amplicon: "Set to true if data is from amplicon sequencing (default: false)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    String sample_name
    File target_fa
    File modified_r1
    File modified_r2
    File untreated_r1
    File untreated_r2
    File? primers_fa
    Int min_depth = 5000
    Boolean is_amplicon = false
    Int cpu_cores = 2
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Create output directory
    mkdir ~{sample_name}_output

    # Build ShapeMapper command
    cmd="shapemapper \
        --name ~{sample_name} \
        --target ~{target_fa} \
        --modified --R1 ~{modified_r1} --R2 ~{modified_r2} \
        --untreated --R1 ~{untreated_r1} --R2 ~{untreated_r2} \
        --min-depth ~{min_depth} \
        --out ~{sample_name}_output"

    # Add optional arguments
    ~{if is_amplicon then "cmd+=' --amplicon'" else ""}
    ~{if defined(primers_fa) then "cmd+=' --primers ${primers_fa}'" else ""}

    # Run ShapeMapper
    eval "$cmd"

    # Zip entire output directory for easier handling
    tar -czf ~{sample_name}_output.tar.gz ~{sample_name}_output
  >>>

  output {
    File output_tar = "~{sample_name}_output.tar.gz"
    File log_file = "~{sample_name}_shapemapper_log.txt"
  }

  runtime {
    docker: "getwilds/shapemapper:2.3"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
