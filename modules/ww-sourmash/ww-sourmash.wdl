# Modular WDL tasks for sourmash subcommands
# This file contains individual tasks that can be imported and used in workflows

version 1.0

task sketch {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Generate a sourmash sketch from BAM or FASTA file"
    outputs: {
      sig: "Sourmash sketch (.sig) file"
    }
  }

  parameter_meta {
    infile: "BAM (.bam) or FASTA (.fasta) file"
    bam_as_input: "Whether the input file is a BAM"
    k_value: "Value of k to use for sourmash sketch"
    scaled: "Scaled value for sourmash sketch"
    track_abundance: "Whether to track k-mer abundance"
    cpu_cores: "Number of CPU cores to use (only used for BAM input)"
    memory_gb: "Memory allocated for the task in GB"
    output_name: "Optional custom output name (defaults to input basename)"
  }

  input {
    File infile
    Boolean bam_as_input
    Int k_value = 31
    Int scaled = 1000
    Boolean track_abundance = true
    Int cpu_cores = 4
    Int memory_gb = 8
    String? output_name
  }

  String base_name = select_first([output_name, if bam_as_input then basename(infile, ".bam") else basename(infile, ".fasta")])
  String abund_param = if track_abundance then ",abund" else ""

  command <<<
    set -eo pipefail

    if ~{bam_as_input}; then
      echo "Running samtools fasta and sourmash sketch for: ~{base_name}.bam"
      samtools fasta -@ ~{cpu_cores} -0 /dev/null -s /dev/null "~{infile}" | \
      sourmash sketch dna -p k=~{k_value},scaled=~{scaled}~{abund_param} - -o "~{base_name}.sig"
    else
      echo "Running sourmash sketch for: ~{base_name}.fasta"
      sourmash sketch dna -p k=~{k_value},scaled=~{scaled}~{abund_param} "~{infile}" -o "~{base_name}.sig"
    fi
  >>>

  output {
    File sig = "~{base_name}.sig"
  }

  runtime {
    docker: "getwilds/sourmash:4.8.2"
    memory: "~{memory_gb} GB"
    disks: "local-disk 200 HDD"
  }
}

task gather {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Run sourmash gather with specific reference signature(s) included"
    outputs: {
      result: "CSV file of sourmash gather results"
    }
  }

  parameter_meta {
    query_sig: "Query sourmash sketch file"
    reference_databases_sigs: "Array of reference database (.zip) and/or signature (.sig) files to search against"
    memory_gb: "Memory allocated for the task in GB"
    output_name: "Optional custom output name (defaults to query basename)"
  }

  input {
    File query_sig
    Array[File] reference_databases_sigs
    Int memory_gb = 8
    String? output_name
  }

  String file_id = select_first([output_name, basename(query_sig, ".sig")])

  command <<<
    set -eo pipefail

    echo "Running sourmash gather on query: ~{file_id}.sig"
    sourmash gather "~{query_sig}" ~{sep=" " reference_databases_sigs} -o "~{file_id}.sourmash_gather.csv"
  >>>

  output {
    File result = "~{file_id}.sourmash_gather.csv"
  }

  runtime {
    docker: "getwilds/sourmash:4.8.2"
    memory: "~{memory_gb} GB"
    disks: "local-disk 200 HDD"
  }
}

task compare {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Generate similarity matrix from signature files using sourmash compare"
    outputs: {
      npy: "Numpy binary matrix file of angular similarity matrix",
      csv: "CSV file of angular similarity matrix"
    }
  }

  parameter_meta {
    sig_inputs: "Array of input signature files"
    save_name: "Name to use for output files"
    k_value: "Value of k used for sourmash sketch"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    Array[File] sig_inputs
    String save_name
    Int k_value
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Create a directory and move all signature files into it
    mkdir -p signatures
    for sig_file in ~{sep=" " sig_inputs}; do
      mv "$sig_file" signatures/
    done

    echo "Running sourmash compare"
    sourmash compare \
      --ksize ~{k_value} \
      --output "~{save_name}.npy" \
      --csv "~{save_name}.csv" \
      signatures/*.sig
  >>>

  output {
    File npy = "~{save_name}.npy"
    File csv = "~{save_name}.csv"
  }

  runtime {
    docker: "getwilds/sourmash:4.8.2"
    memory: "~{memory_gb} GB"
    disks: "local-disk 200 HDD"
  }
}
