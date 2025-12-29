## WILDS WDL Module for DIAMOND
## Sequence alignment using the DIAMOND tool.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task make_database {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Create a DIAMOND database from a FASTA file"
    outputs: {
        diamond_db: "DIAMOND database file (.dmnd) created from input FASTA"
    }
  }

  parameter_meta {
    fasta: "Input FASTA file containing protein sequences"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File fasta
    Int memory_gb = 2
    Int cpu_cores = 1
  }

  String fasta_base = sub(basename(fasta), ".*/", "")

  command <<<
    set -eo pipefail

    diamond makedb --in "~{fasta}" --db "~{fasta_base}.db.dmnd"
  >>>

  output {
    File diamond_db = "~{fasta_base}.db.dmnd"
  }

  runtime {
    docker: "diamond:2.1.16"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task diamond_blastp {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Align protein sequences to DIAMOND database using BLASTP"
    outputs: {
        aln: "Compressed alignment file in tabular format (outfmt 6)"
    }
  }

  parameter_meta {
    diamond_db: "DIAMOND database file (.dmnd)"
    query: "Query FASTA file containing protein sequences to align"
    align_id: "Minimum identity percentage for alignment (default: 50)"
    top_pct: "Top percentage of hits to keep (0 = keep all, default: 0)"
    query_cover: "Minimum query coverage percentage (default: 50)"
    subject_cover: "Minimum subject coverage percentage (default: 50)"
    blocksize: "Block size in billions of sequence letters to be processed at a time."
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File diamond_db
    File query
    String align_id = "50"
    String top_pct = "0"
    String query_cover = "50"
    String subject_cover = "50"
    Float blocksize = 2.0
    Int memory_gb = 2
    Int cpu_cores = 1
  }

  String db_name = basename(diamond_db, ".dmnd")
  String query_name = basename(query)

  command <<<
    set -eo pipefail

    diamond \
        blastp \
        --db "~{diamond_db}" \
        --query "~{query}" \
        --out "~{query_name}.~{db_name}.aln" \
        --outfmt 6 \
        --id "~{align_id}" \
        --top "~{top_pct}" \
        --query-cover "~{query_cover}" \
        --subject-cover "~{subject_cover}" \
        --block-size ~{blocksize} \
        --threads ~{cpu_cores}

    gzip "~{query_name}.~{db_name}.aln"
  >>>

  output {
    File aln = "~{query_name}.~{db_name}.aln.gz"
  }

  runtime {
    docker: "diamond:2.1.16"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}
