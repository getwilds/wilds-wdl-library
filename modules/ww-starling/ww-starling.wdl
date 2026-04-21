## WILDS WDL module for STARLING ensemble generation.
## STARLING predicts coarse-grained structural ensembles of intrinsically
## disordered protein regions using a latent-space probabilistic denoising
## diffusion model with a Vision Transformer architecture.

version 1.0

#### TASK DEFINITIONS ####

task generate_ensemble {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Generate a structural ensemble for an intrinsically disordered protein sequence using STARLING"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-starling/ww-starling.wdl"
    outputs: {
        starling_file: "STARLING ensemble file containing predicted conformations",
        pdb_file: "PDB topology file converted from the STARLING ensemble",
        xtc_file: "XTC trajectory file converted from the STARLING ensemble"
    }
    topic: "proteomics,protein_disordered_structure"
    species: "any"
    operation: "protein_structure_prediction"
    in_sample_req: "none"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "starling_file:structure_prediction:binary_format,pdb_file:topology_data:pdb,xtc_file:trajectory_data:xtc"
    out_ref: "none"
  }

  parameter_meta {
    sequence: "Amino acid sequence string for the disordered protein region"
    sample_name: "Name identifier for the output files"
    num_conformations: "Number of conformations to generate in the ensemble"
    gpu_enabled: "Enable GPU for STARLING inference (sets device to cuda and requests GPU in runtime)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    String sequence
    String sample_name
    Int num_conformations = 400
    Boolean gpu_enabled = true
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  String device = if gpu_enabled then "cuda" else "cpu"

  command <<<
    set -eo pipefail

    # Generate the structural ensemble
    starling "~{sequence}" \
      -c ~{num_conformations} \
      --outname "~{sample_name}" \
      -d ~{device} \
      -r

    # Convert STARLING output to PDB and XTC formats
    starling2pdb "~{sample_name}.starling"
    starling2xtc "~{sample_name}.starling"
  >>>

  output {
    File starling_file = "~{sample_name}.starling"
    File pdb_file = "~{sample_name}_STARLING.pdb"
    File xtc_file = "~{sample_name}_STARLING.xtc"
  }

  runtime {
    docker: "getwilds/starling:2.0.0a3"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
    gpus: if gpu_enabled then "1" else "0"
  }
}

task generate_ensemble_batch {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Generate structural ensembles for multiple protein sequences from a FASTA file using STARLING"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-starling/ww-starling.wdl"
    outputs: {
        starling_files: "Array of STARLING ensemble files, one per sequence in the input FASTA",
        pdb_files: "Array of PDB topology files converted from the STARLING ensembles",
        xtc_files: "Array of XTC trajectory files converted from the STARLING ensembles"
    }
    topic: "proteomics,protein_disordered_structure"
    species: "any"
    operation: "protein_structure_prediction"
    in_sample_req: "fasta_file:protein_sequence:fasta"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "starling_files:structure_prediction:binary_format,pdb_files:topology_data:pdb,xtc_files:trajectory_data:xtc"
    out_ref: "none"
  }

  parameter_meta {
    fasta_file: "FASTA file containing one or more protein sequences"
    num_conformations: "Number of conformations to generate per sequence in the ensemble"
    gpu_enabled: "Enable GPU for STARLING inference (sets device to cuda and requests GPU in runtime)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File fasta_file
    Int num_conformations = 400
    Boolean gpu_enabled = true
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  String device = if gpu_enabled then "cuda" else "cpu"

  command <<<
    set -eo pipefail

    # Generate structural ensembles for all sequences in the FASTA
    starling "~{fasta_file}" \
      -c ~{num_conformations} \
      -d ~{device} \
      -r

    # Convert all STARLING outputs to PDB and XTC formats
    for f in *.starling; do
      starling2pdb "$f"
      starling2xtc "$f"
    done
  >>>

  output {
    Array[File] starling_files = glob("*.starling")
    Array[File] pdb_files = glob("*_STARLING.pdb")
    Array[File] xtc_files = glob("*_STARLING.xtc")
  }

  runtime {
    docker: "getwilds/starling:2.0.0a3"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
    gpus: if gpu_enabled then "1" else "0"
  }
}

task split_fasta {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Split a multi-sequence FASTA file into smaller batch files for parallel processing"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-starling/ww-starling.wdl"
    outputs: {
        batch_files: "Array of FASTA files, each containing up to sequences_per_batch sequences"
    }
    topic: "proteomics,protein_disordered_structure"
    species: "any"
    operation: "splitting"
    in_sample_req: "fasta_file:protein_sequence:fasta"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "batch_files:protein_sequence:fasta"
    out_ref: "none"
  }

  parameter_meta {
    fasta_file: "Input FASTA file containing multiple protein sequences"
    sequences_per_batch: "Maximum number of sequences per output batch file"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File fasta_file
    Int sequences_per_batch = 10
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -eo pipefail

    # Split FASTA into batches of N sequences each
    awk -v n=~{sequences_per_batch} '
      /^>/ {
        seq_count++
        if (seq_count == 1 || (seq_count - 1) % n == 0) {
          batch_num = int((seq_count - 1) / n) + 1
          outfile = sprintf("batch_%04d.fasta", batch_num)
        }
      }
      { print >> outfile }
    ' "~{fasta_file}"
  >>>

  output {
    Array[File] batch_files = glob("batch_*.fasta")
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task ensemble_info {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Query metadata and summary information from a STARLING ensemble file"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-starling/ww-starling.wdl"
    outputs: {
        info_file: "Text file containing ensemble metadata and summary statistics"
    }
    topic: "proteomics,protein_disordered_structure"
    species: "any"
    operation: "statistical_calculation"
    in_sample_req: "starling_file:structure_prediction:binary_format"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "info_file:quality_control_report:textual_format"
    out_ref: "none"
  }

  parameter_meta {
    starling_file: "STARLING ensemble file to query"
    sample_name: "Name identifier for the output file"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File starling_file
    String sample_name
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    starling2info "~{starling_file}" > "~{sample_name}_info.txt"
  >>>

  output {
    File info_file = "~{sample_name}_info.txt"
  }

  runtime {
    docker: "getwilds/starling:2.0.0a3"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
