## WILDS WDL Module: ww-salmon
## Description: Module for RNA-seq quantification using Salmon
## Author: WILDS Team
## Contact: wilds@fredhutch.org

version 1.0

task build_index {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Build Salmon index from reference transcriptome FASTA file"
    outputs: {
        salmon_index: "Compressed tarball containing the Salmon index for quantification"
    }
  }

  parameter_meta {
    transcriptome_fasta: "Reference transcriptome in FASTA format"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File transcriptome_fasta
    Int cpu_cores = 8
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    # Build Salmon index with default parameters and decoys
    # Wait for salmon to fully complete before proceeding
    salmon index \
        -t ~{transcriptome_fasta} \
        --tmpdir=./tmp \
        -i salmon_index \
        -p ~{cpu_cores} \
        --gencode && echo "Salmon index build completed"

    # Ensure all filesystem operations are complete
    # Give extra time for any background writes to finish
    echo "Waiting for filesystem to stabilize..."
    sync
    sleep 5
    sync

    # Create tar archive of the index
    tar -czf salmon_index.tar.gz salmon_index
  >>>

  output {
    File salmon_index = "salmon_index.tar.gz"
  }

  runtime {
    docker: "getwilds/salmon:1.10.3"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task quantify {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Quantify transcript expression from paired-end RNA-seq reads using Salmon"
    outputs: {
        salmon_quant_dir: "Compressed tarball containing Salmon quantification results including abundance estimates",
        output_sample_name: "Sample name used for output files"
    }
  }

  parameter_meta {
    salmon_index_dir: "Compressed tarball containing Salmon genome index"
    sample_name: "Name identifier for the sample"
    fastq_r1: "FASTQ file for read 1"
    fastq_r2: "FASTQ file for read 2"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File salmon_index_dir
    String sample_name
    File fastq_r1
    File fastq_r2
    Int cpu_cores = 8
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    # Extract the Salmon index
    mkdir -p salmon_index
    tar -xzf ~{salmon_index_dir} -C ./

    # Paired-end quantification with best practice parameters
    salmon quant \
        -i salmon_index \
        --libType A \
        -1 ~{fastq_r1} \
        -2 ~{fastq_r2} \
        -o ~{sample_name}_quant \
        -p ~{cpu_cores} \
        --validateMappings \
        --gcBias \
        --seqBias

    # Tar up the output directory
    tar -czf ~{sample_name}_quant.tar.gz ~{sample_name}_quant
  >>>

  output {
    File salmon_quant_dir = "~{sample_name}_quant.tar.gz"
    String output_sample_name = "~{sample_name}"
  }

  runtime {
    docker: "getwilds/salmon:1.10.3"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task merge_results {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Merge Salmon quantification results from multiple samples into count and TPM matrices"
    outputs: {
        tpm_matrix: "Tab-separated matrix of TPM (Transcripts Per Million) values for all samples",
        counts_matrix: "Tab-separated matrix of estimated read counts for all samples",
        sample_list: "Text file listing all sample names included in the matrices"
    }
  }

  parameter_meta {
    salmon_quant_dirs: "Array of compressed Salmon quantification directories from multiple samples"
    sample_names: "Array of sample names corresponding to the quantification directories"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    Array[File] salmon_quant_dirs
    Array[String] sample_names
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Extract all quantification results
    mkdir -p quant_dirs
    for quant_dir in ~{sep=" " salmon_quant_dirs}; do
        base_name=$(basename $quant_dir .tar.gz)
        mkdir -p "quant_dirs/$base_name"
        tar -xzf $quant_dir -C quant_dirs/
    done

    # Create a list of sample names for reference
    echo "~{sep="\n" sample_names}" > sample_names.txt

    # Create a list of quant directories for quantmerge
    quant_dirs_list=()
    for sample in ~{sep=" " sample_names}; do
        quant_dirs_list+=("quant_dirs/${sample}_quant")
    done

    # Merge transcript-level TPM values
    salmon quantmerge \
        --quants "${quant_dirs_list[@]}" \
        --column TPM \
        -o tpm_matrix.tsv

    # Merge transcript-level estimated count values
    salmon quantmerge \
        --quants "${quant_dirs_list[@]}" \
        --column NumReads \
        -o counts_matrix.tsv
  >>>

  output {
    File tpm_matrix = "tpm_matrix.tsv"
    File counts_matrix = "counts_matrix.tsv"
    File sample_list = "sample_names.txt"
  }

  runtime {
    docker: "getwilds/salmon:1.10.3"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}
