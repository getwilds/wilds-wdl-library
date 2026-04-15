## WILDS WDL Module: ww-salmon
## Description: Module for RNA-seq quantification using Salmon
## Author: Taylor Firman
## Contact: tfirman@fredhutch.org

version 1.0

task build_index {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Build Salmon index from reference transcriptome FASTA file"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-salmon/ww-salmon.wdl"
    outputs: {
        salmon_index: "Compressed tarball containing the Salmon index for quantification"
    }
    topic: "transcriptomics,mapping,gene_expression"
    species: "any"
    operation: "indexing"
    in_sample_data: "none"
    in_sample_format: "none"
    in_ref_data: "rna_sequence"
    in_ref_format: "fasta"
    out_sample_data: "none"
    out_sample_format: "none"
    out_ref_data: "data_index"
    out_ref_format: "tar_format"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Quantify transcript expression from RNA-seq reads using Salmon. Supports both paired-end and single-end data."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-salmon/ww-salmon.wdl"
    outputs: {
        salmon_quant_dir: "Compressed tarball containing Salmon quantification results including abundance estimates",
        output_sample_name: "Sample name used for output files"
    }
    topic: "transcriptomics,gene_expression"
    species: "any"
    operation: "rna_seq_quantification"
    in_sample_data: "rna_sequence"
    in_sample_format: "fastq"
    in_ref_data: "data_index"
    in_ref_format: "tar_format"
    out_sample_data: "gene_expression_matrix"
    out_sample_format: "tar_format"
    out_ref_data: "none"
    out_ref_format: "none"
  }

  parameter_meta {
    salmon_index_dir: "Compressed tarball containing Salmon genome index"
    sample_name: "Name identifier for the sample"
    fastq_r1: "FASTQ file for read 1"
    fastq_r2: "Optional FASTQ file for read 2 (omit for single-end data)"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File salmon_index_dir
    String sample_name
    File fastq_r1
    File? fastq_r2
    Int cpu_cores = 8
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail

    # Extract the Salmon index
    mkdir -p salmon_index
    tar -xzf ~{salmon_index_dir} -C ./

    # Build read arguments based on paired-end or single-end data
    R2_FILE="~{if defined(fastq_r2) then select_first([fastq_r2]) else ""}"
    if [ -n "$R2_FILE" ]; then
      READ_ARGS="-1 ~{fastq_r1} -2 $R2_FILE"
    else
      READ_ARGS="-r ~{fastq_r1}"
    fi

    # Quantification with best practice parameters
    salmon quant \
        -i salmon_index \
        --libType A \
        $READ_ARGS \
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Merge Salmon quantification results from multiple samples into count and TPM matrices"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-salmon/ww-salmon.wdl"
    outputs: {
        tpm_matrix: "Tab-separated matrix of TPM (Transcripts Per Million) values for all samples",
        counts_matrix: "Tab-separated matrix of estimated read counts for all samples",
        sample_list: "Text file listing all sample names included in the matrices"
    }
    topic: "transcriptomics,gene_expression"
    species: "any"
    operation: "aggregation"
    in_sample_data: "gene_expression_matrix"
    in_sample_format: "tar_format"
    in_ref_data: "none"
    in_ref_format: "none"
    out_sample_data: "gene_expression_matrix,sample_id"
    out_sample_format: "tsv,textual_format"
    out_ref_data: "none"
    out_ref_format: "none"
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
