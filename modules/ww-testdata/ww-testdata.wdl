## WILDS WDL module for downloading reference data for testing purposes.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

workflow testdata_example {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "WDL workflow for downloading reference data for WILDS WDL tests"
    url: "https://github.com/getwilds/ww-testdata"
    outputs: {
        ref_fasta: "Reference genome FASTA file",
        ref_fasta_index: "Index file for the reference FASTA",
        ref_gtf: "GTF file containing gene annotations for the specified chromosome",
        ref_bed: "BED file covering the entire chromosome"
    }
  }

  parameter_meta {
    chromo: "Chromosome to download (e.g., chr1, chr2, etc.)"
    version: "Reference genome version (e.g., hg38, hg19)"
  }

  input {
    String chromo = "chr1"
    String version = "hg38"
  }

  # Pull down reference genome and index files for the specified chromosome
  call download_ref_data { input:
      chromo = chromo,
      version = version
  }

  output {
    File ref_fasta = download_ref_data.fasta
    File ref_fasta_index = download_ref_data.fasta_index
    File ref_gtf = download_ref_data.gtf
    File ref_bed = download_ref_data.bed
  }
}

task download_ref_data {
  meta {
    description: "Downloads reference genome and index files for WILDS WDL test runs"
    outputs: {
        fasta: "Reference genome FASTA file",
        fasta_index: "Index file for the reference FASTA",
        gtf: "GTF file containing gene annotations for the specified chromosome",
        bed: "BED file covering the entire chromosome"
    }
  }

  parameter_meta {
    chromo: "Chromosome to download (e.g., chr1, chr2, etc.)"
    version: "Reference genome version (e.g., hg38, hg19)"
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    String chromo = "chr1"
    String version = "hg38"
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -euo pipefail

    # Download chromosome fasta
    wget -q -O "~{chromo}.fa.gz" "http://hgdownload.soe.ucsc.edu/goldenPath/~{version}/chromosomes/~{chromo}.fa.gz"
    gunzip "~{chromo}.fa.gz"

    # Create FASTA index file (.fai) for bcftools and other tools
    samtools faidx "~{chromo}.fa"

    # Download chromosome 1 GTF file
    wget -q -O "~{version}.ncbiRefSeq.gtf.gz" "http://hgdownload.soe.ucsc.edu/goldenPath/~{version}/bigZips/genes/~{version}.ncbiRefSeq.gtf.gz"
    gunzip "~{version}.ncbiRefSeq.gtf.gz"
    # Extract only chromosome annotations
    grep "^~{chromo}[[:space:]]" "~{version}.ncbiRefSeq.gtf" > "~{chromo}.gtf"
    rm "~{version}.ncbiRefSeq.gtf"

    # Create a BED file covering the entire chromosome
    # Get chromosome length from the FASTA file
    CHR_LENGTH=$(($(grep -v "^>" "~{chromo}.fa" | tr -d '\n' | wc -c)))
    echo -e "~{chromo}\t0\t${CHR_LENGTH}" > "~{chromo}.bed"
  >>>

  output {
    File fasta = "~{chromo}.fa"
    File fasta_index = "~{chromo}.fa.fai"
    File gtf = "~{chromo}.gtf"
    File bed = "~{chromo}.bed"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
