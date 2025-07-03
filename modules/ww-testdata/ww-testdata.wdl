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
        ref_bed: "BED file covering the entire chromosome",
        ichor_gc_wig: "GC content WIG file for hg38",
        ichor_map_wig: "Mapping quality WIG file for hg38",
        ichor_centromeres: "Centromere locations for hg38",
        ichor_panel_of_norm_rds: "Panel of normals RDS file for hg38",
        annotsv_test_vcf: "Test VCF file for AnnotSV"
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

  call download_ichor_data { }

  call download_annotsv_vcf { }

  output {
    File ref_fasta = download_ref_data.fasta
    File ref_fasta_index = download_ref_data.fasta_index
    File ref_gtf = download_ref_data.gtf
    File ref_bed = download_ref_data.bed
    File ichor_gc_wig = download_ichor_data.wig_gc
    File ichor_map_wig = download_ichor_data.wig_map
    File ichor_centromeres = download_ichor_data.centromeres
    File ichor_panel_of_norm_rds = download_ichor_data.panel_of_norm_rds
    File annotsv_test_vcf = download_annotsv_vcf.test_vcf
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

task download_ichor_data {
  meta {
    description: "Downloads reference data for ichorCNA analysis on hg38"
    outputs: {
        wig_gc: "GC content WIG file for hg38",
        wig_map: "Mapping quality WIG file for hg38",
        centromeres: "Centromere locations for hg38",
        panel_of_norm_rds: "Panel of normals RDS file for hg38"
    }
  }

  parameter_meta {
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -euo pipefail

    # Download ichorCNA reference data files
    # These files are used for normalization and centromere locations in ichorCNA analysis
    wget -q --no-check-certificate -O gc_hg38_500kb.wig https://raw.githubusercontent.com/GavinHaLab/ichorCNA/b2bbce0a9997f31733f0f0ea4278cfba937ded41/inst/extdata/gc_hg38_500kb.wig
    wget -q --no-check-certificate -O map_hg38_500kb.wig https://raw.githubusercontent.com/GavinHaLab/ichorCNA/b2bbce0a9997f31733f0f0ea4278cfba937ded41/inst/extdata/map_hg38_500kb.wig
    wget -q --no-check-certificate -O GRCh38.GCA_000001405.2_centromere_acen.txt https://raw.githubusercontent.com/GavinHaLab/ichorCNA/b2bbce0a9997f31733f0f0ea4278cfba937ded41/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt
    wget -q --no-check-certificate -O nextera_hg38_500kb_median_normAutosome_median.rds_median.n9.gr.rds https://github.com/genome/docker-basespace_chromoseq/raw/refs/heads/master/workflow_files/nextera_hg38_500kb_median_normAutosome_median.rds_median.n9.gr.rds
  >>>

  output {
    File wig_gc = "gc_hg38_500kb.wig"
    File wig_map = "map_hg38_500kb.wig"
    File centromeres = "GRCh38.GCA_000001405.2_centromere_acen.txt"
    File panel_of_norm_rds = "nextera_hg38_500kb_median_normAutosome_median.rds_median.n9.gr.rds"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_annotsv_vcf {
  meta {
    description: "Downloads reference data for ichorCNA analysis on hg38"
    outputs: {
        test_vcf: "Test VCF file for AnnotSV"
    }
  }

  parameter_meta {
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -euo pipefail

    # Download AnnotSV test VCF file
    wget -q -O annotsv_test.vcf https://raw.githubusercontent.com/lgmgeo/AnnotSV/refs/heads/master/share/doc/AnnotSV/Example/test.vcf
  >>>

  output {
    File test_vcf = "annotsv_test.vcf"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
