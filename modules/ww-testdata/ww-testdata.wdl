## WILDS WDL module for downloading reference data for testing purposes.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

workflow testdata_example {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow for downloading reference data for WILDS WDL tests"
    url: "https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-testdata"
    outputs: {
        ref_fasta: "Reference genome FASTA file",
        ref_fasta_index: "Index file for the reference FASTA",
        ref_gtf: "GTF file containing gene annotations for the specified chromosome",
        ref_bed: "BED file covering the entire chromosome",
        r1_fastq: "R1 fastq file downloaded for the sample in question",
        r2_fastq: "R2 fastq file downloaded for the sample in question",
        cram: "CRAM file downloaded for the sample in question",
        crai: "Index file for the CRAM file",
        bam: "BAM file downloaded for the sample in question",
        bai: "Index file for the BAM file",
        ichor_gc_wig: "GC content WIG file for hg38",
        ichor_map_wig: "Mapping quality WIG file for hg38",
        ichor_centromeres: "Centromere locations for hg38",
        ichor_panel_of_norm_rds: "Panel of normals RDS file for hg38",
        dbsnp_vcf: "dbSNP VCF for hg38",
        known_indels_vcf: "Known indels VCF for hg38",
        gnomad_vcf: "Gnomad VCF for hg38",
        annotsv_test_vcf: "Test VCF file for AnnotSV",
        pasilla_counts: "Array of individual count files for each sample from Pasilla dataset",
        pasilla_sample_names: "Array of sample names corresponding to the count files",
        pasilla_sample_conditions: "Array of sample conditions corresponding to the count files",
        pasilla_gene_info: "Gene annotation information including gene symbols and descriptions",
        validation_report: "Validation report summarizing all outputs"
    }
  }


  # Pull down reference genome and index files for chr1
  call download_ref_data { input:
      chromo = "chr1",
      version = "hg38"
  }

  call download_fastq_data { }

  call interleave_fastq { input:
    r1_fq = download_fastq_data.r1_fastq,
    r2_fq = download_fastq_data.r2_fastq
  }

  call download_cram_data { input:
    ref_fasta = download_ref_data.fasta
  }

  call download_bam_data { }

  call download_ichor_data { }

  call download_tritonnp_data { }

  call download_dbsnp_vcf { input:
    region = "NC_000001.11:1-10000000",
    filter_name = "chr1"
  }

  call download_known_indels_vcf { input:
    region = "chr1:1-10000000",
    filter_name = "chr1"
  }

  call download_gnomad_vcf { input:
    region = "chr1:1-10000000",
    filter_name = "chr1"
  }

  call download_annotsv_vcf { }

  call generate_pasilla_counts { }

  call validate_outputs { input:
    ref_fasta = download_ref_data.fasta,
    ref_fasta_index = download_ref_data.fasta_index,
    ref_gtf = download_ref_data.gtf,
    ref_bed = download_ref_data.bed,
    r1_fastq = download_fastq_data.r1_fastq,
    r2_fastq = download_fastq_data.r2_fastq,
    inter_fastq = interleave_fastq.inter_fastq,
    cram = download_cram_data.cram,
    crai = download_cram_data.crai,
    bam = download_bam_data.bam,
    bai = download_bam_data.bai,
    ichor_gc_wig = download_ichor_data.wig_gc,
    ichor_map_wig = download_ichor_data.wig_map,
    ichor_centromeres = download_ichor_data.centromeres,
    ichor_panel_of_norm_rds = download_ichor_data.panel_of_norm_rds,
    dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf,
    known_indels_vcf = download_known_indels_vcf.known_indels_vcf,
    gnomad_vcf = download_gnomad_vcf.gnomad_vcf,
    annotsv_test_vcf = download_annotsv_vcf.test_vcf,
    pasilla_counts = generate_pasilla_counts.individual_count_files,
    pasilla_gene_info = generate_pasilla_counts.gene_info
  }

  output {
    # Outputs from the reference data download
    File ref_fasta = download_ref_data.fasta
    File ref_fasta_index = download_ref_data.fasta_index
    File ref_gtf = download_ref_data.gtf
    File ref_bed = download_ref_data.bed
    # Outputs from the fastq, cram, and bam data downloads
    File r1_fastq = download_fastq_data.r1_fastq
    File r2_fastq = download_fastq_data.r2_fastq
    File cram = download_cram_data.cram
    File crai = download_cram_data.crai
    File bam = download_bam_data.bam
    File bai = download_bam_data.bai
    # Outputs from the ichorCNA data download
    File ichor_gc_wig = download_ichor_data.wig_gc
    File ichor_map_wig = download_ichor_data.wig_map
    File ichor_centromeres = download_ichor_data.centromeres
    File ichor_panel_of_norm_rds = download_ichor_data.panel_of_norm_rds
    # Outputs from VCF downloads
    File dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf
    File known_indels_vcf = download_known_indels_vcf.known_indels_vcf
    File gnomad_vcf = download_gnomad_vcf.gnomad_vcf
    File annotsv_test_vcf = download_annotsv_vcf.test_vcf
    # Outputs from Pasilla DESeq2 count generation
    Array[File] pasilla_counts = generate_pasilla_counts.individual_count_files
    Array[String] pasilla_sample_names = generate_pasilla_counts.sample_names
    Array[String] pasilla_sample_conditions = generate_pasilla_counts.sample_conditions
    File pasilla_gene_info = generate_pasilla_counts.gene_info
    # Validation report summarizing all outputs
    File validation_report = validate_outputs.report
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

task download_fastq_data {
  meta {
    description: "Downloads small example FASTQ files for WILDS WDL test runs"
    outputs: {
        r1_fastq: "R1 fastq file downloaded for the sample in question",
        r2_fastq: "R2 fastq file downloaded for the sample in question"
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
    aws s3 cp --no-sign-request s3://gatk-test-data/wgs_fastq/NA12878_20k/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq .
    aws s3 cp --no-sign-request s3://gatk-test-data/wgs_fastq/NA12878_20k/H06HDADXX130110.1.ATCACGAT.20k_reads_2.fastq .
  >>>

  output {
    File r1_fastq = "H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq"
    File r2_fastq = "H06HDADXX130110.1.ATCACGAT.20k_reads_2.fastq"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task interleave_fastq {
  meta {
    description: "Interleaves a set of R1 and R2 FASTQ files"
    outputs: {
        inter_fastq: "Interleaved FASTQ"
    }
  }

  parameter_meta {
    r1_fq: "Forward (R1) FASTQ file"
    r2_fq: "Reverse (R2) FASTQ file"
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    File r1_fq
    File r2_fq
    Int cpu_cores = 2
    Int memory_gb = 4
  }

  command <<<
    # Read in both files in groups of four lines each
    paste <(gunzip -c ~{r1_fq} | paste - - - -) <(gunzip -c ~{r2_fq} | paste - - - -) | \
    # Interleave lines from each file, modifying quality score characters and including "+" lines
    awk -v OFS="\n" -v FS="\t" '{gsub(/#/, "!", $4); gsub(/#/, "!", $8); print($1,$2,"+",$4,$5,$6,"+",$8)}' | \
    gzip > interleaved.fastq.gz
  >>>

  output {
    File inter_fastq = "interleaved.fastq.gz"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_cram_data {
  meta {
    description: "Downloads small example CRAM files for WILDS WDL test runs"
    outputs: {
        cram: "CRAM file downloaded for the sample in question",
        crai: "Index file for the CRAM file"
    }
  }

  parameter_meta {
    ref_fasta: "Reference genome FASTA file to use for CRAM conversion"
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    File ref_fasta
    Int cpu_cores = 2
    Int memory_gb = 4
  }

  command <<<
    set -euo pipefail

    # Pull down BAM files from GATK test data bucket
    samtools view -@ ~{cpu_cores} -h -b s3://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bam chr1 | \
    samtools view -@ ~{cpu_cores} -s 0.1 -b - > NA12878.bam
    samtools index -@ ~{cpu_cores} NA12878.bam

    # Only keep primary alignments from chr1 (no supplementary alignments)
    samtools view -@ ~{cpu_cores} -h -f 0x2 NA12878.bam chr1 | \
      awk '/^@/ || ($7 == "=" || $7 == "chr1")' | \
      sed 's/\tSA:Z:[^\t]*//' | \
      sed '/^@SQ/d' | \
      sed '1a@SQ\tSN:chr1\tLN:248956422' | \
      samtools view -@ ~{cpu_cores} -b > NA12878_chr1.bam

    # Index the new BAM file
    samtools index -@ ~{cpu_cores} NA12878_chr1.bam

    # Convert BAM to CRAM using the provided reference FASTA
    samtools view -@ ~{cpu_cores} -C -T "~{ref_fasta}" -o NA12878_chr1.cram NA12878_chr1.bam
    samtools index -@ ~{cpu_cores} NA12878_chr1.cram
  >>>

  output {
    File cram = "NA12878_chr1.cram"
    File crai = "NA12878_chr1.cram.crai"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_bam_data {
  meta {
    description: "Downloads small example BAM files for WILDS WDL test runs"
    outputs: {
        bam: "BAM file downloaded for the sample in question",
        bai: "Index file for the BAM file"
    }
  }

  parameter_meta {
    filename: "Filename to save the BAM file as"
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    String filename = "NA12878_chr1.bam"
    Int cpu_cores = 2
    Int memory_gb = 4
  }

  command <<<
    set -euo pipefail

    # Pull down BAM files from GATK test data bucket
    samtools view -@ ~{cpu_cores} -h -b s3://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bam chr1 | \
    samtools view -@ ~{cpu_cores} -s 0.1 -b - > NA12878.bam
    samtools index -@ ~{cpu_cores} NA12878.bam

    # Only keep primary alignments from chr1 (no supplementary alignments)
    samtools view -@ ~{cpu_cores} -h -f 0x2 NA12878.bam chr1 | \
      awk '/^@/ || ($7 == "=" || $7 == "chr1")' | \
      sed 's/\tSA:Z:[^\t]*//' | \
      sed '/^@SQ/d' | \
      sed '1a@SQ\tSN:chr1\tLN:248956422' | \
      samtools view -@ ~{cpu_cores} -b > "~{filename}"

    # Index the new BAM file
    samtools index -@ ~{cpu_cores} "~{filename}"

    # Clean up intermediate files
    rm NA12878.bam NA12878.bam.bai
  >>>

  output {
    File bam = "~{filename}"
    File bai = "~{filename}.bai"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
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
    wget -q --no-check-certificate -O gc_hg38_500kb.wig https://raw.githubusercontent.com/broadinstitute/ichorCNA/refs/heads/master/inst/extdata/gc_hg38_500kb.wig
    wget -q --no-check-certificate -O map_hg38_500kb.wig https://raw.githubusercontent.com/broadinstitute/ichorCNA/refs/heads/master/inst/extdata/map_hg38_500kb.wig
    wget -q --no-check-certificate -O GRCh38.GCA_000001405.2_centromere_acen.txt https://raw.githubusercontent.com/broadinstitute/ichorCNA/refs/heads/master/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt
    wget -q --no-check-certificate -O HD_ULP_PoN_500kb_median_normAutosome_mapScoreFiltered_median.rds https://github.com/broadinstitute/ichorCNA/raw/refs/heads/master/inst/extdata/HD_ULP_PoN_500kb_median_normAutosome_mapScoreFiltered_median.rds
  >>>

  output {
    File wig_gc = "gc_hg38_500kb.wig"
    File wig_map = "map_hg38_500kb.wig"
    File centromeres = "GRCh38.GCA_000001405.2_centromere_acen.txt"
    File panel_of_norm_rds = "HD_ULP_PoN_500kb_median_normAutosome_mapScoreFiltered_median.rds"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_tritonnp_data {
  meta {
    description: "Downloads test data for TritonNP analysis"
    outputs: {
        annotation: "BED annotation file",
        plot_list: "Genes to plot",
        bam: "WGS test file",
        bam_index: "WGS test file index"
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

    # Download TritonNP reference data files
    wget -q --no-check-certificate -O AR.bed https://github.com/caalo/TritonNP/raw/refs/heads/main/reference_data/AR.bed
    wget -q --no-check-certificate -O plot_genes.txt https://github.com/caalo/TritonNP/raw/refs/heads/main/reference_data/plot_genes.txt
    wget -q --no-check-certificate -O NA12878.bam https://github.com/caalo/TritonNP/raw/refs/heads/main/test_data/NA12878.bam
    wget -q --no-check-certificate -O NA12878.bai https://github.com/caalo/TritonNP/raw/refs/heads/main/test_data/NA12878.bai
  >>>

  output {
    File annotation = "AR.bed"
    File plot_list = "plot_genes.txt"
    File bam = "NA12878.bam"
    File bam_index = "NA12878.bai"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_dbsnp_vcf {
  meta {
    description: "Downloads dbSNP VCF files for GATK workflows"
    outputs: {
        dbsnp_vcf: "dbSNP VCF file (filtered down if region specified)"
    }
  }

  parameter_meta {
    region: "Chromosomal region to filter the dbSNP vcf down to, e.g. NC_000001.11:1-10000000"
    filter_name: "Filename tag to save the dbSNP vcf with"
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    String? region
    String filter_name = "hg38"
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    set -euo pipefail

    # Create a full mapping file
    cat > chr_mapping.txt << EOF
    NC_000001.11 chr1
    NC_000002.12 chr2
    NC_000003.12 chr3
    NC_000004.12 chr4
    NC_000005.10 chr5
    NC_000006.12 chr6
    NC_000007.14 chr7
    NC_000008.11 chr8
    NC_000009.12 chr9
    NC_000010.11 chr10
    NC_000011.10 chr11
    NC_000012.12 chr12
    NC_000013.11 chr13
    NC_000014.9 chr14
    NC_000015.10 chr15
    NC_000016.10 chr16
    NC_000017.11 chr17
    NC_000018.10 chr18
    NC_000019.10 chr19
    NC_000020.11 chr20
    NC_000021.9 chr21
    NC_000022.11 chr22
    NC_000023.11 chrX
    NC_000024.10 chrY
    NC_012920.1 chrMT
    EOF

    # Download filtered dbSNP vcf from NCBI and rename chromosomes
    bcftools view ~{if defined(region) then "-r " + region else ""} \
      https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz | \
    bcftools annotate --rename-chrs chr_mapping.txt \
      -O z -o "dbsnp.~{filter_name}.vcf.gz"
  >>>

  output {
    File dbsnp_vcf = "dbsnp.~{filter_name}.vcf.gz"
  }

  runtime {
    docker: "getwilds/bcftools:1.19"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_known_indels_vcf {
  meta {
    description: "Downloads known indel VCF files for GATK workflows"
    outputs: {
        known_indels_vcf: "Known indels VCF file (filtered down if region specified)"
    }
  }

  parameter_meta {
    region: "Chromosomal region to filter the known indels vcf down to, e.g. chr1:1-10000000"
    filter_name: "Filename tag to save the known indels vcf with"
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    String? region
    String filter_name = "hg38"
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    # Download filtered known indels vcf from GATK
    bcftools view ~{if defined(region) then "-r " + region else ""} \
    https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O z -o "mills_1000g_known_indels.~{filter_name}.vcf.gz"
  >>>

  output {
    File known_indels_vcf = "mills_1000g_known_indels.~{filter_name}.vcf.gz"
  }

  runtime {
    docker: "getwilds/bcftools:1.19"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_gnomad_vcf {
  meta {
    description: "Downloads gnomad VCF files for GATK workflows"
    outputs: {
        gnomad_vcf: "Gnomad VCF file (filtered down if region specified)"
    }
  }

  parameter_meta {
    region: "Chromosomal region to filter the gnomad vcf down to, e.g. chr1:1-10000000"
    filter_name: "Filename tag to save the gnomad vcf with"
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    String? region
    String filter_name = "hg38"
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  command <<<
    # Download filtered gnomad vcf from GATK
    bcftools view ~{if defined(region) then "-r " + region else ""} \
    https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz \
    -O z -o "gnomad_af_only.~{filter_name}.vcf.gz"
  >>>

  output {
    File gnomad_vcf = "gnomad_af_only.~{filter_name}.vcf.gz"
  }

  runtime {
    docker: "getwilds/bcftools:1.19"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_annotsv_vcf {
  meta {
    description: "Downloads test VCF files for structural variant annotation workflows"
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
    wget -q --no-check-certificate -O annotsv_test.vcf https://raw.githubusercontent.com/lgmgeo/AnnotSV/refs/heads/master/share/doc/AnnotSV/Example/test2.vcf
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

task generate_pasilla_counts {
  meta {
    description: "Generate DESeq2 test count matrices and metadata using the pasilla Bioconductor dataset raw files"
    outputs: {
        individual_count_files: "Array of individual count files for each sample from Pasilla dataset",
        sample_names: "Array of sample names corresponding to the count files",
        sample_conditions: "Array of sample conditions corresponding to the count files",
        gene_info: "Gene annotation information including gene symbols and descriptions"
    }
  }

  parameter_meta {
    n_samples: "Number of samples to include (default: 7, max: 7 for pasilla dataset)"
    n_genes: "Approximate number of genes to include (will select top expressed genes)"
    condition_name: "Name for the condition column in metadata"
    output_prefix: "Prefix for output files"
    memory_gb: "Memory allocated for the task in GB"
    cpu_cores: "Number of CPU cores allocated for the task"
  }

  input {
    Int n_samples = 7
    Int n_genes = 10000
    String condition_name = "condition"
    String output_prefix = "pasilla"
    Int memory_gb = 4
    Int cpu_cores = 1
  }

  command <<<
    set -eo pipefail

    generate_pasilla_counts.R \
      --nsamples ~{n_samples} \
      --ngenes ~{n_genes} \
      --condition "~{condition_name}" \
      --prefix "~{output_prefix}"

    echo "Pasilla test data generation completed successfully"
  >>>

  output {
    Array[File] individual_count_files = glob("*.ReadsPerGene.out.tab")
    Array[String] sample_names = read_lines("~{output_prefix}_sample_names.txt")
    Array[String] sample_conditions = read_lines("~{output_prefix}_sample_conditions.txt") 
    File gene_info = "~{output_prefix}_gene_info.txt"
  }

  runtime {
    docker: "getwilds/deseq2:1.40.2"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task validate_outputs {
  meta {
    description: "Validates downloaded test data files to ensure they exist and are non-empty"
    outputs: {
        report: "Validation summary reporting file checks and basic statistics"
    }
  }

  parameter_meta {
    ref_fasta: "Reference genome FASTA file to validate"
    ref_fasta_index: "Reference FASTA index file to validate"
    ref_gtf: "GTF annotation file to validate"
    ref_bed: "BED file to validate"
    r1_fastq: "R1 FASTQ file to validate"
    r2_fastq: "R2 FASTQ file to validate"
    inter_fastq: "Interleaved FASTQ to validate"
    cram: "CRAM file to validate"
    crai: "CRAM index file to validate"
    bam: "BAM file to validate"
    bai: "BAM index file to validate"
    ichor_gc_wig: "ichorCNA GC content file to validate"
    ichor_map_wig: "ichorCNA mapping quality file to validate"
    ichor_centromeres: "ichorCNA centromere locations file to validate"
    ichor_panel_of_norm_rds: "ichorCNA panel of normals file to validate"
    dbsnp_vcf: "dbSNP VCF to validate"
    known_indels_vcf: "Known indels VCF to validate"
    gnomad_vcf: "gnomad VCF to validate"
    annotsv_test_vcf: "AnnotSV test VCF file to validate"
    pasilla_counts: "Array of individual count files for each sample from Pasilla dataset to validate"
    pasilla_gene_info: "Pasilla gene annotation information to validate"
    cpu_cores: "Number of CPU cores to use for validation"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    File ref_fasta
    File ref_fasta_index
    File ref_gtf
    File ref_bed
    File r1_fastq
    File r2_fastq
    File inter_fastq
    File cram
    File crai
    File bam
    File bai
    File ichor_gc_wig
    File ichor_map_wig
    File ichor_centromeres
    File ichor_panel_of_norm_rds
    File dbsnp_vcf
    File known_indels_vcf
    File gnomad_vcf
    File annotsv_test_vcf
    Array[File] pasilla_counts
    File pasilla_gene_info
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -euo pipefail

    echo "=== WILDS Test Data Validation Report ===" > validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt

    validation_passed=true

    # Check each file exists and is not empty
    if [[ -f "~{ref_fasta}" && -s "~{ref_fasta}" ]]; then
      echo "Reference FASTA: ~{ref_fasta} - PASSED" >> validation_report.txt
    else
      echo "Reference FASTA: ~{ref_fasta} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{ref_fasta_index}" && -s "~{ref_fasta_index}" ]]; then
      echo "Reference FASTA index: ~{ref_fasta_index} - PASSED" >> validation_report.txt
    else
      echo "Reference FASTA index: ~{ref_fasta_index} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{ref_gtf}" && -s "~{ref_gtf}" ]]; then
      echo "GTF file: ~{ref_gtf} - PASSED" >> validation_report.txt
    else
      echo "GTF file: ~{ref_gtf} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{ref_bed}" && -s "~{ref_bed}" ]]; then
      echo "BED file: ~{ref_bed} - PASSED" >> validation_report.txt
    else
      echo "BED file: ~{ref_bed} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{r1_fastq}" && -s "~{r1_fastq}" ]]; then
      echo "R1 FASTQ: ~{r1_fastq} - PASSED" >> validation_report.txt
    else
      echo "R1 FASTQ: ~{r1_fastq} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{r2_fastq}" && -s "~{r2_fastq}" ]]; then
      echo "R2 FASTQ: ~{r2_fastq} - PASSED" >> validation_report.txt
    else
      echo "R2 FASTQ: ~{r2_fastq} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{inter_fastq}" && -s "~{inter_fastq}" ]]; then
      echo "Interleaved FASTQ: ~{inter_fastq} - PASSED" >> validation_report.txt
    else
      echo "Interleaved FASTQ: ~{inter_fastq} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{cram}" && -s "~{cram}" ]]; then
      echo "CRAM file: ~{cram} - PASSED" >> validation_report.txt
    else
      echo "CRAM file: ~{cram} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{crai}" && -s "~{crai}" ]]; then
      echo "CRAM index: ~{crai} - PASSED" >> validation_report.txt
    else
      echo "CRAM index: ~{crai} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{bam}" && -s "~{bam}" ]]; then
      echo "BAM file: ~{bam} - PASSED" >> validation_report.txt
    else
      echo "BAM file: ~{bam} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{bai}" && -s "~{bai}" ]]; then
      echo "BAM index: ~{bai} - PASSED" >> validation_report.txt
    else
      echo "BAM index: ~{bai} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{ichor_gc_wig}" && -s "~{ichor_gc_wig}" ]]; then
      echo "ichorCNA GC WIG: ~{ichor_gc_wig} - PASSED" >> validation_report.txt
    else
      echo "ichorCNA GC WIG: ~{ichor_gc_wig} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{ichor_map_wig}" && -s "~{ichor_map_wig}" ]]; then
      echo "ichorCNA MAP WIG: ~{ichor_map_wig} - PASSED" >> validation_report.txt
    else
      echo "ichorCNA MAP WIG: ~{ichor_map_wig} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{ichor_centromeres}" && -s "~{ichor_centromeres}" ]]; then
      echo "ichorCNA centromeres: ~{ichor_centromeres} - PASSED" >> validation_report.txt
    else
      echo "ichorCNA centromeres: ~{ichor_centromeres} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{ichor_panel_of_norm_rds}" && -s "~{ichor_panel_of_norm_rds}" ]]; then
      echo "ichorCNA panel of normals: ~{ichor_panel_of_norm_rds} - PASSED" >> validation_report.txt
    else
      echo "ichorCNA panel of normals: ~{ichor_panel_of_norm_rds} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{dbsnp_vcf}" && -s "~{dbsnp_vcf}" ]]; then
      echo "dbSNP VCF: ~{dbsnp_vcf} - PASSED" >> validation_report.txt
    else
      echo "dbSNP VCF: ~{dbsnp_vcf} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{known_indels_vcf}" && -s "~{known_indels_vcf}" ]]; then
      echo "Known Indels VCF: ~{known_indels_vcf} - PASSED" >> validation_report.txt
    else
      echo "Known Indels VCF: ~{known_indels_vcf} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{gnomad_vcf}" && -s "~{gnomad_vcf}" ]]; then
      echo "Gnomad VCF: ~{gnomad_vcf} - PASSED" >> validation_report.txt
    else
      echo "Gnomad VCF: ~{gnomad_vcf} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    if [[ -f "~{annotsv_test_vcf}" && -s "~{annotsv_test_vcf}" ]]; then
      echo "AnnotSV test VCF: ~{annotsv_test_vcf} - PASSED" >> validation_report.txt
    else
      echo "AnnotSV test VCF: ~{annotsv_test_vcf} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    # Check if all pasilla count files exist and are non-empty
    pasilla_count_passed=true
    for count_file in ~{sep=' ' pasilla_counts}; do
      if [[ -f "$count_file" && -s "$count_file" ]]; then
        echo "Pasilla count file: $count_file - PASSED" >> validation_report.txt
      else
        echo "Pasilla count file: $count_file - MISSING OR EMPTY" >> validation_report.txt
        pasilla_count_passed=false
      fi
    done
    if [[ "$pasilla_count_passed" == "false" ]]; then
      validation_passed=false
    fi

    if [[ -f "~{pasilla_gene_info}" && -s "~{pasilla_gene_info}" ]]; then
      echo "Pasilla gene info: ~{pasilla_gene_info} - PASSED" >> validation_report.txt
    else
      echo "Pasilla gene info: ~{pasilla_gene_info} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    {
      echo ""
      echo "=== Validation Summary ==="
      echo "Total files validated: 20"
    } >> validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> validation_report.txt
    else
      echo "Overall Status: FAILED" >> validation_report.txt
      exit 1
    fi

    cat validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
