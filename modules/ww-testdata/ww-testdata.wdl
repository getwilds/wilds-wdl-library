## WILDS WDL module for downloading reference data for testing purposes.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task download_ref_data {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads reference genome and index files for WILDS WDL test runs"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        fasta: "Reference genome FASTA file",
        fasta_index: "Index file for the reference FASTA",
        dict: "Dictionary file for the reference FASTA",
        gtf: "GTF file containing gene annotations for the specified chromosome",
        bed: "BED file covering the entire chromosome"
    }
  }

  parameter_meta {
    chromo: "Chromosome to download (e.g., chr1, chr2, etc.)"
    version: "Reference genome version (e.g., hg38, hg19)"
    region: "Optional region coordinates to extract from chromosome in format '1-30000000'. If not specified, uses entire chromosome."
    output_name: "Optional name for output files (default: uses chromo name)"
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    String chromo = "chr1"
    String version = "hg38"
    String? region
    String? output_name
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  String final_output_name = select_first([output_name, chromo])

  command <<<
    set -euo pipefail

    # Download chromosome fasta
    wget -q -O "~{chromo}.fa.gz" "http://hgdownload.soe.ucsc.edu/goldenPath/~{version}/chromosomes/~{chromo}.fa.gz"
    gunzip "~{chromo}.fa.gz"
    mv "~{chromo}.fa" temp.fa

    # Subset to specified region if provided
    REGION="~{if defined(region) then region else ""}"
    if [ -n "$REGION" ]; then
      samtools faidx temp.fa
      samtools faidx temp.fa "~{chromo}:$REGION" | sed "s/>~{chromo}:$REGION/>~{chromo}/" > "~{final_output_name}.fa"
      rm temp.fa temp.fa.fai
    else
      mv temp.fa "~{final_output_name}.fa"
    fi

    # Create FASTA index file (.fai) for bcftools and other tools
    samtools faidx "~{final_output_name}.fa"

    # Create FASTA dictionary file (.dict) for GATK and other tools
    samtools dict "~{final_output_name}.fa" > "~{final_output_name}.dict"

    # Download GTF file
    wget -q -O "~{version}.ncbiRefSeq.gtf.gz" "http://hgdownload.soe.ucsc.edu/goldenPath/~{version}/bigZips/genes/~{version}.ncbiRefSeq.gtf.gz"
    gunzip "~{version}.ncbiRefSeq.gtf.gz"
    # Extract only chromosome annotations
    grep "^~{chromo}[[:space:]]" "~{version}.ncbiRefSeq.gtf" > temp.gtf
    rm "~{version}.ncbiRefSeq.gtf"

    # If region is specified, filter GTF to only include genes in that region
    if [ -n "$REGION" ]; then
      # Extract start and end positions from region (format: 1-30000000)
      REGION_START=$(echo "$REGION" | cut -d'-' -f1)
      REGION_END=$(echo "$REGION" | cut -d'-' -f2)
      # Filter GTF to only include features within the region
      awk -v start="$REGION_START" -v end="$REGION_END" '$4 <= end && $5 >= start' temp.gtf > "~{final_output_name}.gtf"
      rm temp.gtf
    else
      mv temp.gtf "~{final_output_name}.gtf"
    fi

    # Create a BED file covering the entire chromosome/region
    # Get length from the FASTA file
    CHR_LENGTH=$(($(grep -v "^>" "~{final_output_name}.fa" | tr -d '\n' | wc -c)))
    echo -e "~{chromo}\t0\t${CHR_LENGTH}" > "~{final_output_name}.bed"
  >>>

  output {
    File fasta = "~{final_output_name}.fa"
    File fasta_index = "~{final_output_name}.fa.fai"
    File dict = "~{final_output_name}.dict"
    File gtf = "~{final_output_name}.gtf"
    File bed = "~{final_output_name}.bed"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_fastq_data {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads small example FASTQ files for WILDS WDL test runs. Renames to Illumina naming convention with optional gzip compression."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        r1_fastq: "R1 fastq file downloaded for the sample in question",
        r2_fastq: "R2 fastq file downloaded for the sample in question"
    }
  }

  parameter_meta {
    prefix: "Sample prefix for output filenames (default: 'testdata')"
    gzip_output: "Compress output files with gzip (default: false)"
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    Boolean gzip_output = false
    String prefix = "testdata"
    Int cpu_cores = 1
    Int memory_gb = 4
  }

  # Determine output filenames based on prefix and gzip setting
  String gz_ext = if gzip_output then ".gz" else ""
  String r1_base = "~{prefix}_S1_L001_R1_001.fastq"
  String r2_base = "~{prefix}_S1_L001_R2_001.fastq"
  String r1_output = "~{r1_base}~{gz_ext}"
  String r2_output = "~{r2_base}~{gz_ext}"

  command <<<
    set -eo pipefail

    # Download example FASTQ files from GATK test data bucket
    aws s3 cp --no-sign-request s3://gatk-test-data/wgs_fastq/NA12878_20k/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq "~{r1_base}"
    aws s3 cp --no-sign-request s3://gatk-test-data/wgs_fastq/NA12878_20k/H06HDADXX130110.1.ATCACGAT.20k_reads_2.fastq "~{r2_base}"

    # Optionally gzip the output FASTQ files
    if [ "~{gzip_output}" == "true" ]; then
      gzip "~{r1_base}"
      gzip "~{r2_base}"
    fi
  >>>

  output {
    File r1_fastq = r1_output
    File r2_fastq = r2_output
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task interleave_fastq {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Interleaves a set of R1 and R2 FASTQ files"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads small example CRAM files for WILDS WDL test runs"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
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

    # Download the small BAM file from GATK test data bucket using AWS CLI
    # This approach works in all environments (Docker, Apptainer, local) unlike samtools S3 streaming
    # Suboptimal approach, but necessary to ensure functionality across environments
    aws s3 cp --no-sign-request \
      s3://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bam \
      NA12878_full.bam
    aws s3 cp --no-sign-request \
      s3://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bai \
      NA12878_full.bam.bai

    # Extract chr1 and subsample to reduce size (using lower subsample rate for CRAM)
    samtools view -@ ~{cpu_cores} -h -b NA12878_full.bam chr1 | \
    samtools view -@ ~{cpu_cores} -s 0.05 -b - > NA12878.bam
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

    # Clean up intermediate files
    rm NA12878_full.bam NA12878_full.bam.bai NA12878.bam NA12878.bam.bai NA12878_chr1.bam NA12878_chr1.bam.bai
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads small example BAM files for WILDS WDL test runs"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
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

    # Download the small BAM file from GATK test data bucket using AWS CLI
    # This approach works in all environments (Docker, Apptainer, local) unlike samtools S3 streaming
    # Suboptimal approach, but necessary to ensure functionality across environments
    aws s3 cp --no-sign-request \
      s3://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bam \
      NA12878_full.bam
    aws s3 cp --no-sign-request \
      s3://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bai \
      NA12878_full.bam.bai

    # Extract chr1 and subsample to reduce size
    samtools view -@ ~{cpu_cores} -h -b NA12878_full.bam chr1 | \
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
    rm NA12878_full.bam NA12878_full.bam.bai NA12878.bam NA12878.bam.bai
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads reference data for ichorCNA analysis on hg38"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads test data for TritonNP analysis"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        annotation: "BED annotation file",
        plot_list: "Genes to plot",
        bam: "WGS test file",
        bam_index: "WGS test file index",
        bias: "GC bias",
        reference: "hg19 reference genome fasta",
        reference_index: "Index for hg19 reference genome fasta"
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
    wget -q --no-check-certificate -O NA12878.GC_bias.txt https://github.com/caalo/TritonNP/raw/refs/heads/main/test_data/NA12878.GC_bias.txt
    #Download the entire hg19 reference genome
    wget hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
    tar -zxvf chromFa.tar.gz
    cat chr*.fa > hg19.fa
    samtools faidx hg19.fa

  >>>

  output {
    File annotation = "AR.bed"
    File plot_list = "plot_genes.txt"
    File bam = "NA12878.bam"
    File bam_index = "NA12878.bai"
    File bias = "NA12878.GC_bias.txt"
    File reference = "hg19.fa"
    File reference_index = "hg19.fa.fai"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_dbsnp_vcf {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads dbSNP VCF files for GATK workflows"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        dbsnp_vcf: "dbSNP VCF file (filtered down if region specified)",
        dbsnp_vcf_index: "Index file for the dbSNP VCF"
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

    # Index the filtered VCF
    bcftools index --tbi "dbsnp.~{filter_name}.vcf.gz"
  >>>

  output {
    File dbsnp_vcf = "dbsnp.~{filter_name}.vcf.gz"
    File dbsnp_vcf_index = "dbsnp.~{filter_name}.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/bcftools:1.19"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_known_indels_vcf {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads known indel VCF files for GATK workflows"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        known_indels_vcf: "Known indels VCF file (filtered down if region specified)",
        known_indels_vcf_index: "Index file for the known indels VCF"
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
    # Download filtered known indels vcf from 1000 Genomes EBI FTP
    # Note: The original Google Cloud Storage URL (genomics-public-data) now requires authentication
    bcftools view ~{if defined(region) then "-r " + region else ""} \
    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz \
    -O z -o "mills_1000g_known_indels.~{filter_name}.vcf.gz"

    # Index the filtered VCF
    bcftools index --tbi "mills_1000g_known_indels.~{filter_name}.vcf.gz"
  >>>

  output {
    File known_indels_vcf = "mills_1000g_known_indels.~{filter_name}.vcf.gz"
    File known_indels_vcf_index = "mills_1000g_known_indels.~{filter_name}.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/bcftools:1.19"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_gnomad_vcf {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads gnomad VCF files for GATK workflows"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        gnomad_vcf: "Gnomad VCF file (filtered down if region specified)",
        gnomad_vcf_index: "Index file for the gnomad VCF"
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

    # Index the filtered VCF
    bcftools index --tbi "gnomad_af_only.~{filter_name}.vcf.gz"
  >>>

  output {
    File gnomad_vcf = "gnomad_af_only.~{filter_name}.vcf.gz"
    File gnomad_vcf_index = "gnomad_af_only.~{filter_name}.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/bcftools:1.19"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_annotsv_vcf {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Downloads test VCF files for structural variant annotation workflows"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Generate DESeq2 test count matrices and metadata using the pasilla Bioconductor dataset raw files"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
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

task download_test_transcriptome {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Download a small test transcriptome for RNA-seq quantification testing. NOTE: This uses GENCODE (Ensembl) annotations."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        transcriptome_fasta: "Small test transcriptome FASTA file containing protein-coding transcripts"
    }
  }

  parameter_meta {
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -eo pipefail

    # Download protein-coding transcriptome from GENCODE for testing
    # This is a relatively small file (~46MB compressed) containing all human protein-coding transcripts
    curl -L -o test_transcriptome.fa.gz \
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.pc_transcripts.fa.gz"

    # Decompress the file
    gunzip test_transcriptome.fa.gz
  >>>

  output {
    File transcriptome_fasta = "test_transcriptome.fa"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task create_clean_amplicon_reference {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Extract and clean a reference sequence region for saturation mutagenesis analysis"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        clean_fasta: "Cleaned reference FASTA file with no ambiguous bases",
        clean_fasta_index: "Index file for the cleaned reference FASTA",
        clean_dict: "Dictionary file for the cleaned reference FASTA"
    }
  }

  parameter_meta {
    input_fasta: "Input reference FASTA file"
    region: "Region to extract in format 'chr:start-end' (e.g., 'chr1:1000-2000'). If not specified, uses entire sequence."
    output_name: "Name for the output reference (default: 'amplicon')"
    replace_n_with: "Base to replace N's with (default: 'A'). Use empty string to fail if N's are found."
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File input_fasta
    String? region
    String output_name = "amplicon"
    String replace_n_with = "A"
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -eo pipefail

    # Extract region if specified, otherwise use entire sequence
    if [ -n "~{region}" ]; then
      samtools faidx "~{input_fasta}"
      samtools faidx "~{input_fasta}" "~{region}" > temp_extract.fa

      # Replace the header with just the chromosome name
      sed "s/^>.*/>~{output_name}/" temp_extract.fa > temp.fa
      rm temp_extract.fa
    else
      cp "~{input_fasta}" temp.fa
    fi

    # Check for N bases and handle according to replace_n_with parameter
    n_count=$(grep -v "^>" temp.fa | grep -o "N" | wc -l || true)

    if [ "$n_count" -gt 0 ]; then
      echo "Found $n_count N bases in the reference sequence"

      if [ -z "~{replace_n_with}" ]; then
        echo "ERROR: N bases found and replace_n_with is empty. Cannot proceed."
        exit 1
      else
        echo "Replacing N bases with '~{replace_n_with}'"
        # Replace N's (both upper and lowercase) in the sequence lines only
        # Removes ambiguous bases (N's) that cause issues with GATK AnalyzeSaturationMutagenesis
        awk '/^>/ {print; next} {gsub(/[Nn]/, "~{replace_n_with}"); print}' temp.fa > "~{output_name}.fa"
      fi
    else
      echo "No N bases found in the reference sequence"
      mv temp.fa "~{output_name}.fa"
    fi

    # Verify no ambiguous bases remain
    remaining_n=$(grep -v "^>" "~{output_name}.fa" | grep -o "[^ACGTacgt]" | wc -l || true)
    if [ "$remaining_n" -gt 0 ]; then
      echo "ERROR: Non-ACGT bases still present after cleaning"
      exit 1
    fi

    echo "Reference cleaned successfully: $(basename '~{output_name}.fa')"

    # Create index and dictionary
    samtools faidx "~{output_name}.fa"
    samtools dict "~{output_name}.fa" > "~{output_name}.dict"
  >>>

  output {
    File clean_fasta = "~{output_name}.fa"
    File clean_fasta_index = "~{output_name}.fa.fai"
    File clean_dict = "~{output_name}.dict"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task create_gdc_manifest {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Create a test GDC manifest file with small open-access files for testing gdc-client downloads"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        manifest: "GDC manifest file containing test file UUIDs"
    }
  }

  command <<<
    set -eo pipefail

    # Create a test manifest file with small open-access TCGA files
    # Format: id, filename, md5, size, state (tab-separated)
    cat > gdc_test_manifest.txt <<'EOF'
id	filename	md5	size	state
6e811713-17b0-4413-a756-af178269824f	TARGET_AML_SampleMatrix_Validation_20180914.xlsx	c2070b78d418c134f48d9b8098c9f7ac	172979	released
4e89ba70-022e-48a3-a8f8-04f5720fb2d0	5df0dac1-d9d2-4e2d-b3dc-63279926a402.targeted_sequencing.aliquot_ensemble_raw.maf.gz	aa97da746509a18676adda0f086dedaa	9429	released
52fd584b-9ca3-4f7e-bdd5-fd9dce3d630b	TARGET_AML_CDE_20230524.xlsx	fc5b42f89b1ac84699b8b10dafed02dc	27667	released
EOF

    echo "Created test GDC manifest with 3 small open-access files"
  >>>

  output {
    File manifest = "gdc_test_manifest.txt"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    memory: "2 GB"
    cpu: 1
  }
}

task download_shapemapper_data {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Downloads ShapeMapper example data (TPP riboswitch) from the official repository for testing RNA structure analysis workflows"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        target_fa: "Target RNA FASTA file (TPP riboswitch sequence)",
        modified_r1: "R1 FASTQ file from modified/treated sample (TPPplus)",
        modified_r2: "R2 FASTQ file from modified/treated sample (TPPplus)",
        untreated_r1: "R1 FASTQ file from untreated control sample (TPPminus)",
        untreated_r2: "R2 FASTQ file from untreated control sample (TPPminus)"
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
    set -eo pipefail

    BASE_URL="https://raw.githubusercontent.com/Weeks-UNC/shapemapper2/master/example_data"

    # Download target RNA FASTA file
    echo "Downloading TPP.fa target sequence..."
    wget -q --no-check-certificate -O TPP.fa "${BASE_URL}/TPP.fa"

    # Download and concatenate TPPplus (modified/treated) R1 reads
    echo "Downloading TPPplus (modified) R1 reads..."
    for part in aa ab ac ad; do
      wget -q --no-check-certificate -O "TPPplus_R1_${part}.fastq.gz" \
        "${BASE_URL}/TPPplus/TPPplus_R1_${part}.fastq.gz"
    done
    cat TPPplus_R1_*.fastq.gz > TPPplus_R1.fastq.gz
    rm TPPplus_R1_aa.fastq.gz TPPplus_R1_ab.fastq.gz TPPplus_R1_ac.fastq.gz TPPplus_R1_ad.fastq.gz

    # Download and concatenate TPPplus (modified/treated) R2 reads
    echo "Downloading TPPplus (modified) R2 reads..."
    for part in aa ab ac ad; do
      wget -q --no-check-certificate -O "TPPplus_R2_${part}.fastq.gz" \
        "${BASE_URL}/TPPplus/TPPplus_R2_${part}.fastq.gz"
    done
    cat TPPplus_R2_*.fastq.gz > TPPplus_R2.fastq.gz
    rm TPPplus_R2_aa.fastq.gz TPPplus_R2_ab.fastq.gz TPPplus_R2_ac.fastq.gz TPPplus_R2_ad.fastq.gz

    # Download and concatenate TPPminus (untreated) R1 reads
    echo "Downloading TPPminus (untreated) R1 reads..."
    for part in aa ab ac ad; do
      wget -q --no-check-certificate -O "TPPminus_R1_${part}.fastq.gz" \
        "${BASE_URL}/TPPminus/TPPminus_R1_${part}.fastq.gz"
    done
    cat TPPminus_R1_*.fastq.gz > TPPminus_R1.fastq.gz
    rm TPPminus_R1_aa.fastq.gz TPPminus_R1_ab.fastq.gz TPPminus_R1_ac.fastq.gz TPPminus_R1_ad.fastq.gz

    # Download and concatenate TPPminus (untreated) R2 reads
    echo "Downloading TPPminus (untreated) R2 reads..."
    for part in aa ab ac ad; do
      wget -q --no-check-certificate -O "TPPminus_R2_${part}.fastq.gz" \
        "${BASE_URL}/TPPminus/TPPminus_R2_${part}.fastq.gz"
    done
    cat TPPminus_R2_*.fastq.gz > TPPminus_R2.fastq.gz
    rm TPPminus_R2_aa.fastq.gz TPPminus_R2_ab.fastq.gz TPPminus_R2_ac.fastq.gz TPPminus_R2_ad.fastq.gz

    echo "ShapeMapper example data download complete"
  >>>

  output {
    File target_fa = "TPP.fa"
    File modified_r1 = "TPPplus_R1.fastq.gz"
    File modified_r2 = "TPPplus_R2.fastq.gz"
    File untreated_r1 = "TPPminus_R1.fastq.gz"
    File untreated_r2 = "TPPminus_R2.fastq.gz"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_test_cellranger_ref {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Download a minimal Cell Ranger reference for testing"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        ref_tar: "Cell Ranger reference transcriptome tarball"
    }
  }

  parameter_meta {
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    Int cpu_cores = 2
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    # Download a minimal human reference from Swiss Bioinformatics Institute
    # Only chromosomes 21 and 22
    # https://sib-swiss.github.io/single-cell-training-archived/2023.3/day1/introduction_cellranger/#__tabbed_1_1
    # Emma manually extracted and inspected the files within for any obvious
    # safety issues

    echo "Downloading small test reference (728 MB)..."
    curl -O https://single-cell-transcriptomics.s3.eu-central-1.amazonaws.com/cellranger_index.tar.gz
    echo "Reference download complete"
  >>>

  output {
    File ref_tar = "cellranger_index.tar.gz"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task download_diamond_data {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Download E. coli proteins and create a subset as a test query"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        reference: "Full E. coli proteome FASTA file",
        query: "Subset of first 10 sequences"
    }
  }

  parameter_meta {
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -eo pipefail

    # Download E. coli Swiss-Prot proteins and save as ecoli_proteins.fasta
    curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000625/UP000000625_83333.fasta.gz | gunzip > ecoli_proteins.fasta

    # Create subset with first 10 sequences
    awk '/^>/ {n++} n<=10' ecoli_proteins.fasta > ecoli_subset.fasta
  >>>

  output {
    File reference = "ecoli_proteins.fasta"
    File query = "ecoli_subset.fasta"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_glimpse2_genetic_map {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads genetic map files for GLIMPSE2 imputation from the official GLIMPSE repository"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        genetic_map: "Genetic map file for the specified chromosome"
    }
  }

  parameter_meta {
    chromosome: "Chromosome to download genetic map for (e.g., chr22)"
    genome_build: "Genome build version (b37 or b38)"
    cpu_cores: "Number of CPU cores to use for downloading"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    String chromosome = "chr1"
    String genome_build = "b38"
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -eo pipefail

    # Download genetic map from GLIMPSE repository
    # Note: GLIMPSE genetic maps use chromosome numbers without 'chr' prefix,
    # but 1000 Genomes GRCh38 data uses 'chr' prefix, so we need to convert
    curl -sL -o "original.gmap.gz" \
      "https://raw.githubusercontent.com/odelaneau/GLIMPSE/master/maps/genetic_maps.~{genome_build}/~{chromosome}.~{genome_build}.gmap.gz"

    # Convert chromosome naming from "22" to "chr22" format to match GRCh38 convention
    # The genetic map has format: pos chr cM
    zcat original.gmap.gz | awk 'BEGIN{OFS="\t"} NR==1{print; next} {$2="chr"$2; print}' | gzip > "~{chromosome}.~{genome_build}.gmap.gz"

    rm original.gmap.gz

    echo "Downloaded and converted genetic map for ~{chromosome} (~{genome_build})"
    echo "First few lines:"
    zcat "~{chromosome}.~{genome_build}.gmap.gz" | head -5 || true
  >>>

  output {
    File genetic_map = "~{chromosome}.~{genome_build}.gmap.gz"
  }

  runtime {
    docker: "getwilds/awscli:2.27.49"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_glimpse2_reference_panel {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads and prepares a 1000 Genomes reference panel subset for GLIMPSE2 imputation. Downloads phased data for the specified chromosome and filters to a region."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        reference_vcf: "Reference panel VCF/BCF file for imputation",
        reference_vcf_index: "Index file for the reference panel",
        sites_vcf: "Sites-only VCF for genotype likelihood calculation",
        sites_vcf_index: "Index file for sites VCF"
    }
  }

  parameter_meta {
    chromosome: "Chromosome to download (e.g., chr1, chr22)"
    region: "Genomic region to extract (e.g., chr1:1-10000000). Must match the chromosome parameter."
    exclude_samples: "Comma-separated list of samples to exclude (useful for validation)"
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    String chromosome = "chr1"
    String region = "chr1:1-10000000"
    String exclude_samples = "NA12878"
    Int cpu_cores = 2
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    # Download 1000 Genomes high-coverage phased data for the specified chromosome
    # Files follow the naming pattern: CCDG_14151_B01_GRM_WGS_2020-08-05_chrN.filtered.shapeit2-duohmm-phased.vcf.gz
    echo "Downloading 1000 Genomes ~{chromosome} reference panel..."
    wget -q -O "1000GP.~{chromosome}.vcf.gz" \
      "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_~{chromosome}.filtered.shapeit2-duohmm-phased.vcf.gz"
    wget -q -O "1000GP.~{chromosome}.vcf.gz.tbi" \
      "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_~{chromosome}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi"

    # Create samples to exclude file
    echo "~{exclude_samples}" | tr ',' '\n' > exclude_samples.txt

    # Filter to region, normalize, keep only biallelic SNPs, and exclude validation samples
    echo "Processing reference panel for region ~{region}..."
    bcftools view -r "~{region}" "1000GP.~{chromosome}.vcf.gz" | \
      bcftools norm -m -any | \
      bcftools view -m 2 -M 2 -v snps -S ^exclude_samples.txt | \
      bcftools annotate -x ^INFO/AC,^INFO/AN,^FORMAT/GT -Ob -o reference_panel.bcf

    bcftools index reference_panel.bcf

    # Create sites-only VCF for GL calculation
    bcftools view -G -Oz -o reference_panel.sites.vcf.gz reference_panel.bcf
    bcftools index -t reference_panel.sites.vcf.gz

    # Clean up large intermediate files
    rm "1000GP.~{chromosome}.vcf.gz" "1000GP.~{chromosome}.vcf.gz.tbi"

    echo "Reference panel preparation complete"
    echo "Variants in panel: $(bcftools view -H reference_panel.bcf | wc -l)"
  >>>

  output {
    File reference_vcf = "reference_panel.bcf"
    File reference_vcf_index = "reference_panel.bcf.csi"
    File sites_vcf = "reference_panel.sites.vcf.gz"
    File sites_vcf_index = "reference_panel.sites.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/bcftools:1.19"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_glimpse2_test_gl_vcf {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads low-coverage sequencing data from 1000 Genomes and generates a VCF with genotype likelihoods for GLIMPSE2 imputation testing."
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        gl_vcf: "VCF file with genotype likelihoods (GL field) for imputation",
        gl_vcf_index: "Index file for the GL VCF"
    }
  }

  parameter_meta {
    chromosome: "Chromosome to download (e.g., chr1, chr22). Note: Phase 3 data uses numeric chromosome names (1-22)."
    region: "Genomic region to extract (e.g., chr1:1-10000000). Must match the chromosome parameter."
    sample_name: "Sample to extract from 1000 Genomes (must be in low-coverage dataset)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    String chromosome = "chr1"
    String region = "chr1:1-10000000"
    String sample_name = "NA12878"
    Int cpu_cores = 2
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    # Download 1000 Genomes Phase 3 low-coverage data for the specified chromosome
    # This data includes genotype likelihoods (GL field) from low-coverage sequencing
    # Note: Phase 3 files use numeric chromosome names (e.g., ALL.chr1.phase3...)
    echo "Downloading 1000 Genomes Phase 3 low-coverage data for ~{chromosome}..."

    # Extract region and sample, keeping GL field
    bcftools view \
      -r "~{region}" \
      -s "~{sample_name}" \
      --force-samples \
      "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.~{chromosome}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" | \
    bcftools annotate -x ^INFO/AF,^FORMAT/GT,^FORMAT/GL -Oz -o "~{sample_name}.gl.vcf.gz"

    bcftools index -t "~{sample_name}.gl.vcf.gz"

    echo "GL VCF download complete"
    echo "Variants with GL: $(bcftools view -H ~{sample_name}.gl.vcf.gz | wc -l)"
  >>>

  output {
    File gl_vcf = "~{sample_name}.gl.vcf.gz"
    File gl_vcf_index = "~{sample_name}.gl.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/bcftools:1.19"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task download_jcast_test_data {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Downloads example rMATS output files and Ensembl reference data for JCAST alternative splicing proteomics testing"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl"
    outputs: {
        rmats_output: "Tarball containing rMATS output files for JCAST testing",
        gtf_file: "Ensembl GTF annotation file (human chr15) required by JCAST",
        genome_fasta: "Ensembl genome FASTA file (human chr15) required by JCAST"
    }
  }

  parameter_meta {
    cpu_cores: "Number of CPU cores to use for downloading"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -eo pipefail

    # Base URL for JCAST test data
    BASE_URL="https://raw.githubusercontent.com/ed-lau/jcast/master/tests/data"

    # Create directory for rMATS test data
    mkdir -p rmats_test_output

    echo "Downloading rMATS test data from JCAST repository..."

    # Download each splice type file that JCAST expects
    for splice_type in SE MXE RI A3SS A5SS; do
      echo "Downloading ${splice_type}.MATS.JC.txt..."
      wget -q --no-check-certificate -O "rmats_test_output/${splice_type}.MATS.JC.txt" \
        "${BASE_URL}/rmats/${splice_type}.MATS.JC.txt" || echo "Warning: ${splice_type}.MATS.JC.txt not found"
    done

    # List downloaded rMATS files
    echo "Downloaded rMATS test files:"
    ls -la rmats_test_output/

    # Create tarball of rMATS test data
    tar -czf rmats_test_output.tar.gz rmats_test_output

    # Download Ensembl GTF file (JCAST requires Ensembl format with transcript_type attribute)
    echo "Downloading Ensembl GTF annotation file..."
    wget -q --no-check-certificate -O "test_reference.gtf" \
      "${BASE_URL}/genome/Homo_sapiens.GRCh38.89.chromosome.15.gtf"

    # Download Ensembl genome FASTA file
    echo "Downloading Ensembl genome FASTA file..."
    wget -q --no-check-certificate -O "test_reference.fa.gz" \
      "${BASE_URL}/genome/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz"
    gunzip test_reference.fa.gz

    echo "Test data preparation complete"
    ls -la
  >>>

  output {
    File rmats_output = "rmats_test_output.tar.gz"
    File gtf_file = "test_reference.gtf"
    File genome_fasta = "test_reference.fa"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
