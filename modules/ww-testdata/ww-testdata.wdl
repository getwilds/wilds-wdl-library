## WILDS WDL module for downloading reference data for testing purposes.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task download_ref_data {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads reference genome and index files for WILDS WDL test runs"
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
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
    rm NA12878.bam NA12878.bam.bai NA12878_chr1.bam NA12878_chr1.bam.bai NA12878_24RG_small.hg38.bai
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
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
    # Download filtered known indels vcf from GATK
    bcftools view ~{if defined(region) then "-r " + region else ""} \
    https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
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

task download_test_transcriptome {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Download a small test transcriptome for RNA-seq quantification testing. NOTE: This uses GENCODE (Ensembl) annotations, while other ww-testdata tasks use NCBI RefSeq. For production use, ensure annotation consistency across your pipeline."
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Extract and clean a reference sequence region for saturation mutagenesis analysis"
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Create a test GDC manifest file with small open-access files for testing gdc-client downloads"
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
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Downloads ShapeMapper example data (TPP riboswitch) from the official repository for testing RNA structure analysis workflows"
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
    docker: "getwilds/awscli:2.27.49"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
