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
        annotsv_test_vcf: "Test VCF file for AnnotSV",
        validation_report: "Validation report summarizing all outputs"
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

  call download_annotsv_vcf { }

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
    annotsv_test_vcf = download_annotsv_vcf.test_vcf
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
    # Outputs from the AnnotSV test VCF download
    File annotsv_test_vcf = download_annotsv_vcf.test_vcf
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
    samtools view -@ ~{cpu_cores} -C -T ~{ref_fasta} -o NA12878_chr1.cram NA12878_chr1.bam
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
    cpu_cores: "Number of CPU cores to use for downloading and processing"
    memory_gb: "Memory allocation in GB for the task"
  }

  input {
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
  >>>

  output {
    File bam = "NA12878_chr1.bam"
    File bai = "NA12878_chr1.bam.bai"
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
    annotsv_test_vcf: "AnnotSV test VCF file to validate"
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
    File annotsv_test_vcf
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

    if [[ -f "~{annotsv_test_vcf}" && -s "~{annotsv_test_vcf}" ]]; then
      echo "AnnotSV test VCF: ~{annotsv_test_vcf} - PASSED" >> validation_report.txt
    else
      echo "AnnotSV test VCF: ~{annotsv_test_vcf} - MISSING OR EMPTY" >> validation_report.txt
      validation_passed=false
    fi

    {
      echo ""
      echo "=== Validation Summary ==="
      echo "Total files validated: 15"
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
