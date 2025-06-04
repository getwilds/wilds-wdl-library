## Consensus variant calling workflow for human panel/PCR-based targeted DNA sequencing.
## Input requirements:
## - Pair-end sequencing data in FASTQ format
## - Sample information provided as structs in the input JSON
##
## Output Files:
## - recalibrated bam and it's index
## - GATK vcf
## - samtools/bcftools vcf
## - Annovar annotated vcfs and tabular variant list for each variant caller
## - Basic QC stats from bedtools for mean coverage over regions in panel

version 1.0

struct SampleInfo {
    String omics_sample_name
    String molecular_id
    File r1_fastq
    File r2_fastq
}

workflow ww_leukemia {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Consensus variant calling workflow for human panel/PCR-based targeted DNA sequencing with focus on leukemia analysis"
    url: "https://github.com/getwilds/wilds-wdl-library"
    version: "1.0.0"
    outputs: {
        analysis_ready_bam: "Array of recalibrated BAM files ready for downstream analysis",
        analysis_ready_bai: "Array of index files for the recalibrated BAM files",
        gatk_vcf: "Array of variant calls from GATK HaplotypeCaller",
        sam_vcf: "Array of variant calls from samtools/bcftools mpileup",
        mutect_vcf: "Array of variant calls from GATK Mutect2 (tumor-only mode)",
        mutect_vcf_index: "Array of index files for Mutect2 VCF files",
        mutect_annotated_vcf: "Array of Mutect2 VCF files annotated with Annovar",
        mutect_annotated_table: "Array of Mutect2 variants in tabular format with annotations",
        gatk_annotated_vcf: "Array of GATK HaplotypeCaller VCF files annotated with Annovar",
        gatk_annotated_table: "Array of GATK HaplotypeCaller variants in tabular format with annotations",
        sam_annotated_vcf: "Array of samtools/bcftools VCF files annotated with Annovar",
        sam_annotated: "Array of samtools/bcftools variants in tabular format with annotations",
        panel_qc: "Array of quality control metrics from bedtools showing mean coverage over target regions",
        picard_qc: "Array of hybrid selection metrics from Picard CollectHsMetrics",
        picard_qc_per_target: "Array of per-target coverage metrics from Picard",
        consensus_variants: "Array of consensus variant calls combining results from all three variant callers"
    }
  }

  parameter_meta {
    bed_location: "BED file defining target regions for panel sequencing"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Index file (.fai) for the reference genome FASTA"
    ref_dict: "Sequence dictionary (.dict) for the reference genome"
    ref_amb: "BWA index file (.amb) for the reference genome"
    ref_ann: "BWA index file (.ann) for the reference genome"
    ref_bwt: "BWA index file (.bwt) for the reference genome"
    ref_pac: "BWA index file (.pac) for the reference genome"
    ref_sa: "BWA index file (.sa) for the reference genome"
    dbsnp_vcf: "dbSNP VCF file for base quality score recalibration and variant annotation"
    dbsnp_vcf_index: "Index file for the dbSNP VCF"
    af_only_gnomad: "gnomAD allele frequency VCF file for Mutect2 germline resource"
    af_only_gnomad_index: "Index file for the gnomAD allele frequency VCF"
    known_indels_sites_vcfs: "Array of VCF files containing known indel sites for base quality score recalibration"
    known_indels_sites_indices: "Array of index files for the known indels VCF files"
    samples: "Array of sample information structs containing sample names, molecular IDs, and paired FASTQ files"
    ref_name: "Reference genome build name (e.g., 'hg38', 'hg19') used for Annovar annotation"
    annovar_protocols: "Comma-separated list of Annovar annotation protocols to apply"
    annovar_operation: "Comma-separated list of Annovar operations corresponding to the protocols"
  }

  input {
    File bed_location
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File dbsnp_vcf
    File dbsnp_vcf_index
    File af_only_gnomad
    File af_only_gnomad_index
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices
    Array[SampleInfo] samples
    String ref_name
    String annovar_protocols
    String annovar_operation
  }

  # Docker containers this workflow has been designed for
  String gatk_docker = "getwilds/gatk:4.3.0.0"
  String bwa_docker = "getwilds/bwa:0.7.17"
  String bedtools_docker = "getwilds/bedtools:2.31.1"
  String bcftools_docker = "getwilds/bcftools:1.19"
  String annovar_docker = "getwilds/annovar:~{ref_name}"
  String r_docker = "getwilds/consensus:0.1.1"

  Int bwa_threads = 16

  # Prepare bed file and check sorting
  call sort_bed { input:
      unsorted_bed = bed_location,
      ref_dict = ref_dict,
      docker = gatk_docker
  }

  scatter (sample in samples){
    String sample_name = sample.omics_sample_name
    String molecular_id = sample.molecular_id
    File sample_r1 = sample.r1_fastq
    File sample_r2 = sample.r2_fastq

    String base_file_name = sample_name + "_" + molecular_id + "." + ref_name

    # Map reads to reference directly from paired FASTQ files
    call bwa_mem { input:
        r1_fastq = sample_r1,
        r2_fastq = sample_r2,
        sample_name = sample_name,
        library_name = molecular_id,
        base_file_name = base_file_name,
        ref_fasta = ref_fasta,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        threads = bwa_threads,
        docker = bwa_docker
    }

    # Aggregate aligned+merged flowcell BAM files and mark duplicates
    call mark_duplicates { input:
        raw_bam = bwa_mem.aligned_bam,
        raw_bai = bwa_mem.aligned_bai,
        markdup_bam_basename = base_file_name + ".aligned.duplicates_marked",
        metrics_filename = base_file_name + ".duplicate_metrics",
        docker = gatk_docker
    }

    # Generate the recalibration model by interval and apply it
    call apply_base_recal { input:
        aligned_bam = mark_duplicates.markdup_bam,
        aligned_bam_index = mark_duplicates.markdup_bai,
        base_file_name = base_file_name,
        intervals = sort_bed.intervals,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        known_indels_sites_vcfs = known_indels_sites_vcfs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = gatk_docker
    }

    call bedtools_qc { input:
        aligned_bam = apply_base_recal.recalibrated_bam,
        genome_sort_order = apply_base_recal.sort_order,
        bed_file = sort_bed.sorted_bed,
        base_file_name = base_file_name,
        docker = bedtools_docker
    }

    call collect_hs_metrics { input:
        aligned_bam = apply_base_recal.recalibrated_bam,
        base_file_name = base_file_name,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        intervals = sort_bed.intervals,
        docker = gatk_docker
    }

    # Generate haplotype caller vcf
    call haplotype_caller { input:
        aligned_bam = apply_base_recal.recalibrated_bam,
        aligned_bam_index = apply_base_recal.recalibrated_bai,
        intervals = sort_bed.intervals,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_index = dbsnp_vcf_index,
        docker = gatk_docker
    }

    # Generate mutect2 vcf
    call mutect2_tumoronly { input:
        aligned_bam = apply_base_recal.recalibrated_bam,
        aligned_bam_index = apply_base_recal.recalibrated_bai,
        intervals = sort_bed.intervals,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        genome_ref = af_only_gnomad,
        genome_ref_index = af_only_gnomad_index,
        docker = gatk_docker
    }

    # Generate bcftools vcf
    call bcftools_mpileup { input:
        aligned_bam = apply_base_recal.recalibrated_bam,
        aligned_bam_index = apply_base_recal.recalibrated_bai,
        sorted_bed = sort_bed.sorted_bed,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = bcftools_docker
    }

    # Annotate variants
    call annovar as annotateSAM { input:
        vcf_to_annotate = bcftools_mpileup.mpileup_vcf,
        ref_name = ref_name,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        docker = annovar_docker
    }

    # Annotate variants
    call annovar as annotateMutect { input:
        vcf_to_annotate = mutect2_tumoronly.mutect2_vcf,
        ref_name = ref_name,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        docker = annovar_docker
    }

    # Annotate variants
    call annovar as annotateHaplotype { input:
        vcf_to_annotate = haplotype_caller.haplotype_vcf,
        ref_name = ref_name,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        docker = annovar_docker
    }

    call consensus_processing { input:
        gatk_vars = annotateHaplotype.annotated_table,
        sam_vars = annotateSAM.annotated_table,
        mutect_vars = annotateMutect.annotated_table,
        base_file_name = base_file_name,
        docker = r_docker
    }
  } # End scatter

  # Outputs that will be retained when execution is complete
  output {
    Array[File] analysis_ready_bam = apply_base_recal.recalibrated_bam
    Array[File] analysis_ready_bai = apply_base_recal.recalibrated_bai
    Array[File] gatk_vcf = haplotype_caller.haplotype_vcf
    Array[File] sam_vcf = bcftools_mpileup.mpileup_vcf
    Array[File] mutect_vcf = mutect2_tumoronly.mutect2_vcf
    Array[File] mutect_vcf_index = mutect2_tumoronly.mutect2_vcf_index
    Array[File] mutect_annotated_vcf = annotateMutect.annotated_vcf
    Array[File] mutect_annotated_table = annotateMutect.annotated_table
    Array[File] gatk_annotated_vcf = annotateHaplotype.annotated_vcf
    Array[File] gatk_annotated_table = annotateHaplotype.annotated_table
    Array[File] sam_annotated_vcf = annotateSAM.annotated_vcf
    Array[File] sam_annotated = annotateSAM.annotated_table
    Array[File] panel_qc = bedtools_qc.mean_qc
    Array[File] picard_qc = collect_hs_metrics.picard_metrics
    Array[File] picard_qc_per_target = collect_hs_metrics.picard_pertarget
    Array[File] consensus_variants = consensus_processing.consensus_tsv
  }
} # End workflow

#### TASK DEFINITIONS

# annotate with annovar
task annovar {
  meta {
    description: "Annotate variants using Annovar with customizable protocols and operations"
    outputs: {
        annotated_vcf: "VCF file with Annovar annotations added",
        annotated_table: "Tab-delimited table with variant annotations"
    }
  }

  parameter_meta {
    vcf_to_annotate: "Input VCF file to be annotated"
    ref_name: "Reference genome build name for Annovar annotation (e.g., 'hg38', 'hg19')"
    annovar_protocols: "Comma-separated list of annotation protocols to apply"
    annovar_operation: "Comma-separated list of operations corresponding to the protocols"
    docker: "Docker container image for Annovar"
  }

  input {
    File vcf_to_annotate
    String ref_name
    String annovar_protocols
    String annovar_operation
    String docker
  }

  String base_vcf_name = basename(vcf_to_annotate, ".vcf.gz")

  command <<<
    set -eo pipefail
    perl /annovar/table_annovar.pl "~{vcf_to_annotate}" /annovar/humandb/ \
      -buildver "~{ref_name}" \
      -outfile "~{base_vcf_name}" \
      -remove \
      -protocol "~{annovar_protocols}" \
      -operation "~{annovar_operation}" \
      -nastring . -vcfinput
    sed -i "s/Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\tOtherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10\tOtherinfo11\tOtherinfo12\tOtherinfo13/Otherinfo/g" "~{base_vcf_name}.~{ref_name}_multianno.txt"
  >>>

  output {
    File annotated_vcf = "~{base_vcf_name}.~{ref_name}_multianno.vcf"
    File annotated_table = "~{base_vcf_name}.~{ref_name}_multianno.txt"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2GB"
  }
}

# Generate Base Quality Score Recalibration (BQSR) model and apply it
task apply_base_recal {
  meta {
    description: "Generate Base Quality Score Recalibration (BQSR) model and apply it to improve base quality scores"
    outputs: {
        recalibrated_bam: "BAM file with recalibrated base quality scores",
        recalibrated_bai: "Index file for the recalibrated BAM",
        sort_order: "Text file containing the chromosome sort order of the BAM file"
    }
  }

  parameter_meta {
    aligned_bam: "Input aligned BAM file to be recalibrated"
    aligned_bam_index: "Index file for the input BAM"
    intervals: "Interval list file defining target regions"
    dbsnp_vcf: "dbSNP VCF file for known variant sites"
    dbsnp_vcf_index: "Index file for the dbSNP VCF"
    ref_dict: "Reference genome sequence dictionary"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Index file for the reference FASTA"
    known_indels_sites_vcfs: "Array of VCF files with known indel sites"
    known_indels_sites_indices: "Array of index files for known indels VCFs"
    base_file_name: "Base name for output files"
    docker: "Docker container image for GATK"
  }

  input {
    File aligned_bam
    File aligned_bam_index
    File intervals
    File dbsnp_vcf
    File dbsnp_vcf_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Xms8g -Xmx8g" \
      BaseRecalibrator \
        -R "~{ref_fasta}" \
        -I "~{aligned_bam}" \
        -O "~{base_file_name}.recal_data.csv" \
        --known-sites "~{dbsnp_vcf}" \
        --known-sites ~{sep=" --known-sites " known_indels_sites_vcfs} \
        --intervals "~{intervals}" \
        --interval-padding 100 \
        --verbosity WARNING
    gatk --java-options "-Xms48g -Xmx48g" \
      ApplyBQSR \
        -bqsr "~{base_file_name}.recal_data.csv" \
        -I "~{aligned_bam}" \
        -O "~{base_file_name}.recal.bam" \
        -R "~{ref_fasta}" \
        --intervals "~{intervals}" \
        --interval-padding 100 \
        --verbosity WARNING
    # finds the current sort order of this bam file
    samtools view -H "~{base_file_name}.recal.bam" | \
      grep @SQ|sed 's/@SQ\tSN:\|LN://g' > "~{base_file_name}.sortOrder.txt"
  >>>

  output {
    File recalibrated_bam = "~{base_file_name}.recal.bam"
    File recalibrated_bai = "~{base_file_name}.recal.bai"
    File sort_order = "~{base_file_name}.sortOrder.txt"
  }

  runtime {
    memory: "36 GB"
    cpu: 1
    docker: docker
  }
}

# bcftools Mpileup variant calling
task bcftools_mpileup {
  meta {
    description: "Call variants using bcftools mpileup and call pipeline"
    outputs: {
        mpileup_vcf: "Compressed VCF file containing variant calls from bcftools"
    }
  }

  parameter_meta {
    aligned_bam: "Input aligned BAM file"
    aligned_bam_index: "Index file for the input BAM"
    sorted_bed: "BED file defining target regions for variant calling"
    ref_dict: "Reference genome sequence dictionary"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Index file for the reference FASTA"
    base_file_name: "Base name for output files"
    docker: "Docker container image for bcftools"
  }

  input {
    File aligned_bam
    File aligned_bam_index
    File sorted_bed
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    bcftools mpileup \
      --max-depth 10000 \
      --max-idepth 10000 \
      --annotate "FORMAT/AD,FORMAT/DP" \
      --fasta-ref "~{ref_fasta}" \
      --regions-file "~{sorted_bed}" \
      --ignore-RG \
      --no-BAQ \
      "~{aligned_bam}" | bcftools call -Oz -mv \
          -o "~{base_file_name}.SAM.vcf.gz"
  >>>

  output {
    File mpileup_vcf = "~{base_file_name}.SAM.vcf.gz"
  }

  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 2
  }
}

# use bedtools to find basic QC data
task bedtools_qc {
  meta {
    description: "Generate quality control metrics using bedtools to calculate mean coverage over target regions"
    outputs: {
        mean_qc: "Tab-delimited file with mean coverage statistics for each target region"
    }
  }

  parameter_meta {
    aligned_bam: "Input aligned BAM file"
    bed_file: "BED file defining target regions"
    genome_sort_order: "File containing chromosome sort order information"
    base_file_name: "Base name for output files"
    docker: "Docker container image for bedtools"
  }

  input {
    File aligned_bam
    File bed_file
    File genome_sort_order
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    bedtools sort -g "~{genome_sort_order}" -i "~{bed_file}" > correctly.sorted.bed
    bedtools coverage -mean -sorted -g "~{genome_sort_order}" -a correctly.sorted.bed \
        -b "~{aligned_bam}" > "~{base_file_name}.bedtoolsQCMean.txt"
  >>>

  output {
    File mean_qc = "~{base_file_name}.bedtoolsQCMean.txt"
  }

  runtime {
    docker: docker
    memory: "4 GB"
    cpu: 1
  }
}

# align to genome using paired FASTQ files
task bwa_mem {
  meta {
    description: "Align paired-end reads to reference genome using BWA-MEM algorithm"
    outputs: {
        aligned_bam: "Coordinate-sorted BAM file containing aligned reads",
        aligned_bai: "Index file for the aligned BAM"
    }
  }

  parameter_meta {
    r1_fastq: "Read 1 FASTQ file"
    r2_fastq: "Read 2 FASTQ file"
    ref_fasta: "Reference genome FASTA file"
    ref_amb: "BWA index file (.amb) for the reference genome"
    ref_ann: "BWA index file (.ann) for the reference genome"
    ref_bwt: "BWA index file (.bwt) for the reference genome"
    ref_pac: "BWA index file (.pac) for the reference genome"
    ref_sa: "BWA index file (.sa) for the reference genome"
    sample_name: "Sample name for read group assignment"
    library_name: "Library name for read group assignment"
    base_file_name: "Base name for output files"
    docker: "Docker container image for BWA and samtools"
    threads: "Number of threads to use for alignment"
  }

  input {
    File r1_fastq
    File r2_fastq
    File ref_fasta
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    String sample_name
    String library_name
    String base_file_name
    String docker
    Int threads
  }

  command <<<
    set -eo pipefail
    bwa mem \
      -t ~{threads - 1} \
      -R "@RG\tID:~{library_name}\tSM:~{sample_name}\tLB:~{library_name}\tPU:~{library_name}\tPL:ILLUMINA" \
      "~{ref_fasta}" "~{r1_fastq}" "~{r2_fastq}" | \
      samtools sort -o "~{base_file_name}.aligned.bam"
    samtools index "~{base_file_name}.aligned.bam"
  >>>

  output {
    File aligned_bam = "~{base_file_name}.aligned.bam"
    File aligned_bai = "~{base_file_name}.aligned.bam.bai"
  }

  runtime {
    docker: docker
    memory: "32GB"
    cpu: threads
  }
}

# get hybrid capture based QC metrics via Picard
task collect_hs_metrics {
  meta {
    description: "Collect hybrid selection metrics using Picard to assess target enrichment performance"
    outputs: {
        picard_metrics: "Picard hybrid selection metrics summary file",
        picard_pertarget: "Per-target coverage metrics file"
    }
  }

  parameter_meta {
    aligned_bam: "Input aligned BAM file"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Index file for the reference FASTA"
    intervals: "Interval list file defining bait and target regions"
    base_file_name: "Base name for output files"
    docker: "Docker container image for GATK/Picard"
  }

  input {
    File aligned_bam
    File ref_fasta
    File ref_fasta_index
    File intervals
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Xms64g -Xmx64g" \
      CollectHsMetrics \
      --INPUT "~{aligned_bam}" \
      --OUTPUT "~{base_file_name}.picard.metrics.txt" \
      --REFERENCE_SEQUENCE "~{ref_fasta}" \
      --ALLELE_FRACTION 0.01 \
      --BAIT_INTERVALS "~{intervals}" \
      --TARGET_INTERVALS "~{intervals}" \
      --PER_TARGET_COVERAGE "~{base_file_name}.picard.pertarget.txt" \
      --VERBOSITY WARNING
  >>>

  output {
    File picard_metrics = "~{base_file_name}.picard.metrics.txt"
    File picard_pertarget = "~{base_file_name}.picard.pertarget.txt"
  }

  runtime {
    docker: docker
    cpu: 8
    memory: "64 GB"
  }
}

task consensus_processing {
  meta {
    description: "Generate consensus variant calls by combining results from multiple variant callers"
    outputs: {
        consensus_tsv: "Tab-separated file containing consensus variant calls with evidence from all callers"
    }
  }

  parameter_meta {
    gatk_vars: "Annotated variant table from GATK HaplotypeCaller"
    sam_vars: "Annotated variant table from samtools/bcftools"
    mutect_vars: "Annotated variant table from GATK Mutect2"
    base_file_name: "Base name for output files"
    docker: "Docker container image with R and consensus calling script"
  }

  input {
    File gatk_vars
    File sam_vars
    File mutect_vars
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    Rscript /consensus-trio-unpaired.R \
      "~{gatk_vars}" "~{sam_vars}" "~{mutect_vars}" "~{base_file_name}"
  >>>

  output {
    File consensus_tsv = "~{base_file_name}.consensus.tsv"
  }

  runtime {
    cpu: 1
    memory: "8 GB"
    docker: docker
  }
}

# HaplotypeCaller per-sample
task haplotype_caller {
  meta {
    description: "Call germline variants using GATK HaplotypeCaller"
    outputs: {
        haplotype_vcf: "Compressed VCF file containing germline variant calls",
        haplotype_vcf_index: "Index file for the VCF output"
    }
  }

  parameter_meta {
    aligned_bam: "Input aligned BAM file"
    aligned_bam_index: "Index file for the input BAM"
    intervals: "Interval list file defining target regions"
    ref_dict: "Reference genome sequence dictionary"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Index file for the reference FASTA"
    dbsnp_vcf: "dbSNP VCF file for variant annotation"
    dbsnp_index: "Index file for the dbSNP VCF"
    base_file_name: "Base name for output files"
    docker: "Docker container image for GATK"
  }

  input {
    File aligned_bam
    File aligned_bam_index
    File intervals
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File dbsnp_vcf
    File dbsnp_index
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Xms8g -Xmx8g" \
      HaplotypeCaller \
      -R "~{ref_fasta}" \
      -I "~{aligned_bam}" \
      -O "~{base_file_name}.GATK.vcf.gz" \
      --dbsnp "~{dbsnp_vcf}" \
      --intervals "~{intervals}" \
      --interval-padding 100 \
      --verbosity WARNING 
  >>>

  output {
    File haplotype_vcf = "~{base_file_name}.GATK.vcf.gz"
    File haplotype_vcf_index = "~{base_file_name}.GATK.vcf.gz.tbi"
  }

  runtime {
    docker: docker
    memory: "12 GB"
    cpu: 1
  }
}

# Mutect 2 calling
task mutect2_tumoronly {
  meta {
    description: "Call somatic variants using GATK Mutect2 in tumor-only mode with filtering"
    outputs: {
        mutect2_vcf: "Compressed VCF file containing filtered somatic variant calls",
        mutect2_vcf_index: "Index file for the Mutect2 VCF output"
    }
  }

  parameter_meta {
    aligned_bam: "Input aligned BAM file"
    aligned_bam_index: "Index file for the input BAM"
    intervals: "Interval list file defining target regions"
    ref_dict: "Reference genome sequence dictionary"
    ref_fasta: "Reference genome FASTA file"
    ref_fasta_index: "Index file for the reference FASTA"
    genome_ref: "gnomAD population allele frequency VCF for germline resource"
    genome_ref_index: "Index file for the gnomAD VCF"
    base_file_name: "Base name for output files"
    docker: "Docker container image for GATK"
  }

  input {
    File aligned_bam
    File aligned_bam_index
    File intervals
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File genome_ref
    File genome_ref_index
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Xms16g -Xmx16g" \
      Mutect2 \
        -R "~{ref_fasta}" \
        -I "~{aligned_bam}" \
        -O preliminary.vcf.gz \
        --intervals "~{intervals}" \
        --interval-padding 100 \
        --germline-resource "~{genome_ref}" \
        --verbosity WARNING
    gatk --java-options "-Xms16g -Xmx16g" \
      FilterMutectCalls \
        -V preliminary.vcf.gz \
        -O "~{base_file_name}.mutect2.vcf.gz" \
        -R "~{ref_fasta}" \
        --stats preliminary.vcf.gz.stats \
        --verbosity WARNING
  >>>

  output {
    File mutect2_vcf = "~{base_file_name}.mutect2.vcf.gz"
    File mutect2_vcf_index = "~{base_file_name}.mutect2.vcf.gz.tbi"
  }

  runtime {
    docker: docker
    memory: "24 GB"
    cpu: 1
  }
}

# Prepare bed file and check sorting
task sort_bed {
  meta {
    description: "Sort BED file and convert to GATK interval list format for consistent interval processing"
    outputs: {
        intervals: "GATK interval list file sorted by genomic coordinates",
        sorted_bed: "BED file sorted by genomic coordinates"
    }
  }

  parameter_meta {
    unsorted_bed: "Input BED file defining target regions (may be unsorted)"
    ref_dict: "Reference genome sequence dictionary for interval list conversion"
    docker: "Docker container image for GATK"
  }

  input {
    File unsorted_bed
    File ref_dict
    String docker
  }

  command <<<
    set -eo pipefail
    sort -k1,1V -k2,2n -k3,3n "~{unsorted_bed}" > sorted.bed
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g -Xmx4g" \
      BedToIntervalList \
      --INPUT sorted.bed \
      --OUTPUT sorted.interval_list \
      --SEQUENCE_DICTIONARY "~{ref_dict}"
  >>>

  output {
    File intervals = "sorted.interval_list"
    File sorted_bed = "sorted.bed"
  }

  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 1
  }
}

task mark_duplicates {
  meta {
    description: "Mark duplicate reads in aligned BAM file to improve variant calling accuracy"
    outputs: {
        markdup_bam: "BAM file with duplicate reads marked",
        markdup_bai: "Index file for the duplicate-marked BAM",
        duplicate_metrics: "Metrics file containing duplicate marking statistics"
    }
  }

  parameter_meta {
    raw_bam: "Input aligned BAM file"
    raw_bai: "Index file for the input BAM"
    markdup_bam_basename: "Base name for the output duplicate-marked BAM file"
    metrics_filename: "Name for the duplicate metrics output file"
    docker: "Docker container image for GATK"
  }

  input {
    File raw_bam
    File raw_bai
    String markdup_bam_basename
    String metrics_filename
    String docker
  }

  command <<<
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms16g -Xmx16g" \
      MarkDuplicates \
      --INPUT "~{raw_bam}" \
      --OUTPUT "~{markdup_bam_basename}.bam" \
      --METRICS_FILE "~{metrics_filename}" \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --VERBOSITY WARNING
    samtools index "~{markdup_bam_basename}.bam"
  >>>

  output {
    File markdup_bam = "~{markdup_bam_basename}.bam"
    File markdup_bai = "~{markdup_bam_basename}.bam.bai"
    File duplicate_metrics = "~{metrics_filename}"
  }

  runtime {
    docker: docker
    memory: "24 GB"
    cpu: 1
  }
}
