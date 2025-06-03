version 1.0
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

struct SampleInfo {
    String omics_sample_name
    String molecular_id
    File r1_fastq
    File r2_fastq
}
 
workflow ww_vc_trio {
  input {
    Array[SampleInfo] samples
    File bedLocation
    String ref_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    File af_only_gnomad
    File af_only_gnomad_index
    String annovar_protocols
    String annovar_operation
  }

  # Docker containers this workflow has been designed for
  String GATKDocker = "getwilds/gatk:4.3.0.0"
  String bwaDocker = "getwilds/bwa:0.7.17"
  String bedtoolsDocker = "getwilds/bedtools:2.31.1" 
  String bcftoolsDocker = "getwilds/bcftools:1.19"
  String annovarDocker = "getwilds/annovar:~{ref_name}"
  String RDocker = "getwilds/consensus:0.1.1"

  Int bwaThreads = 16

  # Prepare bed file and check sorting
  call SortBed {
    input:
      unsorted_bed = bedLocation,
      ref_dict = ref_dict,
      docker = GATKDocker
  }
  
  scatter (sample in samples){
    String sampleName = sample.omics_sample_name
    String molecularID = sample.molecular_id
    File sampleR1 = sample.r1_fastq
    File sampleR2 = sample.r2_fastq

    String base_file_name = sampleName + "_" + molecularID + "." + ref_name

    # Map reads to reference directly from paired FASTQ files
    call BwaMem {
      input:
        r1_fastq = sampleR1,
        r2_fastq = sampleR2,
        sample_name = sampleName,
        library_name = molecularID,
        base_file_name = base_file_name,
        ref_fasta = ref_fasta,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        threads = bwaThreads,
        docker = bwaDocker
    }
  
    # Aggregate aligned+merged flowcell BAM files and mark duplicates
    call MarkDuplicates {
      input:
        input_bam = BwaMem.output_bam,
        input_bai = BwaMem.output_bai,
        output_bam_basename = base_file_name + ".aligned.duplicates_marked",
        metrics_filename = base_file_name + ".duplicate_metrics",
        docker = GATKDocker
    }

    # Generate the recalibration model by interval and apply it
    call ApplyBaseRecalibrator {
      input:
        input_bam = MarkDuplicates.output_bam,
        input_bam_index = MarkDuplicates.output_bai,
        base_file_name = base_file_name,
        intervals = SortBed.intervals,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = GATKDocker
    }

    call bedToolsQC {
      input: 
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        genome_sort_order = ApplyBaseRecalibrator.sortOrder,
        bed_file = SortBed.sorted_bed,
        base_file_name = base_file_name,
        docker = bedtoolsDocker
    }

    call CollectHsMetrics {
      input: 
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        base_file_name = base_file_name,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        intervals = SortBed.intervals,
        docker = GATKDocker
    }

    # Generate haplotype caller vcf
    call HaplotypeCaller {
      input:
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        input_bam_index = ApplyBaseRecalibrator.recalibrated_bai,
        intervals = SortBed.intervals,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_index = dbSNP_vcf_index,
        docker = GATKDocker
    }

    # Generate mutect2 vcf
    call Mutect2TumorOnly {
      input:
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        input_bam_index = ApplyBaseRecalibrator.recalibrated_bai,
        intervals = SortBed.intervals,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        genomeReference = af_only_gnomad,
        genomeReferenceIndex = af_only_gnomad_index,
        docker = GATKDocker
    }

    # Generate bcftools vcf
    call bcftoolsMpileup {
      input:
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        input_bam_index = ApplyBaseRecalibrator.recalibrated_bai,
        sorted_bed = SortBed.sorted_bed,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = bcftoolsDocker
    }

    # Annotate variants
    call annovar as annotateSAM {
      input:
        input_vcf = bcftoolsMpileup.output_vcf,
        ref_name = ref_name,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        docker = annovarDocker
    }

    # Annotate variants
    call annovar as annotateMutect {
      input:
        input_vcf = Mutect2TumorOnly.output_vcf,
        ref_name = ref_name,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        docker = annovarDocker
    }
    
    # Annotate variants
    call annovar as annotateHaplotype {
      input:
        input_vcf = HaplotypeCaller.output_vcf,
        ref_name = ref_name,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        docker = annovarDocker
    }

    call consensusProcessingR {
      input:
        GATKVars = annotateHaplotype.output_annotated_table,
        MutectVars = annotateMutect.output_annotated_table,
        SAMVars = annotateSAM.output_annotated_table,
        base_file_name = base_file_name,
        docker = RDocker
    }
  } # End scatter 

  # Outputs that will be retained when execution is complete
  output {
    Array[File] analysis_ready_bam = ApplyBaseRecalibrator.recalibrated_bam 
    Array[File] analysis_ready_bai = ApplyBaseRecalibrator.recalibrated_bai
    Array[File] GATK_vcf = HaplotypeCaller.output_vcf
    Array[File] SAM_vcf = bcftoolsMpileup.output_vcf
    Array[File] Mutect_Vcf = Mutect2TumorOnly.output_vcf
    Array[File] Mutect_VcfIndex = Mutect2TumorOnly.output_vcf_index
    Array[File] Mutect_AnnotatedVcf = annotateMutect.output_annotated_vcf
    Array[File] Mutect_AnnotatedTable = annotateMutect.output_annotated_table
    Array[File] GATK_annotated_vcf = annotateHaplotype.output_annotated_vcf
    Array[File] GATK_annotated = annotateHaplotype.output_annotated_table
    Array[File] SAM_annotated_vcf = annotateSAM.output_annotated_vcf
    Array[File] SAM_annotated = annotateSAM.output_annotated_table
    Array[File] panelQC = bedToolsQC.meanQC
    Array[File] PicardQC = CollectHsMetrics.picardMetrics
    Array[File] PicardQCpertarget = CollectHsMetrics.picardPerTarget
    Array[File] consensusVariants = consensusProcessingR.consensusTSV
  }
} # End workflow

#### TASK DEFINITIONS

# annotate with annovar
task annovar {
  input {
    File input_vcf
    String ref_name
    String annovar_protocols
    String annovar_operation
    String docker
  }

  String base_vcf_name = basename(input_vcf, ".vcf.gz")

  command <<<
    set -eo pipefail
    perl /annovar/table_annovar.pl "~{input_vcf}" /annovar/humandb/ \
      -buildver "~{ref_name}" \
      -outfile "~{base_vcf_name}" \
      -remove \
      -protocol "~{annovar_protocols}" \
      -operation "~{annovar_operation}" \
      -nastring . -vcfinput
    sed -i "s/Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\tOtherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10\tOtherinfo11\tOtherinfo12\tOtherinfo13/Otherinfo/g" "~{base_vcf_name}.~{ref_name}_multianno.txt"
  >>>

  output {
    File output_annotated_vcf = "~{base_vcf_name}.~{ref_name}_multianno.vcf"
    File output_annotated_table = "~{base_vcf_name}.~{ref_name}_multianno.txt"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2GB"
  }
}

# Generate Base Quality Score Recalibration (BQSR) model and apply it
task ApplyBaseRecalibrator {
  input {
    File input_bam
    File intervals 
    File input_bam_index
    String base_file_name
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Xms8g -Xmx8g" \
      BaseRecalibrator \
        -R "~{ref_fasta}" \
        -I "~{input_bam}" \
        -O "~{base_file_name}.recal_data.csv" \
        --known-sites "~{dbSNP_vcf}" \
        --known-sites ~{sep=" --known-sites " known_indels_sites_VCFs} \
        --intervals ~{intervals} \
        --interval-padding 100 \
        --verbosity WARNING
    gatk --java-options "-Xms48g -Xmx48g" \
      ApplyBQSR \
        -bqsr "~{base_file_name}.recal_data.csv" \
        -I "~{input_bam}" \
        -O "~{base_file_name}.recal.bam" \
        -R "~{ref_fasta}" \
        --intervals ~{intervals} \
        --interval-padding 100 \
        --verbosity WARNING
    # finds the current sort order of this bam file
    samtools view -H "~{base_file_name}.recal.bam" | grep @SQ|sed 's/@SQ\tSN:\|LN://g' > "~{base_file_name}.sortOrder.txt"
  >>>

  output {
    File recalibrated_bam = "~{base_file_name}.recal.bam"
    File recalibrated_bai = "~{base_file_name}.recal.bai"
    File sortOrder = "~{base_file_name}.sortOrder.txt"
  }

  runtime {
    memory: "36 GB"
    cpu: 1
    docker: docker
  }
}

# bcftools Mpileup variant calling
task bcftoolsMpileup {
  input {
    File input_bam
    File input_bam_index
    String base_file_name
    File sorted_bed
    File ref_dict
    File ref_fasta
    File ref_fasta_index
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
      "~{input_bam}" | bcftools call -Oz -mv \
          -o "~{base_file_name}.SAM.vcf.gz"
  >>>

  output {
    File output_vcf = "~{base_file_name}.SAM.vcf.gz"
  }

  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 2
  }
}

# use bedtools to find basic QC data
task bedToolsQC {
  input {
    File input_bam
    File bed_file
    File genome_sort_order
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    bedtools sort -g "~{genome_sort_order}" -i "~{bed_file}" > correctly.sorted.bed
    bedtools coverage -mean -sorted -g "~{genome_sort_order}" -a correctly.sorted.bed \
        -b "~{input_bam}" > "~{base_file_name}.bedtoolsQCMean.txt"
  >>>

  output {
    File meanQC = "~{base_file_name}.bedtoolsQCMean.txt"
  }

  runtime {
    docker: docker
    memory: "4 GB"
    cpu: 1
  }
}

# align to genome using paired FASTQ files
task BwaMem {
  input {
    File r1_fastq
    File r2_fastq
    String sample_name
    String library_name
    String base_file_name
    File ref_fasta
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    Int threads
    String docker
  }

  command <<<
    set -eo pipefail
    bwa mem \
      -t ~{threads - 1} \
      -R "@RG\tID:~{library_name}\tSM:~{sample_name}\tLB:~{library_name}\tPU:~{library_name}\tPL:ILLUMINA" \
      "~{ref_fasta}" "~{r1_fastq}" "~{r2_fastq}" | samtools sort -o "~{base_file_name}.aligned.bam"
    samtools index "~{base_file_name}.aligned.bam"
  >>>

  output {
    File output_bam = "~{base_file_name}.aligned.bam"
    File output_bai = "~{base_file_name}.aligned.bam.bai"
  }

  runtime {
    docker: docker
    memory: "32GB"
    cpu: "~{threads}"
  }
}

# get hybrid capture based QC metrics via Picard
task CollectHsMetrics {
  input {
    File input_bam
    String base_file_name
    File ref_fasta
    File ref_fasta_index
    File intervals
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Xms64g -Xmx64g" \
      CollectHsMetrics \
      --INPUT "~{input_bam}" \
      --OUTPUT "~{base_file_name}.picard.metrics.txt" \
      --REFERENCE_SEQUENCE "~{ref_fasta}" \
      --ALLELE_FRACTION 0.01 \
      --BAIT_INTERVALS "~{intervals}" \
      --TARGET_INTERVALS "~{intervals}" \
      --PER_TARGET_COVERAGE "~{base_file_name}.picard.pertarget.txt" \
      --VERBOSITY WARNING
  >>>

  output {
    File picardMetrics = "~{base_file_name}.picard.metrics.txt"
    File picardPerTarget = "~{base_file_name}.picard.pertarget.txt"
  }

  runtime {
    docker: docker
    cpu: 8
    memory: "64 GB"
  }
}

task consensusProcessingR {
  input {
    File GATKVars
    File SAMVars
    File MutectVars
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    Rscript /consensus-trio-unpaired.R "~{GATKVars}" "~{SAMVars}" "~{MutectVars}" "~{base_file_name}"
  >>>

  output {
    File consensusTSV = "~{base_file_name}.consensus.tsv"
  }

  runtime {
    cpu: 1
    memory: "8 GB"
    docker: docker
  }
}

# HaplotypeCaller per-sample
task HaplotypeCaller {
  input {
    File input_bam
    File input_bam_index
    String base_file_name
    File intervals
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File dbSNP_vcf
    File dbSNP_index
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Xms8g -Xmx8g" \
      HaplotypeCaller \
      -R "~{ref_fasta}" \
      -I "~{input_bam}" \
      -O "~{base_file_name}.GATK.vcf.gz" \
      --dbsnp "~{dbSNP_vcf}" \
      --intervals "~{intervals}" \
      --interval-padding 100 \
      --verbosity WARNING 
  >>>

  output {
    File output_vcf = "~{base_file_name}.GATK.vcf.gz"
    File output_vcf_index = "~{base_file_name}.GATK.vcf.gz.tbi"
  }

  runtime {
    docker: docker
    memory: "12 GB"
    cpu: 1
  }
}

# Mutect 2 calling
task Mutect2TumorOnly {
  input {
    File input_bam
    File input_bam_index
    String base_file_name
    File intervals
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File genomeReference
    File genomeReferenceIndex
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Xms16g -Xmx16g" \
      Mutect2 \
        -R "~{ref_fasta}" \
        -I "~{input_bam}" \
        -O preliminary.vcf.gz \
        --intervals "~{intervals}" \
        --interval-padding 100 \
        --germline-resource "~{genomeReference}" \
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
    File output_vcf = "~{base_file_name}.mutect2.vcf.gz"
    File output_vcf_index = "~{base_file_name}.mutect2.vcf.gz.tbi"
  }

  runtime {
    docker: docker
    memory: "24 GB"
    cpu: 1
  }
}

# Read unmapped BAM, convert to FASTQ
task SamToFastq {
  input {
    File input_bam
    String base_file_name
    String docker
  }

  # this now sorts first in case the original bam is not queryname sorted or mark dup spark will complain
  command <<<
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms12g -Xmx12g" \
      SortSam \
     --INPUT "~{input_bam}" \
     --OUTPUT sorted.bam \
     --SORT_ORDER queryname \
     --VERBOSITY WARNING
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms8g -Xmx8g" \
      SamToFastq \
      --INPUT sorted.bam \
      --FASTQ "~{base_file_name}.fastq" \
      --INTERLEAVE true \
      --INCLUDE_NON_PF_READS true \
      --VERBOSITY WARNING
  >>>

  output {
    File output_fastq = "~{base_file_name}.fastq"
  }

  runtime {
    docker: docker
    memory: "16 GB"
    cpu: 1
  }
}

# Prepare bed file and check sorting
task SortBed {
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

task MarkDuplicates {
  input {
    File input_bam
    File input_bai
    String output_bam_basename
    String metrics_filename
    String docker
  }

  command <<<
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms16g -Xmx16g" \
      MarkDuplicates \
      --INPUT "~{input_bam}" \
      --OUTPUT "~{output_bam_basename}.bam" \
      --METRICS_FILE "~{metrics_filename}" \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --VERBOSITY WARNING
    samtools index "~{output_bam_basename}.bam"
  >>>

  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bai = "~{output_bam_basename}.bam.bai"
    File duplicate_metrics = "~{metrics_filename}"
  }

  runtime {
    docker: docker
    memory: "24 GB"
    cpu: 1
  }
}

# Sort to queryname when needed for a dataset so that markduplicatesSpark can be used. 
task SortSam {
  input {
    File input_bam
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms12g -Xmx12g" \
      SortSam \
     --INPUT "~{input_bam}" \
     --OUTPUT "~{base_file_name}.sorted.bam" \
     --SORT_ORDER queryname \
     --VERBOSITY WARNING
  >>>

  output {
    File output_bam = "~{base_file_name}.sorted.bam"
  }

  runtime {
    docker: docker
    memory: "16 GB"
    cpu: 1
  }
}

