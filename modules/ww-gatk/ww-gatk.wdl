## WILDS WDL module for GATK variant calling and analysis tasks.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/add-gatk/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct GatkSample {
    String name
    File bam
    File bai
}

workflow gatk_example {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow demonstrating GATK variant calling tasks"
    url: "https://github.com/getwilds/wilds-wdl-library/modules/ww-gatk"
    outputs: {
        recalibrated_bams: "Array of base quality recalibrated BAM files",
        recalibrated_bais: "Array of index files for recalibrated BAMs",
        haplotype_vcfs: "Array of germline variant calls from HaplotypeCaller",
        mutect2_vcfs: "Array of somatic variant calls from Mutect2",
        wgs_metrics: "Array of whole genome sequencing metrics files",
        validation_report: "Validation report summarizing all outputs"
    }
  }

  parameter_meta {
    samples: "Array of sample information objects containing BAM files and indices"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Reference genome FASTA index (.fai)"
    dbsnp_vcf: "dbSNP VCF file for variant annotation and base recalibration"
    known_indels_sites_vcfs: "Array of VCF files with known indel sites for BQSR"
    gnomad_vcf: "gnomAD population allele frequency VCF for Mutect2"
    intervals: "Optional interval list file defining target regions"
  }

  input {
    Array[GatkSample]? samples
    File? reference_fasta
    File? reference_fasta_index
    File? dbsnp_vcf
    Array[File]? known_indels_sites_vcfs
    File? gnomad_vcf
    File? intervals
  }

  # Determine which genome files to use
  if (!defined(reference_fasta) || !defined(reference_fasta_index)) {
    call ww_testdata.download_ref_data { }
  }
  File genome_fasta = select_first([reference_fasta, download_ref_data.fasta])
  File genome_fasta_index = select_first([reference_fasta_index, download_ref_data.fasta_index])

  # Create reference dictionary
  call create_sequence_dictionary { input:
      reference_fasta = genome_fasta
  }

  # Determine which dbSNP vcf to use
  if (!defined(dbsnp_vcf)) {
    call ww_testdata.download_dbsnp_vcf { input:
      region = "NC_000001.11:1-10000000",
      filter_name = "chr1"
    }
  }
  File final_dbsnp = select_first([dbsnp_vcf, download_dbsnp_vcf.dbsnp_vcf])

  # Determine which known indels vcf to use
  if (!defined(known_indels_sites_vcfs)) {
    call ww_testdata.download_known_indels_vcf { input:
      region = "chr1:1-10000000",
      filter_name = "chr1"
    }
  }
  Array[File] final_known_indels = select_first([known_indels_sites_vcfs, select_all([download_known_indels_vcf.known_indels_vcf])])

  # Determine which gnomad vcf to use
  if (!defined(gnomad_vcf)) {
    call ww_testdata.download_gnomad_vcf { input:
      region = "chr1:1-10000000",
      filter_name = "chr1"
    }
  }
  File final_gnomad = select_first([gnomad_vcf, download_gnomad_vcf.gnomad_vcf])

  # If no samples provided, download demonstration data
  if (!defined(samples)) {
    call ww_testdata.download_bam_data { }
  }

  # Create samples array - either from input or from BWA alignment
  Array[GatkSample] final_samples = if defined(samples) then select_first([samples]) else [
    {
      "name": "demo_sample",
      "bam": select_first([download_bam_data.bam]),
      "bai": select_first([download_bam_data.bai])
    }
  ]

  scatter (sample in final_samples) {
    call base_recalibrator { input:
        aligned_bam = sample.bam,
        aligned_bam_index = sample.bai,
        intervals = intervals,
        dbsnp_vcf = final_dbsnp,
        reference_fasta = genome_fasta,
        reference_fasta_index = genome_fasta_index,
        reference_dict = create_sequence_dictionary.sequence_dict,
        known_indels_sites_vcfs = final_known_indels,
        output_basename = sample.name + ".recalibrated"
    }

    call haplotype_caller { input:
        bam = base_recalibrator.recalibrated_bam,
        bam_index = base_recalibrator.recalibrated_bai,
        intervals = intervals,
        reference_fasta = genome_fasta,
        reference_fasta_index = genome_fasta_index,
        reference_dict = create_sequence_dictionary.sequence_dict,
        dbsnp_vcf = final_dbsnp,
        output_basename = sample.name + ".haplotypecaller"
    }

    call mutect2 { input:
        bam = base_recalibrator.recalibrated_bam,
        bam_index = base_recalibrator.recalibrated_bai,
        intervals = intervals,
        reference_fasta = genome_fasta,
        reference_fasta_index = genome_fasta_index,
        reference_dict = create_sequence_dictionary.sequence_dict,
        gnomad_vcf = final_gnomad,
        output_basename = sample.name + ".mutect2"
    }

    call collect_wgs_metrics { input:
        bam = base_recalibrator.recalibrated_bam,
        bam_index = base_recalibrator.recalibrated_bai,
        reference_fasta = genome_fasta,
        reference_fasta_index = genome_fasta_index,
        intervals = intervals,
        output_basename = sample.name + ".wgsmetrics"
    }
  }

  call validate_outputs { input:
      recalibrated_bams = base_recalibrator.recalibrated_bam,
      recalibrated_bais = base_recalibrator.recalibrated_bai,
      haplotype_vcfs = haplotype_caller.vcf,
      mutect2_vcfs = mutect2.vcf,
      wgs_metrics = collect_wgs_metrics.metrics_file
  }

  output {
    Array[File] recalibrated_bams = base_recalibrator.recalibrated_bam
    Array[File] recalibrated_bais = base_recalibrator.recalibrated_bai
    Array[File] haplotype_vcfs = haplotype_caller.vcf
    Array[File] mutect2_vcfs = mutect2.vcf
    Array[File] wgs_metrics = collect_wgs_metrics.metrics_file
    File validation_report = validate_outputs.report
  }
}

task base_recalibrator {
  meta {
    description: "Generate Base Quality Score Recalibration (BQSR) model and apply it to improve base quality scores"
    outputs: {
        recalibrated_bam: "BAM file with recalibrated base quality scores",
        recalibrated_bai: "Index file for the recalibrated BAM",
        recalibration_report: "Base recalibration report table"
    }
  }

  parameter_meta {
    aligned_bam: "Input aligned BAM file to be recalibrated"
    aligned_bam_index: "Index file for the input BAM"
    intervals: "Optional interval list file defining target regions"
    dbsnp_vcf: "dbSNP VCF file for known variant sites"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    reference_dict: "Reference genome sequence dictionary"
    known_indels_sites_vcfs: "Array of VCF files with known indel sites"
    output_basename: "Base name for output files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File aligned_bam
    File aligned_bam_index
    File? intervals
    File dbsnp_vcf
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] known_indels_sites_vcfs
    String output_basename
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail

    # Add local symbolic link for reference fasta and dict
    # If soft links aren't allowed on your HPC system, copy them locally instead
    ln -s "~{reference_fasta}" "~{basename(reference_fasta)}"
    ln -s "~{reference_fasta_index}" "~{basename(reference_fasta_index)}"
    ln -s "~{reference_dict}" "~{basename(reference_dict)}"
    
    # Generate vcf index files using GATK
    gatk IndexFeatureFile -I "~{dbsnp_vcf}"
    known_vcfs=(~{sep=" " known_indels_sites_vcfs})
    for known_vcf in "${known_vcfs[@]}"; do
      gatk IndexFeatureFile -I "${known_vcf}"
    done

    # Generate Base Recalibration Table
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      BaseRecalibrator \
      -R "~{basename(reference_fasta)}" \
      -I "~{aligned_bam}" \
      -O "~{output_basename}.recal_data.table" \
      --known-sites "~{dbsnp_vcf}" \
      --known-sites ~{sep=" --known-sites " known_indels_sites_vcfs} \
      ~{if defined(intervals) then "--intervals " + intervals else ""} \
      --verbosity WARNING

    # Apply Base Quality Score Recalibration
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      ApplyBQSR \
      -R "~{basename(reference_fasta)}" \
      -I "~{aligned_bam}" \
      -bqsr "~{output_basename}.recal_data.table" \
      -O "~{output_basename}.bam" \
      ~{if defined(intervals) then "--intervals " + intervals else ""} \
      --verbosity WARNING

    # Index the recalibrated BAM
    samtools index "~{output_basename}.bam"
  >>>

  output {
    File recalibrated_bam = "~{output_basename}.bam"
    File recalibrated_bai = "~{output_basename}.bam.bai"
    File recalibration_report = "~{output_basename}.recal_data.table"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task haplotype_caller {
  meta {
    description: "Call germline variants using GATK HaplotypeCaller"
    outputs: {
        vcf: "Compressed VCF file containing germline variant calls",
        vcf_index: "Index file for the VCF output"
    }
  }

  parameter_meta {
    bam: "Input aligned BAM file"
    bam_index: "Index file for the input BAM"
    intervals: "Optional interval list file defining target regions"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    reference_dict: "Reference genome sequence dictionary"
    dbsnp_vcf: "dbSNP VCF file for variant annotation"
    output_basename: "Base name for output files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File bam
    File bam_index
    File? intervals
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File dbsnp_vcf
    String output_basename
    Int memory_gb = 16
    Int cpu_cores = 4
  }

  command <<<
    set -eo pipefail
    
    # Add local symbolic link for reference fasta and dict
    # If soft links aren't allowed on your HPC system, copy them locally instead
    ln -s "~{reference_fasta}" "~{basename(reference_fasta)}"
    ln -s "~{reference_fasta_index}" "~{basename(reference_fasta_index)}"
    ln -s "~{reference_dict}" "~{basename(reference_dict)}"

    # Create index for dbSNP vcf
    gatk IndexFeatureFile -I "~{dbsnp_vcf}"

    # Run HaplotypeCaller    
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      HaplotypeCaller \
      -R "~{basename(reference_fasta)}" \
      -I "~{bam}" \
      -O "~{output_basename}.vcf.gz" \
      --dbsnp "~{dbsnp_vcf}" \
      ~{if defined(intervals) then "--intervals " + intervals else ""} \
      ~{if defined(intervals) then "--interval-padding 100" else ""} \
      --verbosity WARNING
  >>>

  output {
    File vcf = "~{output_basename}.vcf.gz"
    File vcf_index = "~{output_basename}.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task mutect2 {
  meta {
    description: "Call somatic variants using GATK Mutect2 in tumor-only mode with filtering"
    outputs: {
        vcf: "Compressed VCF file containing filtered somatic variant calls",
        vcf_index: "Index file for the Mutect2 VCF output",
        stats_file: "Mutect2 statistics file",
        f1r2_counts: "F1R2 counts for filtering"
    }
  }

  parameter_meta {
    bam: "Input aligned BAM file"
    bam_index: "Index file for the input BAM"
    intervals: "Optional interval list file defining target regions"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    reference_dict: "Reference genome sequence dictionary"
    gnomad_vcf: "gnomAD population allele frequency VCF for germline resource"
    output_basename: "Base name for output files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File bam
    File bam_index
    File? intervals
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File gnomad_vcf
    String output_basename
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail
    
    # Add local symbolic link for reference fasta and dict
    # If soft links aren't allowed on your HPC system, copy them locally instead
    ln -s "~{reference_fasta}" "~{basename(reference_fasta)}"
    ln -s "~{reference_fasta_index}" "~{basename(reference_fasta_index)}"
    ln -s "~{reference_dict}" "~{basename(reference_dict)}"

    # Index gnomad VCF
    gatk IndexFeatureFile -I "~{gnomad_vcf}"

    # Run Mutect2
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      Mutect2 \
      -R "~{basename(reference_fasta)}" \
      -I "~{bam}" \
      -O "~{output_basename}.unfiltered.vcf.gz" \
      ~{if defined(intervals) then "--intervals " + intervals else ""} \
      ~{if defined(intervals) then "--interval-padding 100" else ""} \
      --germline-resource "~{gnomad_vcf}" \
      --f1r2-tar-gz "~{output_basename}.f1r2.tar.gz" \
      --verbosity WARNING

    # Filter Mutect2 calls
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      FilterMutectCalls \
      -V "~{output_basename}.unfiltered.vcf.gz" \
      -R "~{basename(reference_fasta)}" \
      -O "~{output_basename}.vcf.gz" \
      --stats "~{output_basename}.unfiltered.vcf.gz.stats" \
      --verbosity WARNING
  >>>

  output {
    File vcf = "~{output_basename}.vcf.gz"
    File vcf_index = "~{output_basename}.vcf.gz.tbi"
    File stats_file = "~{output_basename}.unfiltered.vcf.gz.stats"
    File f1r2_counts = "~{output_basename}.f1r2.tar.gz"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task create_sequence_dictionary {
  meta {
    description: "Create a sequence dictionary file from a reference FASTA file"
    outputs: {
        sequence_dict: "Sequence dictionary file (.dict) for the reference genome"
    }
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File reference_fasta
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  String dict_basename = basename(reference_fasta, ".fa")

  command <<<
    set -eo pipefail
    
    gatk --java-options "-Xms~{memory_gb - 2}g -Xmx~{memory_gb - 1}g" \
      CreateSequenceDictionary \
      -R "~{reference_fasta}" \
      -O "~{dict_basename}.dict" \
      --VERBOSITY WARNING
  >>>

  output {
    File sequence_dict = "~{dict_basename}.dict"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task collect_wgs_metrics {
  meta {
    description: "Collect whole genome sequencing metrics using GATK CollectWgsMetrics"
    outputs: {
        metrics_file: "Comprehensive WGS metrics file with coverage and quality statistics"
    }
  }

  parameter_meta {
    bam: "Input aligned BAM file"
    bam_index: "Index file for the input BAM"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    intervals: "Optional interval list file defining target regions"
    output_basename: "Base name for output files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
    minimum_mapping_quality: "Minimum mapping quality for reads to be included"
    minimum_base_quality: "Minimum base quality for bases to be included"
    coverage_cap: "Maximum coverage depth to analyze"
  }

  input {
    File bam
    File bam_index
    File reference_fasta
    File reference_fasta_index
    File? intervals
    String output_basename
    Int memory_gb = 8
    Int cpu_cores = 2
    Int minimum_mapping_quality = 20
    Int minimum_base_quality = 20
    Int coverage_cap = 250
  }

  command <<<
    set -eo pipefail
    
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      CollectWgsMetrics \
      -I "~{bam}" \
      -O "~{output_basename}.wgs_metrics.txt" \
      -R "~{reference_fasta}" \
      --MINIMUM_MAPPING_QUALITY ~{minimum_mapping_quality} \
      --MINIMUM_BASE_QUALITY ~{minimum_base_quality} \
      --COVERAGE_CAP ~{coverage_cap} \
      ~{if defined(intervals) then "--INTERVALS " + intervals else ""} \
      --VERBOSITY WARNING
  >>>

  output {
    File metrics_file = "~{output_basename}.wgs_metrics.txt"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task validate_outputs {
  meta {
    description: "Validate GATK outputs and generate comprehensive statistics report"
    outputs: {
        report: "Validation summary with file checks and basic statistics"
    }
  }

  parameter_meta {
    recalibrated_bams: "Array of recalibrated BAM files"
    recalibrated_bais: "Array of BAM index files"
    haplotype_vcfs: "Array of HaplotypeCaller VCF files"
    mutect2_vcfs: "Array of Mutect2 VCF files"
    wgs_metrics: "Array of WGS metrics files"
  }

  input {
    Array[File] recalibrated_bams
    Array[File] recalibrated_bais
    Array[File] haplotype_vcfs
    Array[File] mutect2_vcfs
    Array[File] wgs_metrics
  }

  command <<<
    set -eo pipefail
    
    echo "GATK Pipeline Validation Report" > validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "======================================" >> validation_report.txt
    echo "" >> validation_report.txt
    
    echo "Sample Summary:" >> validation_report.txt
    echo "Number of samples processed: ~{length(recalibrated_bams)}" >> validation_report.txt
    echo "" >> validation_report.txt
    
    echo "File Validation:" >> validation_report.txt
    echo "Checking file existence and basic properties..." >> validation_report.txt
    
    # Check BAM files
    for bam in ~{sep=" " recalibrated_bams}; do
      if [[ -f "$bam" ]]; then
        size=$(stat -f%z "$bam" 2>/dev/null || stat -c%s "$bam" 2>/dev/null || echo "unknown")
        echo "BAM: $(basename $bam) (${size} bytes)" >> validation_report.txt
      else
        echo "Missing BAM: $bam" >> validation_report.txt
      fi
    done

    # Check BAI files
    for bai in ~{sep=" " recalibrated_bais}; do
      if [[ -f "$bai" ]]; then
        size=$(stat -f%z "$bai" 2>/dev/null || stat -c%s "$bai" 2>/dev/null || echo "unknown")
        echo "BAM: $(basename $bai) (${size} bytes)" >> validation_report.txt
      else
        echo "Missing BAI: $bai" >> validation_report.txt
      fi
    done
    
    # Check VCF files
    for vcf in ~{sep=" " haplotype_vcfs}; do
      if [[ -f "$vcf" ]]; then
        variants=$(zcat "$vcf" | grep -v '^#' | wc -l || echo "0")
        echo "HaplotypeCaller VCF: $(basename $vcf) (${variants} variants)" >> validation_report.txt
      else
        echo "Missing HaplotypeCaller VCF: $vcf" >> validation_report.txt
      fi
    done
    
    # Check Mutect2 VCFs if present
    for vcf in ~{sep=" " mutect2_vcfs}; do
      if [[ -f "$vcf" ]]; then
        variants=$(zcat "$vcf" | grep -v "^#" | wc -l || echo "0")
        echo "Mutect2 VCF: $(basename $vcf) (${variants} variants)" >> validation_report.txt
      else
        echo "Missing Mutect2 VCF: $vcf" >> validation_report.txt
      fi
    done
    
    # Check metrics files
    for metrics in ~{sep=" " wgs_metrics}; do
      if [[ -f "$metrics" ]]; then
        echo "WGS Metrics: $(basename $metrics)" >> validation_report.txt
      else
        echo "Missing WGS Metrics: $metrics" >> validation_report.txt
      fi
    done
    
    echo "" >> validation_report.txt
    echo "Validation completed successfully." >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "4 GB"
    cpu: 1
  }
}
