version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gatk/ww-gatk.wdl" as ww_gatk
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bwa/ww-bwa.wdl" as ww_bwa

struct GatkSample {
    String name
    File bam
    File bai
}

workflow gatk_example {
  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_bam_data { }
  call ww_testdata.download_fastq_data { }

  # Create reference dictionary
  call ww_gatk.create_sequence_dictionary { input:
      reference_fasta = download_ref_data.fasta
  }

  # Download dbSNP VCF
  call ww_testdata.download_dbsnp_vcf { input:
    region = "NC_000001.11:1-10000000",
    filter_name = "chr1"
  }

  # Download known indels VCF
  call ww_testdata.download_known_indels_vcf { input:
    region = "chr1:1-10000000",
    filter_name = "chr1"
  }

  # Download gnomAD VCF
  call ww_testdata.download_gnomad_vcf { input:
    region = "chr1:1-10000000",
    filter_name = "chr1"
  }

  # Create test samples array
  Array[GatkSample] final_samples = [
    {
      "name": "demo_sample",
      "bam": download_bam_data.bam,
      "bai": download_bam_data.bai
    }
  ]

  # Test fastq_to_bam utility task
  call ww_gatk.fastq_to_bam { input:
      r1_fastq = [download_fastq_data.r1_fastq],
      r2_fastq = [download_fastq_data.r2_fastq],
      base_file_name = "test_fastq_to_bam",
      sample_name = "NA12878_test"
  }

  # Test validate_sam_file utility task on the unmapped BAM
  call ww_gatk.validate_sam_file as validate_unmapped_bam { input:
      input_file = fastq_to_bam.unmapped_bam,
      base_file_name = "test_fastq_to_bam_validation"
  }

  # Splitting intervals for parallel processing
  call ww_gatk.split_intervals { input:
      reference_fasta = download_ref_data.fasta,
      reference_fasta_index = download_ref_data.fasta_index,
      reference_dict = create_sequence_dictionary.sequence_dict,
      scatter_count = 2
  }

  # Scattering across the samples provided
  scatter (sample in final_samples) {
    # Mark duplicates in the raw BAM file
    call ww_gatk.mark_duplicates { input:
        bam = sample.bam,
        bam_index = sample.bai,
        base_file_name = sample.name
    }

    # Base recalibration using GATK BaseRecalibrator
    call ww_gatk.base_recalibrator { input:
        bam = mark_duplicates.markdup_bam,
        bam_index = mark_duplicates.markdup_bai,
        dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf,
        reference_fasta = download_ref_data.fasta,
        reference_fasta_index = download_ref_data.fasta_index,
        reference_dict = create_sequence_dictionary.sequence_dict,
        known_indels_sites_vcfs = [download_known_indels_vcf.known_indels_vcf],
        base_file_name = sample.name
    }

    # Collect WGS metrics after base recalibration
    call ww_gatk.collect_wgs_metrics { input:
        bam = base_recalibrator.recalibrated_bam,
        bam_index = base_recalibrator.recalibrated_bai,
        reference_fasta = download_ref_data.fasta,
        reference_fasta_index = download_ref_data.fasta_index,
        base_file_name = sample.name
    }

    # Perform all three operations in a single task
    call ww_gatk.markdup_recal_metrics { input:
        bam = sample.bam,
        bam_index = sample.bai,
        dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf,
        reference_fasta = download_ref_data.fasta,
        reference_fasta_index = download_ref_data.fasta_index,
        reference_dict = create_sequence_dictionary.sequence_dict,
        known_indels_sites_vcfs = [download_known_indels_vcf.known_indels_vcf],
        base_file_name = sample.name + ".combined",
        minimum_mapping_quality = 20,
        minimum_base_quality = 20,
        coverage_cap = 250
    }

    # Split BAM by intervals for scatter-gather approach
    call ww_gatk.print_reads { input:
      bam = base_recalibrator.recalibrated_bam,
      bam_index = base_recalibrator.recalibrated_bai,
      intervals = split_intervals.interval_files,
      reference_fasta = download_ref_data.fasta,
      reference_fasta_index = download_ref_data.fasta_index,
      reference_dict = create_sequence_dictionary.sequence_dict,
      output_basename = sample.name + ".recalibrated"
    }

    # Scatter HaplotypeCaller and Mutect2 across interval-specific BAMs
    scatter (i in range(length(split_intervals.interval_files))) {
      call ww_gatk.haplotype_caller { input:
        bam = print_reads.interval_bams[i],
        bam_index = print_reads.interval_bam_indices[i],
        intervals = split_intervals.interval_files[i],
        reference_fasta = download_ref_data.fasta,
        reference_fasta_index = download_ref_data.fasta_index,
        reference_dict = create_sequence_dictionary.sequence_dict,
        dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf,
        base_file_name = sample.name + "." + basename(split_intervals.interval_files[i], ".interval_list")
      }

      call ww_gatk.mutect2 { input:
        bam = print_reads.interval_bams[i],
        bam_index = print_reads.interval_bam_indices[i],
        intervals = split_intervals.interval_files[i],
        reference_fasta = download_ref_data.fasta,
        reference_fasta_index = download_ref_data.fasta_index,
        reference_dict = create_sequence_dictionary.sequence_dict,
        gnomad_vcf = download_gnomad_vcf.gnomad_vcf,
        base_file_name = sample.name + "." + basename(split_intervals.interval_files[i], ".interval_list")
      }
    }

    # Merge HaplotypeCaller results
    call ww_gatk.merge_vcfs as merge_haplotype_vcfs { input:
        vcfs = haplotype_caller.vcf,
        vcf_indices = haplotype_caller.vcf_index,
        base_file_name = sample.name + ".haplotypecaller",
        reference_dict = create_sequence_dictionary.sequence_dict
    }

    # Merge Mutect2 VCFs
    call ww_gatk.merge_vcfs as merge_mutect2_vcfs { input:
        vcfs = mutect2.vcf,
        vcf_indices = mutect2.vcf_index,
        base_file_name = sample.name + ".mutect2",
        reference_dict = create_sequence_dictionary.sequence_dict
    }

    # Merge Mutect2 stats
    call ww_gatk.merge_mutect_stats { input:
        stats = mutect2.stats_file,
        base_file_name = sample.name + ".mutect2"
    }

    # Run HaplotypeCaller with internal parallelization
    call ww_gatk.haplotype_caller_parallel { input:
        bam = base_recalibrator.recalibrated_bam,
        bam_index = base_recalibrator.recalibrated_bai,
        intervals = split_intervals.interval_files,
        reference_fasta = download_ref_data.fasta,
        reference_fasta_index = download_ref_data.fasta_index,
        reference_dict = create_sequence_dictionary.sequence_dict,
        dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf,
        base_file_name = sample.name
    }

    # Run Mutect2 with internal parallelization
    call ww_gatk.mutect2_parallel { input:
        bam = base_recalibrator.recalibrated_bam,
        bam_index = base_recalibrator.recalibrated_bai,
        intervals = split_intervals.interval_files,
        reference_fasta = download_ref_data.fasta,
        reference_fasta_index = download_ref_data.fasta_index,
        reference_dict = create_sequence_dictionary.sequence_dict,
        gnomad_vcf = download_gnomad_vcf.gnomad_vcf,
        base_file_name = sample.name
    }

    # Test validate_sam_file on recalibrated BAM
    call ww_gatk.validate_sam_file as validate_recalibrated_bam { input:
        input_file = base_recalibrator.recalibrated_bam,
        base_file_name = sample.name + ".recal_validation"
    }
  }

  # Test saturation mutagenesis analysis
  # Download a small region of chr1 (5KB to fit in CI/CD memory limits)
  call ww_testdata.download_ref_data as download_saturation_ref {
    input:
      chromo = "chr1",
      version = "hg38",
      region = "1000000-1005000"
  }

  # Create clean amplicon reference (removes N bases required by GATK AnalyzeSaturationMutagenesis)
  call ww_testdata.create_clean_amplicon_reference {
    input:
      input_fasta = download_saturation_ref.fasta,
      output_name = "chr1_saturation_amplicon",
      replace_n_with = "A"
  }

  # Build BWA index for the clean amplicon reference
  call ww_bwa.bwa_index as bwa_index_saturation {
    input:
      reference_fasta = create_clean_amplicon_reference.clean_fasta,
      cpu_cores = 2,
      memory_gb = 8
  }

  # Align test FASTQ to the clean amplicon reference to create a BAM with matching headers
  call ww_bwa.bwa_mem as bwa_mem_saturation {
    input:
      bwa_genome_tar = bwa_index_saturation.bwa_index_tar,
      reference_fasta = create_clean_amplicon_reference.clean_fasta,
      reads = download_fastq_data.r1_fastq,
      mates = download_fastq_data.r2_fastq,
      name = "saturation_test_sample",
      cpu_cores = 2,
      memory_gb = 8
  }

  # Analyze saturation mutagenesis with properly aligned BAM
  call ww_gatk.analyze_saturation_mutagenesis {
    input:
      bam = bwa_mem_saturation.sorted_bam,
      bam_index = bwa_mem_saturation.sorted_bai,
      reference_fasta = create_clean_amplicon_reference.clean_fasta,
      reference_fasta_index = create_clean_amplicon_reference.clean_fasta_index,
      reference_dict = create_clean_amplicon_reference.clean_dict,
      orf_range = "1-99",
      base_file_name = "saturation_test",
      memory_gb = 8,
      cpu_cores = 2
  }

  # Validate outputs to ensure all tasks completed successfully
  call validate_outputs { input:
      markdup_bams = mark_duplicates.markdup_bam,
      markdup_bais = mark_duplicates.markdup_bai,
      recalibrated_bams = base_recalibrator.recalibrated_bam,
      recalibrated_bais = base_recalibrator.recalibrated_bai,
      sequential_bams = markdup_recal_metrics.recalibrated_bam,
      sequential_bais = markdup_recal_metrics.recalibrated_bai,
      haplotype_vcfs = merge_haplotype_vcfs.merged_vcf,
      mutect2_vcfs = merge_mutect2_vcfs.merged_vcf,
      mutect2_stats = merge_mutect_stats.merged_stats,
      parallel_haplotype_vcfs = haplotype_caller_parallel.vcf,
      parallel_mutect2_vcfs = mutect2_parallel.vcf,
      wgs_metrics = collect_wgs_metrics.metrics_file,
      unmapped_bam = fastq_to_bam.unmapped_bam,
      unmapped_bam_validation = validate_unmapped_bam.validation_report,
      recalibrated_bam_validations = validate_recalibrated_bam.validation_report,
      saturation_variant_counts = analyze_saturation_mutagenesis.variant_counts,
      saturation_aa_counts = analyze_saturation_mutagenesis.aa_counts
  }

  output {
    Array[File] recalibrated_bams = base_recalibrator.recalibrated_bam
    Array[File] recalibrated_bais = base_recalibrator.recalibrated_bai
    Array[File] haplotype_vcfs = merge_haplotype_vcfs.merged_vcf
    Array[File] mutect2_vcfs = merge_mutect2_vcfs.merged_vcf
    Array[File] mutect2_stats = merge_mutect_stats.merged_stats
    Array[File] wgs_metrics = collect_wgs_metrics.metrics_file
    File saturation_variant_counts = analyze_saturation_mutagenesis.variant_counts
    File saturation_aa_counts = analyze_saturation_mutagenesis.aa_counts
    File validation_report = validate_outputs.report
  }
}

task validate_outputs {
  # TODO: Add validation for GATK dictionary, intervals, and duplicate metrics files
  # TODO: Ensure that the recalibrated BAM files match between the separate and all-in-one tasks
  meta {
    description: "Validate GATK outputs and generate comprehensive statistics report"
    outputs: {
        report: "Validation summary with file checks and basic statistics"
    }
  }

  parameter_meta {
    markdup_bams: "Array of MarkDuplicates BAM files"
    markdup_bais: "Array of MarkDuplicates BAM index files"
    recalibrated_bams: "Array of recalibrated BAM files"
    recalibrated_bais: "Array of recalibrated BAM index files"
    sequential_bams: "Array of sequential Markdup-Recal-Metrics BAM files"
    sequential_bais: "Array of sequential Markdup-Recal-Metrics BAM index files"
    haplotype_vcfs: "Array of HaplotypeCaller VCF files called via scatter-gather parallelization"
    mutect2_vcfs: "Array of Mutect2 VCF files called via scatter-gather parallelization"
    mutect2_stats: "Array of merged Mutect2 statistics files"
    parallel_haplotype_vcfs: "Array of HaplotypeCaller VCF files called via internal parallelization"
    parallel_mutect2_vcfs: "Array of Mutect2 VCF files called via internal parallelization"
    wgs_metrics: "Array of WGS metrics files"
    unmapped_bam: "Unmapped BAM file created from FASTQ using fastq_to_bam"
    unmapped_bam_validation: "Validation report for the unmapped BAM"
    recalibrated_bam_validations: "Array of validation reports for recalibrated BAMs"
    saturation_variant_counts: "Variant count table from saturation mutagenesis analysis"
    saturation_aa_counts: "Amino acid count table from saturation mutagenesis analysis"
  }

  input {
    Array[File] markdup_bams
    Array[File] markdup_bais
    Array[File] recalibrated_bams
    Array[File] recalibrated_bais
    Array[File] sequential_bams
    Array[File] sequential_bais
    Array[File] haplotype_vcfs
    Array[File] mutect2_vcfs
    Array[File] mutect2_stats
    Array[File] parallel_haplotype_vcfs
    Array[File] parallel_mutect2_vcfs
    Array[File] wgs_metrics
    File unmapped_bam
    File unmapped_bam_validation
    Array[File] recalibrated_bam_validations
    File saturation_variant_counts
    File saturation_aa_counts
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

    # Check MarkDuplicates BAM files
    for bam in ~{sep=" " markdup_bams}; do
      if [[ -f "$bam" ]]; then
        size=$(stat -f%z "$bam" 2>/dev/null || stat -c%s "$bam" 2>/dev/null || echo "unknown")
        echo "MarkDuplicates BAM: $(basename $bam) (${size} bytes)" >> validation_report.txt
      else
        echo "Missing MarkDuplicates BAM: $bam" >> validation_report.txt
      fi
    done

    # Check MarkDuplicates BAM index files
    for bai in ~{sep=" " markdup_bais}; do
      if [[ -f "$bai" ]]; then
        size=$(stat -f%z "$bai" 2>/dev/null || stat -c%s "$bai" 2>/dev/null || echo "unknown")
        echo "MarkDuplicates BAM Index: $(basename $bai) (${size} bytes)" >> validation_report.txt
      else
        echo "Missing MarkDuplicates BAM Index: $bai" >> validation_report.txt
      fi
    done

    # Check Recalibrated BAM files
    for bam in ~{sep=" " recalibrated_bams}; do
      if [[ -f "$bam" ]]; then
        size=$(stat -f%z "$bam" 2>/dev/null || stat -c%s "$bam" 2>/dev/null || echo "unknown")
        echo "Recalibrated BAM: $(basename $bam) (${size} bytes)" >> validation_report.txt
      else
        echo "Missing Recalibrated BAM: $bam" >> validation_report.txt
      fi
    done

    # Check Sequential Markdup-Recal-Metrics BAM files
    for bam in ~{sep=" " sequential_bams}; do
      if [[ -f "$bam" ]]; then
        size=$(stat -f%z "$bam" 2>/dev/null || stat -c%s "$bam" 2>/dev/null || echo "unknown")
        echo "Sequential BAM: $(basename $bam) (${size} bytes)" >> validation_report.txt
      else
        echo "Missing Sequential BAM: $bam" >> validation_report.txt
      fi
    done

    # Check Sequential Markdup-Recal-Metrics BAI files
    for bai in ~{sep=" " sequential_bais}; do
      if [[ -f "$bai" ]]; then
        size=$(stat -f%z "$bai" 2>/dev/null || stat -c%s "$bai" 2>/dev/null || echo "unknown")
        echo "Sequential BAI: $(basename $bai) (${size} bytes)" >> validation_report.txt
      else
        echo "Missing Sequential BAI: $bai" >> validation_report.txt
      fi
    done

    # Check Recalibrated BAM index files
    for bai in ~{sep=" " recalibrated_bais}; do
      if [[ -f "$bai" ]]; then
        size=$(stat -f%z "$bai" 2>/dev/null || stat -c%s "$bai" 2>/dev/null || echo "unknown")
        echo "Recalibrated BAM Index: $(basename $bai) (${size} bytes)" >> validation_report.txt
      else
        echo "Missing Recal BAM Index: $bai" >> validation_report.txt
      fi
    done

    # Check HaplotypeCaller VCFs if present
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

    # Check Parallel HaplotypeCaller VCFs if present
    for vcf in ~{sep=" " parallel_haplotype_vcfs}; do
      if [[ -f "$vcf" ]]; then
        variants=$(zcat "$vcf" | grep -v '^#' | wc -l || echo "0")
        echo "Parallel HaplotypeCaller VCF: $(basename $vcf) (${variants} variants)" >> validation_report.txt
      else
        echo "Missing Parallel HaplotypeCaller VCF: $vcf" >> validation_report.txt
      fi
    done

    # Check Parallel Mutect2 VCFs if present
    for vcf in ~{sep=" " parallel_mutect2_vcfs}; do
      if [[ -f "$vcf" ]]; then
        variants=$(zcat "$vcf" | grep -v "^#" | wc -l || echo "0")
        echo "Parallel Mutect2 VCF: $(basename $vcf) (${variants} variants)" >> validation_report.txt
      else
        echo "Missing Parallel Mutect2 VCF: $vcf" >> validation_report.txt
      fi
    done

    # Check Mutect2 stats files
    for stats in ~{sep=" " mutect2_stats}; do
      if [[ -f "$stats" ]]; then
        size=$(stat -f%z "$stats" 2>/dev/null || stat -c%s "$stats" 2>/dev/null || echo "unknown")
        echo "Mutect2 Merged Stats: $(basename $stats) (${size} bytes)" >> validation_report.txt
      else
        echo "Missing Mutect2 Merged Stats: $stats" >> validation_report.txt
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

    # Check unmapped BAM from fastq_to_bam
    if [[ -f "~{unmapped_bam}" ]]; then
      size=$(stat -f%z "~{unmapped_bam}" 2>/dev/null || stat -c%s "~{unmapped_bam}" 2>/dev/null || echo "unknown")
      echo "Unmapped BAM from FASTQ: $(basename ~{unmapped_bam}) (${size} bytes)" >> validation_report.txt
    else
      echo "Missing Unmapped BAM from FASTQ: ~{unmapped_bam}" >> validation_report.txt
    fi

    # Check unmapped BAM validation report
    if [[ -f "~{unmapped_bam_validation}" ]]; then
      echo "Unmapped BAM Validation Report: $(basename ~{unmapped_bam_validation})" >> validation_report.txt
      echo "  Contents:" >> validation_report.txt
      cat "~{unmapped_bam_validation}" | head -10 | sed 's/^/    /' >> validation_report.txt
    else
      echo "Missing Unmapped BAM Validation Report: ~{unmapped_bam_validation}" >> validation_report.txt
    fi

    # Check recalibrated BAM validation reports
    for validation in ~{sep=" " recalibrated_bam_validations}; do
      if [[ -f "$validation" ]]; then
        echo "Recalibrated BAM Validation: $(basename $validation)" >> validation_report.txt
      else
        echo "Missing Recalibrated BAM Validation: $validation" >> validation_report.txt
      fi
    done

    # Check saturation mutagenesis variant counts
    if [[ -f "~{saturation_variant_counts}" ]]; then
      size=$(stat -f%z "~{saturation_variant_counts}" 2>/dev/null || stat -c%s "~{saturation_variant_counts}" 2>/dev/null || echo "unknown")
      lines=$(wc -l < "~{saturation_variant_counts}" 2>/dev/null || echo "unknown")
      echo "Saturation Variant Counts: $(basename ~{saturation_variant_counts}) (${size} bytes, ${lines} lines)" >> validation_report.txt
    else
      echo "Missing Saturation Variant Counts: ~{saturation_variant_counts}" >> validation_report.txt
    fi

    # Check saturation mutagenesis amino acid counts
    if [[ -f "~{saturation_aa_counts}" ]]; then
      size=$(stat -f%z "~{saturation_aa_counts}" 2>/dev/null || stat -c%s "~{saturation_aa_counts}" 2>/dev/null || echo "unknown")
      lines=$(wc -l < "~{saturation_aa_counts}" 2>/dev/null || echo "unknown")
      echo "Saturation AA Counts: $(basename ~{saturation_aa_counts}) (${size} bytes, ${lines} lines)" >> validation_report.txt
    else
      echo "Missing Saturation AA Counts: ~{saturation_aa_counts}" >> validation_report.txt
    fi

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
