version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-gatk/ww-gatk.wdl" as ww_gatk

struct GatkSample {
    String name
    File bam
    File bai
}

workflow gatk_example {
  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_bam_data { }

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
      parallel_haplotype_vcfs = haplotype_caller_parallel.vcf,
      parallel_mutect2_vcfs = mutect2_parallel.vcf,
      wgs_metrics = collect_wgs_metrics.metrics_file
  }

  output {
    Array[File] recalibrated_bams = base_recalibrator.recalibrated_bam
    Array[File] recalibrated_bais = base_recalibrator.recalibrated_bai
    Array[File] haplotype_vcfs = merge_haplotype_vcfs.merged_vcf
    Array[File] mutect2_vcfs = merge_mutect2_vcfs.merged_vcf
    Array[File] wgs_metrics = collect_wgs_metrics.metrics_file
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
    parallel_haplotype_vcfs: "Array of HaplotypeCaller VCF files called via internal parallelization"
    parallel_mutect2_vcfs: "Array of Mutect2 VCF files called via internal parallelization"
    wgs_metrics: "Array of WGS metrics files"
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
    Array[File] parallel_haplotype_vcfs
    Array[File] parallel_mutect2_vcfs
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
