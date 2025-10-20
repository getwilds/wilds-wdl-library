## WILDS WDL module for GATK variant calling and analysis tasks.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct GatkSample {
    String name
    File bam
    File bai
}

workflow gatk_example {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow demonstrating GATK variant calling tasks with test data"
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

  # Download test data
  call ww_testdata.download_ref_data { }
  call ww_testdata.download_bam_data { }

  # Create reference dictionary
  call create_sequence_dictionary { input:
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
  call split_intervals { input:
      reference_fasta = download_ref_data.fasta,
      reference_fasta_index = download_ref_data.fasta_index,
      reference_dict = create_sequence_dictionary.sequence_dict,
      scatter_count = 2
  }

  # Scattering across the samples provided
  scatter (sample in final_samples) {
    # Mark duplicates in the raw BAM file
    call mark_duplicates { input:
        bam = sample.bam,
        bam_index = sample.bai,
        base_file_name = sample.name
    }

    # Base recalibration using GATK BaseRecalibrator
    call base_recalibrator { input:
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
    call collect_wgs_metrics { input:
        bam = base_recalibrator.recalibrated_bam,
        bam_index = base_recalibrator.recalibrated_bai,
        reference_fasta = download_ref_data.fasta,
        reference_fasta_index = download_ref_data.fasta_index,
        base_file_name = sample.name
    }

    # Perform all three operations in a single task
    call markdup_recal_metrics { input:
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
    call print_reads { input:
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
      call haplotype_caller { input:
        bam = print_reads.interval_bams[i],
        bam_index = print_reads.interval_bam_indices[i],
        intervals = split_intervals.interval_files[i],
        reference_fasta = download_ref_data.fasta,
        reference_fasta_index = download_ref_data.fasta_index,
        reference_dict = create_sequence_dictionary.sequence_dict,
        dbsnp_vcf = download_dbsnp_vcf.dbsnp_vcf,
        base_file_name = sample.name + "." + basename(split_intervals.interval_files[i], ".interval_list")
      }

      call mutect2 { input:
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
    call merge_vcfs as merge_haplotype_vcfs { input:
        vcfs = haplotype_caller.vcf,
        vcf_indices = haplotype_caller.vcf_index,
        base_file_name = sample.name + ".haplotypecaller",
        reference_dict = create_sequence_dictionary.sequence_dict
    }

    # Merge Mutect2 VCFs
    call merge_vcfs as merge_mutect2_vcfs { input:
        vcfs = mutect2.vcf,
        vcf_indices = mutect2.vcf_index,
        base_file_name = sample.name + ".mutect2",
        reference_dict = create_sequence_dictionary.sequence_dict
    }

    # Merge Mutect2 stats
    call merge_mutect_stats { input:
        stats = mutect2.stats_file,
        base_file_name = sample.name + ".mutect2"
    }

    # Run HaplotypeCaller with internal parallelization
    call haplotype_caller_parallel { input:
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
    call mutect2_parallel { input:
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

  String dict_basename = basename(basename(reference_fasta, ".fa"), ".fasta")

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
    bam: "Aligned input BAM file"
    bam_index: "Index file for the aligned input BAM"
    base_file_name: "Base name for the output files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File bam
    File bam_index
    String base_file_name
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail

    gatk --java-options "-Dsamjdk.compression_level=5 -Xms16g -Xmx16g" \
      MarkDuplicates \
      --INPUT "~{bam}" \
      --OUTPUT "~{base_file_name}.markdup.bam" \
      --METRICS_FILE "~{base_file_name}.duplicate_metrics" \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --VERBOSITY WARNING
    samtools index "~{base_file_name}.markdup.bam"
  >>>

  output {
    File markdup_bam = "~{base_file_name}.markdup.bam"
    File markdup_bai = "~{base_file_name}.markdup.bam.bai"
    File duplicate_metrics = "~{base_file_name}.duplicate_metrics"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
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
    bam: "Input aligned BAM file to be recalibrated"
    bam_index: "Index file for the input BAM"
    dbsnp_vcf: "dbSNP VCF file for known variant sites"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    reference_dict: "Reference genome sequence dictionary"
    known_indels_sites_vcfs: "Array of VCF files with known indel sites"
    base_file_name: "Base name for output files"
    intervals: "Optional interval list file defining target regions"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File bam
    File bam_index
    File dbsnp_vcf
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] known_indels_sites_vcfs
    String base_file_name
    File? intervals
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
      -I "~{bam}" \
      -O "~{base_file_name}.recal_data.table" \
      --known-sites "~{dbsnp_vcf}" \
      --known-sites ~{sep=" --known-sites " known_indels_sites_vcfs} \
      ~{if defined(intervals) then "--intervals " + intervals else ""} \
      --verbosity WARNING

    # Apply Base Quality Score Recalibration
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      ApplyBQSR \
      -R "~{basename(reference_fasta)}" \
      -I "~{bam}" \
      -bqsr "~{base_file_name}.recal_data.table" \
      -O "~{base_file_name}.recal.bam" \
      ~{if defined(intervals) then "--intervals " + intervals else ""} \
      --verbosity WARNING

    # Index resulting bam file
    samtools index "~{base_file_name}.recal.bam"
  >>>

  output {
    File recalibrated_bam = "~{base_file_name}.recal.bam"
    File recalibrated_bai = "~{base_file_name}.recal.bam.bai"
    File recalibration_report = "~{base_file_name}.recal_data.table"
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
    base_file_name: "Base name for output files"
    intervals: "Optional interval list file defining target regions"
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
    String base_file_name
    File? intervals
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
      -O "~{base_file_name}.wgs_metrics.txt" \
      -R "~{reference_fasta}" \
      --MINIMUM_MAPPING_QUALITY ~{minimum_mapping_quality} \
      --MINIMUM_BASE_QUALITY ~{minimum_base_quality} \
      --COVERAGE_CAP ~{coverage_cap} \
      ~{if defined(intervals) then "--INTERVALS " + intervals else ""} \
      --VERBOSITY WARNING
  >>>

  output {
    File metrics_file = "~{base_file_name}.wgs_metrics.txt"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task markdup_recal_metrics {
  meta {
    description: "Performs duplicate marking, base recalibration, and WGS metrics in a single task to avoid data duplication"
    outputs: {
        recalibrated_bam: "BAM file with recalibrated base quality scores",
        recalibrated_bai: "Index file for the recalibrated BAM",
        recalibration_report: "Base recalibration report table",
        duplicate_metrics: "Metrics file containing duplicate marking statistics",
        wgs_metrics: "Comprehensive WGS metrics file with coverage and quality statistics"
    }
  }

  parameter_meta {
    bam: "Aligned input BAM file"
    bam_index: "Index file for the aligned input BAM"
    dbsnp_vcf: "dbSNP VCF file for known variant sites"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    reference_dict: "Reference genome sequence dictionary"
    known_indels_sites_vcfs: "Array of VCF files with known indel sites"
    base_file_name: "Base name for output files"
    intervals: "Optional interval list file defining target regions"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
    minimum_mapping_quality: "Minimum mapping quality for reads to be included"
    minimum_base_quality: "Minimum base quality for bases to be included"
    coverage_cap: "Maximum coverage depth to analyze"
  }

  input {
    File bam
    File bam_index
    File dbsnp_vcf
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] known_indels_sites_vcfs
    String base_file_name
    File? intervals
    Int memory_gb = 8
    Int cpu_cores = 2
    Int minimum_mapping_quality = 20
    Int minimum_base_quality = 20
    Int coverage_cap = 250
  }

  command <<<
    set -eo pipefail

    # Mark duplicates in the aligned BAM file
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms16g -Xmx16g" \
      MarkDuplicates \
      --INPUT "~{bam}" \
      --OUTPUT "~{base_file_name}.markdup.bam" \
      --METRICS_FILE "~{base_file_name}.duplicate_metrics" \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --VERBOSITY WARNING

    # Index resulting bam file
    samtools index "~{base_file_name}.markdup.bam"

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
      -I "~{bam}" \
      -O "~{base_file_name}.recal_data.table" \
      --known-sites "~{dbsnp_vcf}" \
      --known-sites ~{sep=" --known-sites " known_indels_sites_vcfs} \
      ~{if defined(intervals) then "--intervals " + intervals else ""} \
      --verbosity WARNING

    # Apply Base Quality Score Recalibration
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      ApplyBQSR \
      -R "~{basename(reference_fasta)}" \
      -I "~{bam}" \
      -bqsr "~{base_file_name}.recal_data.table" \
      -O "~{base_file_name}.recal.bam" \
      ~{if defined(intervals) then "--intervals " + intervals else ""} \
      --verbosity WARNING

    # Index resulting bam file
    samtools index "~{base_file_name}.recal.bam"

    # Cleaning up intermediate MarkDuplicates BAM file
    rm "~{base_file_name}.markdup.bam" "~{base_file_name}.markdup.bam.bai"

    # Collect WGS metrics
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      CollectWgsMetrics \
      -I "~{bam}" \
      -O "~{base_file_name}.wgs_metrics.txt" \
      -R "~{basename(reference_fasta)}" \
      --MINIMUM_MAPPING_QUALITY ~{minimum_mapping_quality} \
      --MINIMUM_BASE_QUALITY ~{minimum_base_quality} \
      --COVERAGE_CAP ~{coverage_cap} \
      ~{if defined(intervals) then "--INTERVALS " + intervals else ""} \
      --VERBOSITY WARNING
  >>>

  output {
    File recalibrated_bam = "~{base_file_name}.recal.bam"
    File recalibrated_bai = "~{base_file_name}.recal.bam.bai"
    File recalibration_report = "~{base_file_name}.recal_data.table"
    File duplicate_metrics = "~{base_file_name}.duplicate_metrics"
    File wgs_metrics = "~{base_file_name}.wgs_metrics.txt"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task split_intervals {
  meta {
    description: "Split intervals into smaller chunks for parallelization using GATK SplitIntervals"
    outputs: {
        interval_files: "Array of interval files optimized for parallel processing"
    }
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Reference genome FASTA index file"
    reference_dict: "Reference genome sequence dictionary"
    intervals: "Optional interval list file defining target regions to split"
    scatter_count: "Number of interval files to create (default: 24)"
    filter_to_canonical_chromosomes: "Whether to restrict analysis to canonical chromosomes (chr1-22,X,Y,M) (default: true)"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File? intervals
    Int scatter_count = 24
    Boolean filter_to_canonical_chromosomes = true
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail

    # Add local symbolic link for reference files
    ln -s "~{reference_fasta}" "~{basename(reference_fasta)}"
    ln -s "~{reference_fasta_index}" "~{basename(reference_fasta_index)}"
    ln -s "~{reference_dict}" "~{basename(reference_dict)}"

    # Create output directories
    mkdir -p scattered_intervals

    # Create canonical chromosome intervals if filtering is enabled and no custom intervals provided
    if [[ "~{filter_to_canonical_chromosomes}" == "true" && ! -f "~{intervals}" ]]; then
      echo "Creating canonical chromosome intervals..."
      echo "@HD	VN:1.0	SO:coordinate" > canonical_chromosomes.interval_list

      # Add sequence dictionary headers for canonical chromosomes only
      grep -E "^@SQ.*SN:(chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY|chrM)\s" "~{basename(reference_dict)}" >> canonical_chromosomes.interval_list || true

      # Add interval entries for canonical chromosomes
      grep -E "^@SQ.*SN:(chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY|chrM)\s" "~{basename(reference_dict)}" | \
        awk '{
          # Extract sequence name (SN:)
          for(i=1; i<=NF; i++) {
            if($i ~ /^SN:/) {
              gsub(/^SN:/, "", $i);
              seq_name = $i;
            }
            if($i ~ /^LN:/) {
              gsub(/^LN:/, "", $i);
              seq_length = $i;
            }
          }
          if(seq_name && seq_length) {
            print seq_name "\t1\t" seq_length "\t+\t.";
          }
          seq_name=""; seq_length="";
        }' >> canonical_chromosomes.interval_list

      INTERVAL_ARG="--intervals canonical_chromosomes.interval_list"
    elif [[ -f "~{intervals}" ]]; then
      INTERVAL_ARG="--intervals ~{intervals}"
    else
      INTERVAL_ARG=""
    fi

    # Run SplitIntervals
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      SplitIntervals \
      -R "~{basename(reference_fasta)}" \
      --scatter-count ~{scatter_count} \
      ${INTERVAL_ARG} \
      -O scattered_intervals/ \
      --verbosity WARNING

    # List all created files for output
    find scattered_intervals/ -name "*.interval_list" | sort -V > interval_files.txt
  >>>

  output {
    Array[File] interval_files = read_lines("interval_files.txt")
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task print_reads {
  meta {
    description: "Extract reads from specific intervals using GATK PrintReads"
    outputs: {
        interval_bams: "Array of BAM files containing reads from specified intervals",
        interval_bam_indices: "Array of index files for the interval BAMs"
    }
  }

  parameter_meta {
    bam: "Input BAM file to extract reads from"
    bam_index: "Index file for the input BAM"
    intervals: "Array of interval files defining regions to extract"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    reference_dict: "Reference genome sequence dictionary"
    output_basename: "Base name for output files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File bam
    File bam_index
    Array[File] intervals
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    String output_basename
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail

    # Add local symbolic link for reference files
    ln -s "~{reference_fasta}" "~{basename(reference_fasta)}"
    ln -s "~{reference_fasta_index}" "~{basename(reference_fasta_index)}"
    ln -s "~{reference_dict}" "~{basename(reference_dict)}"

    # Create arrays to store output filenames
    bam_files=()
    bai_files=()

    # Process each interval file
    counter=0
    for interval_file in ~{sep=" " intervals}; do
      interval_name=$(basename "$interval_file" .interval_list)
      output_bam="~{output_basename}.${counter}.${interval_name}.bam"

      gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
        PrintReads \
        -R "~{basename(reference_fasta)}" \
        -I "~{bam}" \
        -L "$interval_file" \
        -O "$output_bam" \
        --verbosity WARNING
      samtools index "$output_bam"

      bam_files+=("$output_bam")
      bai_files+=("${output_bam}.bai")
      counter=$((counter + 1))
    done

    # Write output file lists
    printf '%s\n' "${bam_files[@]}" > bam_files.txt
    printf '%s\n' "${bai_files[@]}" > bai_files.txt
  >>>

  output {
    Array[File] interval_bams = read_lines("bam_files.txt")
    Array[File] interval_bam_indices = read_lines("bai_files.txt")
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
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    reference_dict: "Reference genome sequence dictionary"
    dbsnp_vcf: "dbSNP VCF file for variant annotation"
    base_file_name: "Base name for output files"
    intervals: "Optional interval list file defining target regions"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File bam
    File bam_index
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File dbsnp_vcf
    String base_file_name
    File? intervals
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

    # Create index for dbSNP vcf
    gatk IndexFeatureFile -I "~{dbsnp_vcf}"

    # Run HaplotypeCaller
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      HaplotypeCaller \
      -R "~{basename(reference_fasta)}" \
      -I "~{bam}" \
      -O "~{base_file_name}.haplotypecaller.vcf.gz" \
      --dbsnp "~{dbsnp_vcf}" \
      ~{if defined(intervals) then "--intervals " + intervals else ""} \
      ~{if defined(intervals) then "--interval-padding 100" else ""} \
      --verbosity WARNING
  >>>

  output {
    File vcf = "~{base_file_name}.haplotypecaller.vcf.gz"
    File vcf_index = "~{base_file_name}.haplotypecaller.vcf.gz.tbi"
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
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    reference_dict: "Reference genome sequence dictionary"
    gnomad_vcf: "gnomAD population allele frequency VCF for germline resource"
    base_file_name: "Base name for output files"
    intervals: "Optional interval list file defining target regions"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File bam
    File bam_index
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File gnomad_vcf
    String base_file_name
    File? intervals
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
      -O "~{base_file_name}.unfiltered.vcf.gz" \
      ~{if defined(intervals) then "--intervals " + intervals else ""} \
      ~{if defined(intervals) then "--interval-padding 100" else ""} \
      --germline-resource "~{gnomad_vcf}" \
      --f1r2-tar-gz "~{base_file_name}.f1r2.tar.gz" \
      --verbosity WARNING

    # Filter Mutect2 calls
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      FilterMutectCalls \
      -V "~{base_file_name}.unfiltered.vcf.gz" \
      -R "~{basename(reference_fasta)}" \
      -O "~{base_file_name}.mutect2.vcf.gz" \
      --stats "~{base_file_name}.unfiltered.vcf.gz.stats" \
      --verbosity WARNING
  >>>

  output {
    File vcf = "~{base_file_name}.mutect2.vcf.gz"
    File vcf_index = "~{base_file_name}.mutect2.vcf.gz.tbi"
    File stats_file = "~{base_file_name}.unfiltered.vcf.gz.stats"
    File f1r2_counts = "~{base_file_name}.f1r2.tar.gz"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task merge_vcfs {
  meta {
    description: "Merge multiple VCF files into a single VCF"
    outputs: {
        merged_vcf: "Merged VCF file",
        merged_vcf_index: "Index for merged VCF file"
    }
  }

  parameter_meta {
    vcfs: "Array of VCF files to merge"
    vcf_indices: "Array of VCF index files"
    base_file_name: "Base name for output files"
    reference_dict: "Reference sequence dictionary"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices
    String base_file_name
    File reference_dict
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail

    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      MergeVcfs \
      -I ~{sep=" -I " vcfs} \
      -D "~{reference_dict}" \
      -O "~{base_file_name}.merged.vcf.gz" \
      --VERBOSITY WARNING
  >>>

  output {
    File merged_vcf = "~{base_file_name}.merged.vcf.gz"
    File merged_vcf_index = "~{base_file_name}.merged.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task merge_mutect_stats {
  meta {
    description: "Merge Mutect2 statistics files"
    outputs: {
        merged_stats: "Merged Mutect2 statistics file"
    }
  }

  parameter_meta {
    stats: "Array of Mutect2 stats files to merge"
    base_file_name: "Base name for output files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    Array[File] stats
    String base_file_name
    Int memory_gb = 4
    Int cpu_cores = 1
  }

  command <<<
    set -eo pipefail

    gatk --java-options "-Xms~{memory_gb - 2}g -Xmx~{memory_gb - 1}g" \
      MergeMutectStats \
      --stats ~{sep=" --stats " stats} \
      -O "~{base_file_name}.merged.stats" \
      --verbosity WARNING
  >>>

  output {
    File merged_stats = "~{base_file_name}.merged.stats"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task haplotype_caller_parallel {
  meta {
    description: "Call germline variants using GATK HaplotypeCaller with internal parallelization for reduced data duplication"
    outputs: {
        vcf: "Compressed VCF file containing germline variant calls",
        vcf_index: "Index file for the VCF output"
    }
  }

  parameter_meta {
    bam: "Input aligned BAM file"
    bam_index: "Index file for the input BAM"
    intervals: "Array of interval files for parallel processing"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    reference_dict: "Reference genome sequence dictionary"
    dbsnp_vcf: "dbSNP VCF file for variant annotation"
    base_file_name: "Base name for output files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use (should match number of intervals)"
  }

  input {
    File bam
    File bam_index
    Array[File] intervals
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File dbsnp_vcf
    String base_file_name
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail

    # Add local symbolic link for reference fasta and dict
    ln -s "~{reference_fasta}" "~{basename(reference_fasta)}"
    ln -s "~{reference_fasta_index}" "~{basename(reference_fasta_index)}"
    ln -s "~{reference_dict}" "~{basename(reference_dict)}"

    # Create index for dbSNP vcf
    gatk IndexFeatureFile -I "~{dbsnp_vcf}"

    # Calculate memory per parallel job
    mem_per_job=$(( (~{memory_gb} - 4) / ~{length(intervals)} ))
    if [ $mem_per_job -lt 2 ]; then
      mem_per_job=2
    fi

    # Create array of interval files
    intervals=(~{sep=" " intervals})

    # Function to run HaplotypeCaller on a single interval
    run_haplotypecaller() {
      local interval_file=$1
      local interval_name=$(basename "$interval_file" .interval_list)

      gatk --java-options "-Xms${mem_per_job}g -Xmx${mem_per_job}g" \
        HaplotypeCaller \
        -R "~{basename(reference_fasta)}" \
        -I "~{bam}" \
        -O "~{base_file_name}.${interval_name}.vcf.gz" \
        --intervals "$interval_file" \
        --interval-padding 100 \
        --dbsnp "~{dbsnp_vcf}" \
        --verbosity WARNING
    }

    # Export the function so it can be used by parallel
    export -f run_haplotypecaller
    export reference_fasta="~{basename(reference_fasta)}"
    export bam="~{bam}"
    export dbsnp_vcf="~{dbsnp_vcf}"
    export base_file_name="~{base_file_name}"
    export mem_per_job

    # Run HaplotypeCaller in parallel across intervals
    printf '%s\n' "${intervals[@]}" | \
      parallel -j ~{cpu_cores} run_haplotypecaller {}

    # Collect all VCF files for merging
    find . -name "~{base_file_name}.*.vcf.gz" | sort -V > vcf_list.txt

    # Merge VCFs using MergeVcfs with input list file
    gatk --java-options "-Xms4g -Xmx6g" \
      MergeVcfs \
      --INPUT vcf_list.txt \
      -D "~{basename(reference_dict)}" \
      -O "~{base_file_name}.haplotypecaller.vcf.gz" \
      --VERBOSITY WARNING
  >>>

  output {
    File vcf = "~{base_file_name}.haplotypecaller.vcf.gz"
    File vcf_index = "~{base_file_name}.haplotypecaller.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task mutect2_parallel {
  meta {
    description: "Call somatic variants using GATK Mutect2 with internal parallelization for reduced data duplication"
    outputs: {
        vcf: "Compressed VCF file containing filtered somatic variant calls",
        vcf_index: "Index file for the Mutect2 VCF output",
        stats_file: "Merged Mutect2 statistics file"
    }
  }

  parameter_meta {
    bam: "Input aligned BAM file"
    bam_index: "Index file for the input BAM"
    intervals: "Array of interval files for parallel processing"
    reference_fasta: "Reference genome FASTA file"
    reference_fasta_index: "Index file for the reference FASTA"
    reference_dict: "Reference genome sequence dictionary"
    gnomad_vcf: "gnomAD population allele frequency VCF for germline resource"
    base_file_name: "Base name for output files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File bam
    File bam_index
    Array[File] intervals
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File gnomad_vcf
    String base_file_name
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

    # Calculate memory per parallel job
    mem_per_job=$(( (~{memory_gb} - 4) / ~{length(intervals)} ))
    if [ $mem_per_job -lt 2 ]; then
      mem_per_job=2
    fi

    # Create array of interval files
    intervals=(~{sep=" " intervals})

    # Function to run Mutect2 on a single interval
    run_mutect2() {
      local interval_file=$1
      local interval_name=$(basename "$interval_file" .interval_list)

      # Run Mutect2
      gatk --java-options "-Xms${mem_per_job}g -Xmx${mem_per_job}g" \
        Mutect2 \
        -R "~{basename(reference_fasta)}" \
        -I "~{bam}" \
        -O "~{base_file_name}.${interval_name}.unfiltered.vcf.gz" \
        --intervals "$interval_file" \
        --interval-padding 100 \
        --germline-resource "~{gnomad_vcf}" \
        --f1r2-tar-gz "~{base_file_name}.${interval_name}.f1r2.tar.gz" \
        --verbosity WARNING

      # Filter Mutect2 calls
      gatk --java-options "-Xms${mem_per_job}g -Xmx${mem_per_job}g" \
        FilterMutectCalls \
        -V "~{base_file_name}.${interval_name}.unfiltered.vcf.gz" \
        -R "~{basename(reference_fasta)}" \
        -O "~{base_file_name}.${interval_name}.vcf.gz" \
        --stats "~{base_file_name}.${interval_name}.unfiltered.vcf.gz.stats" \
        --verbosity WARNING
    }

    # Export the function and variables
    export -f run_mutect2
    export reference_fasta="~{basename(reference_fasta)}"
    export bam="~{bam}"
    export gnomad_vcf="~{gnomad_vcf}"
    export base_file_name="~{base_file_name}"
    export mem_per_job

    # Run Mutect2 in parallel across intervals
    printf '%s\n' "${intervals[@]}" | \
      parallel -j ~{cpu_cores} run_mutect2 {}

    # Collect VCF files and stats files for merging
    find . -name "~{base_file_name}.*.vcf.gz" -not -name "*.unfiltered.*" | sort -V > vcf_list.txt
    find . -name "~{base_file_name}.*.unfiltered.vcf.gz.stats" | sort -V > stats_list.txt

    # Merge VCFs using input list file
    gatk --java-options "-Xms4g -Xmx6g" \
      MergeVcfs \
      --INPUT vcf_list.txt \
      -D "~{basename(reference_dict)}" \
      -O "~{base_file_name}.mutect2.vcf.gz" \
      --VERBOSITY WARNING

    # Merge stats files - build argument list properly
    stats_args=""
    while IFS= read -r stats_file; do
      stats_args="${stats_args} --stats ${stats_file}"
    done < stats_list.txt

    gatk --java-options "-Xms2g -Xmx3g" \
      MergeMutectStats \
      ${stats_args} \
      -O "~{base_file_name}.mutect2.stats" \
      --verbosity WARNING
  >>>

  output {
    File vcf = "~{base_file_name}.mutect2.vcf.gz"
    File vcf_index = "~{base_file_name}.mutect2.vcf.gz.tbi"
    File stats_file = "~{base_file_name}.mutect2.stats"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
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
