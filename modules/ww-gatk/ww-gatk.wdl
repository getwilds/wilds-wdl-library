## WILDS WDL module for GATK variant calling and analysis tasks.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task create_sequence_dictionary {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
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
    max_mnp_distance: "Distance at which to merge MNPs (default: 1)"
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
    Int max_mnp_distance = 1
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
        --max-mnp-distance "~{max_mnp_distance}" \
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

task fastq_to_bam {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Convert paired FASTQ files to unmapped BAM using GATK FastqToSam"
    outputs: {
        unmapped_bam: "Unmapped BAM file containing reads from input FASTQ files"
    }
  }

  parameter_meta {
    r1_fastq: "Array of R1 FASTQ files for the library"
    r2_fastq: "Array of R2 FASTQ files for the library"
    base_file_name: "Base name for output file"
    sample_name: "Sample name to insert into the read group header"
    library_name: "Library name to place into the LB attribute in the read group header (defaults to sample_name if not provided)"
    platform: "Sequencing platform (default: illumina)"
    sequencing_center: "Location where the sample was sequenced (defaults to '.' if not provided)"
    read_group_name: "Read group name (if not provided, defaults to sample_name)"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    String base_file_name
    String sample_name
    String? library_name
    String platform = "illumina"
    String? sequencing_center
    String? read_group_name
    Int memory_gb = 8
    Int cpu_cores = 4
  }

  String rg_name = select_first([read_group_name, sample_name])
  String lib_name = select_first([library_name, sample_name])
  String seq_center = select_first([sequencing_center, "."])

  command <<<
    set -eo pipefail

    gatk --java-options "-Dsamjdk.compression_level=5 -Xms~{memory_gb - 2}g" \
      FastqToSam \
      --FASTQ ~{sep=" " r1_fastq} \
      --FASTQ2 ~{sep=" " r2_fastq} \
      --OUTPUT "~{base_file_name}.unmapped.bam" \
      --READ_GROUP_NAME "~{rg_name}" \
      --SAMPLE_NAME "~{sample_name}" \
      --LIBRARY_NAME "~{lib_name}" \
      --PLATFORM "~{platform}" \
      --SEQUENCING_CENTER "~{seq_center}"
  >>>

  output {
    File unmapped_bam = "~{base_file_name}.unmapped.bam"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task validate_sam_file {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Validate BAM/CRAM/SAM files for formatting issues using GATK ValidateSamFile"
    outputs: {
        validation_report: "Text file containing validation statistics and any errors/warnings"
    }
  }

  parameter_meta {
    input_file: "BAM/CRAM/SAM file to validate"
    base_file_name: "Base name for output validation file"
    mode: "Validation mode: VERBOSE (detailed), SUMMARY (summary only)"
    ignore_warnings: "Whether to ignore warnings (default: false)"
    reference_fasta: "Reference genome FASTA (required for CRAM files)"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File input_file
    String base_file_name
    String mode = "SUMMARY"
    Boolean ignore_warnings = false
    File? reference_fasta
    Int memory_gb = 4
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail

    gatk --java-options "-Dsamjdk.compression_level=5 -Xms~{memory_gb - 2}g" \
      ValidateSamFile \
      --INPUT "~{input_file}" \
      ~{if defined(reference_fasta) then "--REFERENCE_SEQUENCE " + reference_fasta else ""} \
      --MODE ~{mode} \
      --IGNORE_WARNINGS ~{ignore_warnings} > "~{base_file_name}.validation.txt"
  >>>

  output {
    File validation_report = "~{base_file_name}.validation.txt"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task analyze_saturation_mutagenesis {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Analyze saturation mutagenesis data using GATK AnalyzeSaturationMutagenesis"
    outputs: {
        aa_counts: "Amino acid count table",
        aa_fractions: "Amino acid fraction table",
        codon_counts: "Codon count table",
        codon_fractions: "Codon fraction table",
        cov_length_counts: "Coverage length count table",
        read_counts: "Read count table",
        ref_coverage: "Reference coverage table",
        variant_counts: "Variant count table"
    }
  }

  parameter_meta {
    bam: "Input aligned BAM file (will be automatically sorted by queryname)"
    bam_index: "Index file for the input BAM"
    reference_fasta: "Reference genome FASTA file (must contain only A, C, G, T bases - no N's or other IUPAC codes)"
    reference_fasta_index: "Index file for the reference FASTA"
    reference_dict: "Reference genome sequence dictionary"
    orf_range: "Open reading frame range to analyze (e.g., '1-100')"
    base_file_name: "Base name for output files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File bam
    File bam_index
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    String orf_range
    String base_file_name
    Int memory_gb = 16
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail

    # Add local symbolic link for reference fasta and dict
    # If soft links aren't allowed on your HPC system, copy them locally instead
    ln -s "~{reference_fasta}" "~{basename(reference_fasta)}"
    ln -s "~{reference_fasta_index}" "~{basename(reference_fasta_index)}"
    ln -s "~{reference_dict}" "~{basename(reference_dict)}"

    # Sort BAM by queryname so mates are adjacent (required by AnalyzeSaturationMutagenesis)
    samtools sort -n -@ ~{cpu_cores} -o "~{base_file_name}.sorted.bam" "~{bam}"

    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      AnalyzeSaturationMutagenesis \
      -R "~{basename(reference_fasta)}" \
      -I "~{base_file_name}.sorted.bam" \
      --orf ~{orf_range} \
      -O "~{base_file_name}" \
      --verbosity WARNING
  >>>

  output {
    File aa_counts = "~{base_file_name}.aaCounts"
    File aa_fractions = "~{base_file_name}.aaFractions"
    File codon_counts = "~{base_file_name}.codonCounts"
    File codon_fractions = "~{base_file_name}.codonFractions"
    File cov_length_counts = "~{base_file_name}.coverageLengthCounts"
    File read_counts = "~{base_file_name}.readCounts"
    File ref_coverage = "~{base_file_name}.refCoverage"
    File variant_counts = "~{base_file_name}.variantCounts"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task create_somatic_pon {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Create a somatic panel of normals (PON) from VCF files"
    outputs: {
        pon_vcf: "Gzipped VCF file containing the panel of normals",
        pon_vcf_index: "Index file for the panel of normals VCF"
    }
  }

  parameter_meta {
    normal_vcfs: "Array of Mutect2 gzipped VCF files generated from normal samples"
    base_file_name: "Base name for output files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    Array[File] normal_vcfs
    String base_file_name
    Int memory_gb = 4
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail

    # Run CreateSomaticPanelOfNormals directly with VCFs
    gatk --java-options "-Xms~{memory_gb - 4}g -Xmx~{memory_gb - 2}g" \
      CreateSomaticPanelOfNormals \
      -V ~{sep=" -V " normal_vcfs} \
      -O "~{base_file_name}.pon.vcf.gz" \
      --create-output-variant-index true \
      --QUIET true \
      --verbosity WARNING
  >>>

  output {
    File pon_vcf = "~{base_file_name}.pon.vcf.gz"
    File pon_vcf_index = "~{base_file_name}.pon.vcf.gz.tbi"
  }

  runtime {
    docker: "getwilds/gatk:4.6.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}
