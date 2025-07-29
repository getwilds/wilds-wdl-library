## WILDS WDL for processing genomic files with Samtools.
## Intended for use alone or as a modular component in the WILDS ecosystem.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct SampleInfo {
  String name
  Array[File] cram_files
}

workflow samtools_example {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "WDL workflow for processing genomic files with Samtools"
    url: "https://github.com/getwilds/wilds-wdl-library/modules/ww-samtools"
    outputs: {
        r1_fastqs: "Array of R1 FASTQ files generated from CRAM/BAM/SAM files",
        r2_fastqs: "Array of R2 FASTQ files generated from CRAM/BAM/SAM files",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    samples: "Array of sample objects, each containing name and paths to CRAM/BAM/SAM files"
    reference_fasta: "Reference genome FASTA file"
    interval_bed: "Optional BED file defining genomic intervals for splitting BAM files"
    cpus: "Number of CPU cores allocated for each non-validation task"
    memory_gb: "Memory in GB allocated for each non-validation task"
  }

  input {
    Array[SampleInfo]? samples
    File? reference_fasta
    File? interval_bed
    Int cpus = 2
    Int memory_gb = 8
  }

  # If no reference genome provided, download test data
  if (!defined(reference_fasta) || !defined(interval_bed)) {
    call ww_testdata.download_ref_data { }
  }

  # Determine which genome and interval files to use
  File genome_fasta = select_first([reference_fasta, download_ref_data.fasta])
  File intervals = select_first([interval_bed, download_ref_data.bed])

  # If no samples provided, download demonstration data
  if (!defined(samples)) {
    call ww_testdata.download_cram_data { input:
        ref_fasta = genome_fasta
    }
  }

  # Create test sample array when no samples provided
  Array[SampleInfo] test_samples = if defined(download_cram_data.cram) then [
    object {
      name: "test_sample",
      cram_files: [
        select_first([download_cram_data.cram]),
        select_first([download_cram_data.cram]),
        select_first([download_cram_data.cram])
      ]
    }
  ] else []

  # Create the samples array - either from input or from test data
  Array[SampleInfo] final_samples = select_first([samples, test_samples])

  scatter (sample in final_samples) {
    call crams_to_fastq { input:
        cram_files = sample.cram_files,
        ref = genome_fasta,
        name = sample.name,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  # Testing splitting BAM files by intervals
  call ww_testdata.download_bam_data { }
  call split_bam_by_intervals { input:
      input_bam = download_bam_data.bam,
      input_bam_index = download_bam_data.bai,
      bed_files = [intervals],
      output_basename = "split_intervals_test",
      memory_gb = memory_gb,
      cpu_cores = cpus
  }

  call validate_outputs { input:
      r1_fastqs = crams_to_fastq.r1_fastq,
      r2_fastqs = crams_to_fastq.r2_fastq,
      sample_names = crams_to_fastq.sample_name
  }

  output {
    Array[File] r1_fastqs = crams_to_fastq.r1_fastq
    Array[File] r2_fastqs = crams_to_fastq.r2_fastq
    Array[File] interval_bams = split_bam_by_intervals.interval_bams
    Array[File] interval_bam_indices = split_bam_by_intervals.interval_bam_indices
    File validation_report = validate_outputs.report
  }
}

task crams_to_fastq {
  meta {
    description: "Merge CRAM/BAM/SAM files and convert to FASTQ's using samtools."
    outputs: {
        r1_fastq: "R1 FASTQ file generated from merged CRAM/BAM/SAM file",
        r2_fastq: "R2 FASTQ file generated from merged CRAM/BAM/SAM file",
        sample_name: "Sample name that was processed"
    }
  }

  parameter_meta {
    cram_files: "List of CRAM/BAM/SAM files for a sample"
    ref: "Reference genome FASTA file"
    name: "Name of the sample (used for file naming)"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    Array[File] cram_files
    File ref
    String name
    Int cpu_cores = 2
    Int memory_gb = 16
  }

  command <<<
    # Merge CRAM/BAM/SAM files if more than one, then convert to FASTQ
    samtools merge -@ ~{cpu_cores} --reference "~{ref}" -f "~{name}.merged.cram" ~{sep=" " cram_files} && \
    samtools collate -u -O "~{name}.merged.cram" | \
    samtools fastq --reference "~{ref}" -1 "~{name}_R1.fastq.gz" -2 "~{name}_R2.fastq.gz" -0 /dev/null -s /dev/null -

    # Cleaning up merged CRAM file to save space
    rm -f "~{name}.merged.cram"
  >>>

  output {
    File r1_fastq = "~{name}_R1.fastq.gz"
    File r2_fastq = "~{name}_R2.fastq.gz"
    String sample_name = "~{name}"
  }

  runtime {
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
    docker: "getwilds/samtools:1.11"
  }
}

task split_bam_by_intervals {
  meta {
    description: "Split BAM file by genomic intervals using samtools view for scatter-gather parallelization"
    outputs: {
        interval_bams: "Array of BAM files split by intervals",
        interval_bam_indices: "Array of BAM index files for the split BAMs"
    }
  }

  parameter_meta {
    input_bam: "Input BAM file to split"
    input_bam_index: "Index file for the input BAM"
    bed_files: "Array of BED files defining intervals to split by"
    output_basename: "Base name for output BAM files"
    memory_gb: "Memory allocation in GB"
    cpu_cores: "Number of CPU cores to use"
  }

  input {
    File input_bam
    File input_bam_index
    Array[File] bed_files
    String output_basename
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail
    
    # Create array to track output files
    declare -a bam_files
    declare -a bai_files
    
    # Process each BED file
    for bed_file in ~{sep=" " bed_files}; do
      # Extract interval name for output filename
      interval_name=$(basename "$bed_file" .bed)
      output_bam="~{output_basename}.${interval_name}.bam"
      
      # Use samtools view with BED file directly (-L flag)
      samtools view -b -h -L "$bed_file" "~{input_bam}" > "$output_bam"
      
      # Index the resulting BAM
      samtools index "$output_bam"
      
      # Add to arrays for output
      bam_files+=("$output_bam")
      bai_files+=("${output_bam}.bai")
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
    docker: "getwilds/samtools:1.11"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task validate_outputs {
  meta {
    description: "Validates that FASTQ files exist and are non-empty after CRAM-to-FASTQ conversion."
    outputs: {
        report: "Validation report with pass/fail summary"
    }
  }

  parameter_meta {
    r1_fastqs: "Array of R1 FASTQ files to check"
    r2_fastqs: "Array of R2 FASTQ files to check"
    sample_names: "Array of sample names corresponding to FASTQ files"
  }

  input {
    Array[File] r1_fastqs
    Array[File] r2_fastqs
    Array[String] sample_names
  }

  command <<<
    set -eo pipefail

    echo "=== Samtools FASTQ Validation Report ===" > samtools_validation_report.txt
    echo "" >> samtools_validation_report.txt

    r1_fastqs=(~{sep=" " r1_fastqs})
    r2_fastqs=(~{sep=" " r2_fastqs})
    sample_names=(~{sep=" " sample_names})

    validation_passed=true

    for i in "${!sample_names[@]}"; do
      sample="${sample_names[$i]}"
      r1="${r1_fastqs[$i]}"
      r2="${r2_fastqs[$i]}"

      echo "--- Sample: $sample ---" >> samtools_validation_report.txt

      if [[ -f "$r1" && -s "$r1" ]]; then
        size=$(stat -c%s "$r1")
        lines=$(zcat "$r1" | wc -l)
        echo "  R1 FASTQ: PASS (${size} bytes, ${lines} lines)" >> samtools_validation_report.txt
      else
        echo "  R1 FASTQ: FAIL - MISSING OR EMPTY" >> samtools_validation_report.txt
        validation_passed=false
      fi

      if [[ -f "$r2" && -s "$r2" ]]; then
        size=$(stat -c%s "$r2")
        lines=$(zcat "$r2" | wc -l)
        echo "  R2 FASTQ: PASS (${size} bytes, ${lines} lines)" >> samtools_validation_report.txt
      else
        echo "  R2 FASTQ: FAIL - MISSING OR EMPTY" >> samtools_validation_report.txt
        validation_passed=false
      fi

      echo "" >> samtools_validation_report.txt
    done

    echo "=== Summary ===" >> samtools_validation_report.txt
    echo "Samples processed: ${#sample_names[@]}" >> samtools_validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
      echo "Status: PASSED" >> samtools_validation_report.txt
    else
      echo "Status: FAILED" >> samtools_validation_report.txt
      exit 1
    fi

    cat samtools_validation_report.txt
  >>>

  output {
    File report = "samtools_validation_report.txt"
  }

  runtime {
    docker: "getwilds/samtools:1.11"
    cpu: 1
    memory: "2 GB"
  }
}
