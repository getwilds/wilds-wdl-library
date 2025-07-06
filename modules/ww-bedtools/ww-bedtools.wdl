## WILDS WDL for working with genomic intervals using BEDTools.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bwa/ww-bwa.wdl" as ww_bwa
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-sra/ww-sra.wdl" as ww_sra
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct BedtoolsSample {
    String name
    File bam
    File bam_index
}

#### WORKFLOW DEFINITION

workflow bedtools_example {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "WDL workflow for working with genomic intervals using BEDTools."
    url: "https://github.com/getwilds/wilds-wdl-library/modules/ww-bedtools"
    version: "1.0"
    outputs: {
        intersect_results: "BEDTools intersect output files for each sample",
        coverage_results: "Read coverage results across BED intervals",
        window_counts: "Tarballs of per-chromosome BED files of read counts",
        validation_report: "Validation report summarizing file check results"
    }
  }

  parameter_meta {
    bed_file: "BED file containing genomic intervals of interest"
    samples: "Array of sample information containing name, BAM file, and BAM index"
    reference_fasta: "Reference genome FASTA file used for analysis"
    reference_index: "Index file for the reference genome"
    intersect_flags: "Flags for BEDTools intersect command"
    demo_sra_id: "SRA accession ID to use for demonstration when no samples are provided"
    chromosomes: "List of chromosomes to analyze for window-based counting"
    tmp_dir : "Path to a temporary directory"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    memory_gb: "Memory allocated for each task in the workflow in GB"
  }

  input {
    File? bed_file
    Array[BedtoolsSample]? samples
    File? reference_fasta
    File? reference_index
    String demo_sra_id = "ERR1258306"
    Array[String] chromosomes = [
      "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
      "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
      "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
      ]
    String intersect_flags = "-header -wo"
    String tmp_dir = "/tmp"
    Int cpus = 2
    Int memory_gb = 16
  }

  # If no reference genome provided, download test data
  if (!defined(reference_fasta) || !defined(reference_index) || !defined(bed_file)) {
    call ww_testdata.download_ref_data { }
  }

  # Determine which genome files to use
  File genome_fasta = select_first([reference_fasta, download_ref_data.fasta])
  File genome_fasta_index = select_first([reference_index, download_ref_data.fasta_index])
  File bed_file_final = select_first([bed_file, download_ref_data.bed])

  # If no samples provided, download demonstration data from SRA and align with BWA
  if (!defined(samples)) {
    call ww_sra.fastqdump { input:
        sra_id = demo_sra_id,
        ncpu = cpus
    }

    # Build BWA index for alignment
    call ww_bwa.bwa_index { input:
        reference_fasta = genome_fasta,
        cpu_cores = cpus,
        memory_gb = memory_gb * 2
    }

    # Align the SRA sample using BWA
    call ww_bwa.bwa_mem { input:
        bwa_genome_tar = bwa_index.bwa_index_tar,
        reference_fasta = genome_fasta,
        sample_data = {
          "name": demo_sra_id,
          "r1": fastqdump.r1_end,
          "r2": fastqdump.r2_end
        },
        cpu_cores = cpus,
        memory_gb = memory_gb * 2
    }
  }

  # Create samples array - either from input or from BWA alignment
  Array[BedtoolsSample] final_samples = if defined(samples) then select_first([samples]) else [
    {
      "name": demo_sra_id,
      "bam": select_first([bwa_mem.sorted_bam]),
      "bai": select_first([bwa_mem.sorted_bai])
    }
  ]

  scatter (sample in final_samples) {
    call coverage { input:
        bed_file = bed_file_final,
        aligned_bam = sample.bam,
        sample_name = sample.name,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }

    call intersect { input:
        bed_file = bed_file_final,
        aligned_bam = sample.bam,
        sample_name = sample.name,
        flags = intersect_flags,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }

    call makewindows { input:
        bed_file = bed_file_final,
        aligned_bam = sample.bam,
        bam_index = sample.bam_index,
        sample_name = sample.name,
        reference_fasta = genome_fasta,
        reference_index = genome_fasta_index,
        list_chr = chromosomes,
        tmp_dir = tmp_dir,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  call validate_outputs { input:
      intersect_files = intersect.intersect_output,
      coverage_files = coverage.mean_coverage,
      window_count_files = makewindows.counts_bed,
      sample_names = coverage.name
  }

  output {
    Array[File] intersect_results = intersect.intersect_output
    Array[File] coverage_results = coverage.mean_coverage
    Array[File] window_counts = makewindows.counts_bed
    File validation_report = validate_outputs.report
  }
}

#### TASK DEFINITIONS

task coverage {
  meta {
    description: "Task that gets BAM coverage across BED intervals using BEDTools"
    outputs: {
        name: "Sample name that was processed",
        mean_coverage: "File containing mean read coverage across BED intervals"
    }
  }

  parameter_meta {
    bed_file:  "BED file containing genomic intervals"
    aligned_bam: "Input aligned and indexed BAM file"
    sample_name: "Name of the sample provided for output files"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bed_file
    File aligned_bam
    String sample_name
    Int cpu_cores = 2
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail && \
    # Find the current sort order of this bam file
    samtools view -H "~{aligned_bam}" | grep @SQ | sed 's/@SQ\tSN:\|LN://g' > "~{sample_name}.sortOrder.txt" && \
    # Sort the BED file according to BAM sort order
    bedtools sort -g "~{sample_name}.sortOrder.txt" -i "~{bed_file}" > correctly.sorted.bed && \
    # Calculate mean coverage
    bedtools coverage -mean -sorted -g "~{sample_name}.sortOrder.txt" \
        -a correctly.sorted.bed -b "~{aligned_bam}" > "~{sample_name}.bedtools_mean_covg.txt"
  >>>

  output {
    String name = sample_name
    File mean_coverage = "~{sample_name}.bedtools_mean_covg.txt"
  }

  runtime {
    docker: "getwilds/bedtools:2.31.1"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task intersect {
  meta {
    description: "Task that performs BEDTools intersect between two genomic files"
    outputs: {
        name: "Sample name that was processed",
        intersect_output: "BEDTools intersect results file"
    }
  }

  parameter_meta {
    bed_file: "BED file containing genomic intervals"
    aligned_bam: "Input aligned and indexed BAM file"
    sample_name: "Name of the sample provided for output files"
    flags: "BEDTools intersect command flags"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bed_file
    File aligned_bam
    String sample_name
    String flags = "-header -wo"
    Int cpu_cores = 2
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail && \
    bedtools intersect ~{flags} -a "~{bed_file}" -b "~{aligned_bam}" > "~{sample_name}.intersect.txt"
  >>>

  output {
    String name = sample_name
    File intersect_output = "~{sample_name}.intersect.txt"
  }

  runtime {
    docker: "getwilds/bedtools:2.31.1"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task makewindows {
  meta {
    description: "Task that creates genomic windows and counts reads per window across chromosomes"
    outputs: {
        name: "Sample name that was processed",
        counts_bed: "Tarball of per-chromosome BED files of read counts"
    }
  }

  parameter_meta {
    bed_file: "BED file containing genomic intervals"
    reference_fasta: "Reference genome FASTA file"
    reference_index: "Reference genome index file"
    aligned_bam: "Input aligned BAM file"
    bam_index: "Index of aligned BAM file"
    list_chr: "Array of chromosome names to process"
    sample_name: "Name of the sample provided for output files"
    tmp_dir: "Path to a temporary directory"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bed_file
    File reference_fasta
    File reference_index
    File aligned_bam
    File bam_index
    Array[String] list_chr
    String sample_name
    String tmp_dir
    Int cpu_cores = 10
    Int memory_gb = 24
  }

 command <<<
    set -eo pipefail
    mkdir -p "~{sample_name}"

    for Chrom in ~{sep=" " list_chr}
    do
      # Create windows for this chromosome
      bedtools makewindows -b "~{bed_file}" -w 500000 | \
        awk -v OFS="\t" -v C="${Chrom}" '$1==C && NF==3' > "~{tmp_dir}"/"${Chrom}".windows.bed

      # Count reads in windows for this chromosome (run in background for parallelization)
      samtools view -@ 5 -b -f 0x2 -F 0x400 -q 20 -T "~{reference_fasta}" "~{aligned_bam}" "${Chrom}" | \
        bedtools intersect -sorted -c -a "~{tmp_dir}"/"${Chrom}".windows.bed -b stdin > "~{sample_name}/${Chrom}.counts.bed" &
    done

    wait

    tar -czf "~{sample_name}.tar.gz" "~{sample_name}"
  >>>

  output {
    String name = sample_name
    File counts_bed = "${sample_name}.tar.gz"
  }

  runtime {
    docker: "getwilds/bedtools:2.31.1"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task validate_outputs {
  # TODO: Do a basic check of the file contents too
  meta {
    description: "Validate that BEDTools output files exist and are non-empty"
    outputs: {
        report: "Validation report summarizing file check results"
    }
  }

  parameter_meta {
    intersect_files: "BEDTools intersect output files for each sample"
    coverage_files: "Coverage analysis results across BED intervals"
    window_count_files: "Tarballs of per-chromosome BED files of read counts"
    sample_names: "Array of sample names that were processed"
  }

  input {
    Array[File] intersect_files
    Array[File] coverage_files
    Array[File] window_count_files
    Array[String] sample_names
  }

  command <<<
    set -eo pipefail

    echo "=== BEDTools Workflow Validation Report ===" > bedtools_validation_report.txt
    echo "" >> bedtools_validation_report.txt

    intersect_files=("~{sep=" " intersect_files}")
    coverage_files=("~{sep=" " coverage_files}")
    window_count_files=("~{sep=" " window_count_files}")
    sample_names=("~{sep=" " sample_names}")

    validation_passed=true

    for i in "${!sample_names[@]}"; do
      sample="${sample_names[$i]}"
      intersect="${intersect_files[$i]}"
      coverage="${coverage_files[$i]}"
      tarball="${window_count_files[$i]}"

      echo "--- Sample: $sample ---" >> bedtools_validation_report.txt

      # Intersect file
      if [[ -f "$intersect" && -s "$intersect" ]]; then
        size=$(stat -c%s "$intersect")
        lines=$(wc -l < "$intersect")
        echo "  Intersect: PASS (${size} bytes, ${lines} lines)" >> bedtools_validation_report.txt
      else
        echo "  Intersect: FAIL - MISSING OR EMPTY" >> bedtools_validation_report.txt
        validation_passed=false
      fi

      # Coverage file
      if [[ -f "$coverage" && -s "$coverage" ]]; then
        size=$(stat -c%s "$coverage")
        lines=$(wc -l < "$coverage")
        echo "  Coverage: PASS (${size} bytes, ${lines} lines)" >> bedtools_validation_report.txt
      else
        echo "  Coverage: FAIL - MISSING OR EMPTY" >> bedtools_validation_report.txt
        validation_passed=false
      fi

      # Window count tarball
      if [[ -f "$tarball" && -s "$tarball" ]]; then
        size=$(stat -c%s "$tarball")
        echo "  Windows BED file Tarball: PASS (${size} bytes)" >> bedtools_validation_report.txt
      else
        echo "  Window BED file Tarball: FAIL - MISSING OR EMPTY" >> bedtools_validation_report.txt
        validation_passed=false
      fi

      echo "" >> bedtools_validation_report.txt
    done

    echo "=== Validation Summary ===" >> bedtools_validation_report.txt
    echo "Total samples processed: ${#sample_names[@]}" >> bedtools_validation_report.txt
    if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> bedtools_validation_report.txt
    else
      echo "Overall Status: FAILED" >> bedtools_validation_report.txt
      exit 1
    fi

    # Also output to stdout for immediate feedback
    cat bedtools_validation_report.txt
  >>>

  output {
    File report = "bedtools_validation_report.txt"
  }

  runtime {
    docker: "getwilds/bedtools:2.31.1"
    cpu: 1
    memory: "2 GB"
  }
}
