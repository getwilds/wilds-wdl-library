## WILDS WDL for working with genomics files using BEDTools.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

struct SampleInfo {
    String name
    File bam
    File bam_index
}

#### WORKFLOW DEFINITION

workflow bedtools_example {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "WDL workflow for manipulating genomic files via BEDTools."
    url: "https://github.com/getwilds/wilds-wdl-library/modules/ww-bedtools"
    version: "1.0"
    outputs: {
        intersect_results: "BEDTools intersect output files for each sample",
        coverage_results: "Read coverage results across BED intervals",
        window_counts: "Tarballs of per-chromosome BED files of read counts"
    }
  }

  parameter_meta {
    bed_file: "BED file containing genomic intervals of interest"
    samples: "Array of sample information containing name, BAM file, and BAM index"
    reference_fasta: "Reference genome FASTA file used for analysis"
    reference_index: "Index file for the reference genome"
    intersect_flags: "Flags for BEDTools intersect command"
    chromosomes: "List of chromosomes to analyze for window-based counting"
    tmp_dir : "Path to a temporary directory"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    memory_gb: "Memory allocated for each task in the workflow in GB"
  }

  input {
    File bed_file
    Array[SampleInfo] samples
    File reference_fasta
    File reference_index
    String intersect_flags = "-header -wo"
    Array[String] chromosomes = [
      "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
      "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
      "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    String tmp_dir = "/tmp"
    Int cpus = 8
    Int memory_gb = 32
  }

  scatter (sample in samples) {
    call coverage { 
      input:
        bed_file = bed_file,
        sample_data = sample,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }

    call intersect { 
      input:
        bed_file = bed_file,
        sample_data = sample,
        flags = intersect_flags,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }

    call makewindows { 
      input:
        bed_file = bed_file,
        sample_data = sample,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        list_chr = chromosomes,
        tmp_dir = tmp_dir,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  call validate_outputs { 
    input:
      intersect_files = intersect.intersect_output,
      coverage_files = coverage.mean_coverage,
      window_count_files = makewindows.counts_bed,
      sample_names = coverage.name
  }

  output {
    Array[File] intersect_results = intersect.intersect_output
    Array[File] coverage_results = coverage.mean_coverage
    Array[File] window_count_results = makewindows.counts_bed
    File bedtools_validation_report = validate_outputs.report
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
    sample_data: "Sample object containing sample name, BAM, and BAM index (.bai) files"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bed_file
    SampleInfo sample_data
    Int cpu_cores = 8
    Int memory_gb = 16
  }

  command <<<
    set -eo pipefail && \
    # Find the current sort order of this bam file
    samtools view -H "~{sample_data.bam}" | grep @SQ | sed 's/@SQ\tSN:\|LN://g' > "~{sample_data.name}.sortOrder.txt" && \
    # Sort the BED file according to BAM sort order
    bedtools sort -g "~{sample_data.name}.sortOrder.txt" -i "~{bed_file}" > correctly.sorted.bed && \
    # Calculate mean coverage
    bedtools coverage -mean -sorted -g "~{sample_data.name}.sortOrder.txt" \
        -a correctly.sorted.bed -b "~{sample_data.bam}" > "~{sample_data.name}.bedtools_mean_covg.txt"
  >>>

  output {
    String name = sample_data.name
    File mean_coverage = "~{sample_data.name}.bedtools_mean_covg.txt"
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
    sample_data: "Sample object containing sample name, BAM, and BAM index (.bai) files"
    flags: "BEDTools intersect command flags"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bed_file
    SampleInfo sample_data
    String flags = "-header -wo"
    Int cpu_cores = 6
    Int memory_gb = 32
  }

  command <<<
    set -eo pipefail && \
    bedtools intersect ~{flags} -a "~{bed_file}" -b "~{sample_data.bam}" > "~{sample_data.name}.intersect.txt"
  >>>

  output {
    String name = sample_data.name
    File intersect_output = "~{sample_data.name}.intersect.txt"
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
    sample_data: "Sample object containing sample name, BAM, and BAM index (.bai) files"
    reference_fasta: "Reference genome FASTA file"
    reference_index: "Reference genome index file"
    list_chr: "Array of chromosome names to process"
    tmp_dir: "Path to a temporary directory"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bed_file
    SampleInfo sample_data
    File reference_fasta
    File reference_index
    Array[String] list_chr
    String tmp_dir
    Int cpu_cores = 10
    Int memory_gb = 24
  }
  
 command <<<
    set -eo pipefail
    mkdir -p ~{sample_data.name}
    
    for Chrom in ~{sep=' ' list_chr}
    do
      # Create windows for this chromosome
      bedtools makewindows -b ~{bed_file} -w 500000 | \
        awk -v OFS="\t" -v C="${Chrom}" '$1==C && NF==3' > ~{tmp_dir}/${Chrom}.windows.bed
      
      # Count reads in windows for this chromosome (run in background for parallelization)
      samtools view -@ 5 -b -f 0x2 -F 0x400 -q 20 -T ~{reference_fasta} ~{sample_data.bam} ${Chrom} | \
        bedtools intersect -sorted -c -a ~{tmp_dir}/${Chrom}.windows.bed -b stdin > "~{sample_data.name}/${Chrom}.counts.bed" &
    done

    wait

    tar -czf "~{sample_data.name}.tar.gz" "~{sample_data.name}"
  >>>

  runtime {
    docker: "getwilds/bedtools:2.31.1"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }

  output {
    String name = sample_data.name
    File counts_bed = "${sample_data.name}.tar.gz"
  }
}

task validate_outputs {
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
