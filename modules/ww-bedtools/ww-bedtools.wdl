## WILDS WDL for working with genomic intervals using BEDTools.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task coverage {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
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
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
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
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
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
