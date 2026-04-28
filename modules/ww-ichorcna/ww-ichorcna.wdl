## WILDS WDL module for tumor fraction estimation in cfDNA using ichorCNA.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task readcounter_wig {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Generate tumor WIG file from aligned bam files using HMMcopy's readCounter"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-ichorcna/ww-ichorcna.wdl"
    outputs: {
        wig_file: "WIG file created from binned read count data within input BED files"
    }
  }

  parameter_meta {
    bam_file: "Aligned bam file containing reads to be analyzed"
    bam_index: "Index for the bam file"
    sample_name: "Name of the sample being analyzed"
    window_size: "Window size in base pairs for WIG format"
    chromosomes: "Chromosomes to include in WIG file"
    memory_gb: "Memory allocated for the task in GB"
    cpus: "Number of CPU cores allocated for the task"
  }

  input {
    File bam_file
    File bam_index
    String sample_name
    Array[String] chromosomes
    Int window_size = 500000
    Int memory_gb = 8
    Int cpus = 2
  }

  command <<<
    set -eo pipefail
    
    # Use readCounter from HMMcopy to generate WIG directly from BAM
    readCounter \
      --window ~{window_size} \
      --quality 20 \
      --chromosome ~{sep="," chromosomes} \
      "~{bam_file}" > "~{sample_name}.wig"
  >>>

  output {
    File wig_file = "${sample_name}.wig"
  }

  runtime {
    docker: "getwilds/hmmcopy:1.0.0"
    cpu: cpus
    memory: "~{memory_gb} GB"
  }
}

task ichorcna_call {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "Estimate cfDNA tumor fraction using ichorCNA"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-ichorcna/ww-ichorcna.wdl"
    outputs: {
        params: "Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions",
        seg: "Segments called by the Viterbi algorithm, including subclonal status of segments (0=clonal, 1=subclonal), and filtered fo exclude Y chromosome segments if not male",
        genomewide_pdf: "Genome wide plot of data annotated for estimated copy number, tumor fraction, and ploidy for the optimal solution",
        allgenomewide_pdf: "Combined PDF of all solutions",
        correct_pdf:  "Genome wide correction comparisons",
        rdata: "Saved R image after ichorCNA has finished. Results for all solutions will be included"
    }
  }

  parameter_meta {
    wig_tumor: "Tumor WIG file being analyzed"
    wig_gc: "GC-content WIG file"
    wig_map: "Mappability score WIG file"
    panel_of_norm_rds: "RDS file of median corrected depth from panel of normals"
    centromeres: "Text file containing Centromere locations"
    name: "Sample ID"
    sex: "User-specified: male or female"
    genome: "Genome build (e.g. hg19)"
    genome_style: "Chromosome naming convention (use UCSC if desired output is to have 'chr' string): NCBI or UCSC"
    chrs: "Chromosomes to analyze (default: chr 1-22, X, and Y)"
    memory_gb: "Memory allocated for each task in the workflow in GB"
    cpus: "Number of CPU cores allocated for each task in the workflow"
  }

  input {
    File wig_tumor
    File wig_gc
    File wig_map
    File panel_of_norm_rds
    File centromeres
    String name
    String sex
    String genome
    String genome_style
    String chrs = "c(1:22, 'X', 'Y')"
    Int memory_gb = 16
    Int cpus = 6
  }

  command <<<
    set -eo pipefail

    Rscript /usr/local/bin/ichorCNA/scripts/runIchorCNA.R \
    --id "~{name}" \
    --WIG "~{wig_tumor}" \
    --ploidy "c(2)" \
    --normal "c(0.1,0.5,.85)" \
    --maxCN 3 \
    --gcWig "~{wig_gc}" \
    --mapWig "~{wig_map}" \
    --centromere "~{centromeres}" \
    --normalPanel "~{panel_of_norm_rds}" \
    --genomeBuild "~{genome}" \
    --sex "~{sex}" \
    --chrs "~{chrs}" \
    --fracReadsInChrYForMale 0.0005 \
    --txnE 0.999999 \
    --txnStrength 1000000 \
    --genomeStyle "~{genome_style}" \
    --libdir /usr/local/bin/ichorCNA/ && \
    # Keep only biologically relevant segments by keeping all Y segments if
    # male and only non-Y segments if female
    awk -v G="~{sex}" '$2!~/Y/ || G=="male"' "~{name}.seg.txt" \
    > "~{name}.ichor.segs.txt" && \
    mv "~{name}"/*.pdf .
  >>>

  output {
    File params = "${name}.params.txt"
    File seg = "${name}.ichor.segs.txt"
    File genomewide_pdf = "${name}_genomeWide.pdf"
    File allgenomewide_pdf = "${name}_genomeWide_all_sols.pdf"
    File correct_pdf = "${name}_genomeWideCorrection.pdf"
    File rdata = "${name}.RData"
  }

  runtime {
    docker: "getwilds/ichorcna:0.2.0"
    cpu: cpus
    memory: "~{memory_gb} GB"
  }
}
