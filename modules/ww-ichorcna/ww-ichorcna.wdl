## WILDS WDL module for tumor fraction estimation in cfDNA using ichorCNA.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

struct IchorSample {
    String name
    File bam
    File bam_index
}

workflow ichorcna_example {
  meta {
    author: "Emma Bishop"
    email: "ebishop@fredhutch.org"
    description: "WDL workflow for cfDNA tumor fraction estimation via ichorCNA"
    url: "https://github.com/getwilds/modules/ww-ichorcna"
    outputs: {
        params: "Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions",
        seg: "Segments called by the Viterbi algorithm, including subclonal status of segments (0=clonal, 1=subclonal), and filtered fo exclude Y chromosome segments if not male",
        genomewide_pdf: "Genome wide plot of data annotated for estimated copy number, tumor fraction, and ploidy for the optimal solution",
        allgenomewide_pdf: "Combined PDF of all solutions",
        correct_pdf:  "Genome wide correction comparisons",
        rdata: "Saved R image after ichorCNA has finished. Results for all solutions will be included",
        wig: "WIG file created from binned read count data within input BED files",
        validation_report: "Validation report summarizing file check results"
    }
  }

  parameter_meta {
    samples: "Array of sample information containing name and tarball of per-chromosome BED files of read counts"
    wig_gc: "GC-content WIG file"
    wig_map: "Mappability score WIG file"
    panel_of_norm_rds: "RDS file of median corrected depth from panel of normals"
    centromeres: "Text file containing Centromere locations"
    sex: "User-specified: male or female"
    genome: "Genome build (e.g. hg38)"
    genome_style: "Chromosome naming convention (use UCSC if desired output is to have 'chr' string): NCBI or UCSC"
    memory_gb: "Memory allocated for each task in the workflow in GB"
    cpus: "Number of CPU cores allocated for each task in the workflow"
    chrs_list: "Chromosomes to analyze as an array of strings (default: chr 1-22, X, and Y)"
    chrs_vec: "Chromosomes to analyze as an R vector (default: chr 1-22, X, and Y)"
  }

  input {
    Array[IchorSample]? samples
    File? wig_gc
    File? wig_map
    File? panel_of_norm_rds
    File? centromeres
    String sex = "male"
    String genome = "hg38"
    String genome_style = "UCSC"
    Int memory_gb = 8
    Int cpus = 2
    Array[String] chrs_list = ["chr1"] # Limiting to chr1 for test workflow
    String chrs_vec = "c(1)" # Limiting to chr1 for test workflow
  }

  # Determine which ichorCNA data files to use
  if (!defined(wig_gc) || !defined(wig_map) || !defined(panel_of_norm_rds) || !defined(centromeres)) {
    call ww_testdata.download_ichor_data { }
  }
  File final_wig_gc = select_first([wig_gc, download_ichor_data.wig_gc])
  File final_wig_map = select_first([wig_map, download_ichor_data.wig_map])
  File final_panel_of_norm = select_first([panel_of_norm_rds, download_ichor_data.panel_of_norm_rds])
  File final_centromeres = select_first([centromeres, download_ichor_data.centromeres])

  # If no samples provided, download demonstration data
  if (!defined(samples)) {
    call ww_testdata.download_bam_data { }
  }

  # Create samples array - either from input or from test data download
  Array[IchorSample] final_samples = if defined(samples) then select_first([samples]) else [
    {
      "name": "demo_sample",
      "bam": select_first([download_bam_data.bam]),
      "bam_index": select_first([download_bam_data.bai])
    }
  ]

  scatter (sample in final_samples) {
    call readcounter_wig { input:
      bam_file = sample.bam,
      bam_index = sample.bam_index,
      sample_name = sample.name,
      window_size = 500000,
      chromosomes = chrs_list,
      memory_gb = memory_gb,
      cpus = cpus
    }

    call ichorcna_call { input:
      wig_tumor = readcounter_wig.wig_file,
      wig_gc = final_wig_gc,
      wig_map = final_wig_map,
      panel_of_norm_rds = final_panel_of_norm,
      centromeres = final_centromeres,
      name = sample.name,
      sex = sex,
      chrs = chrs_vec,
      genome = genome,
      genome_style = genome_style,
      cpus = cpus,
      memory_gb = memory_gb
    }
  }

  call validate_outputs { input:
      params_files = ichorcna_call.params,
      seg_files = ichorcna_call.seg,
      genome_pdfs = ichorcna_call.genomewide_pdf,
      allgenome_pdfs = ichorcna_call.allgenomewide_pdf,
      correct_pdfs = ichorcna_call.correct_pdf,
      rdata_files = ichorcna_call.rdata,
      wig_files = readcounter_wig.wig_file,
  }

  output {
    Array[File] params = ichorcna_call.params
    Array[File] seg = ichorcna_call.seg
    Array[File] genomewide_pdf = ichorcna_call.genomewide_pdf
    Array[File] allgenomewide_pdf = ichorcna_call.allgenomewide_pdf
    Array[File] correct_pdf = ichorcna_call.correct_pdf
    Array[File] rdata = ichorcna_call.rdata
    Array[File] wig = readcounter_wig.wig_file
    File validation_report = validate_outputs.report
  }
}

task readcounter_wig {
  meta {
    description: "Generate tumor WIG file from aligned bam files using HMMcopy's readCounter"
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
    description: "Estimate cfDNA tumor fraction using ichorCNA"
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
    # TODO: Take a user-generated WIG file as input, instead of making one from BED files
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

task validate_outputs {
  # TODO: Do a basic check of the file contents too
  meta {
    description: "Validate that ichorCNA outputs are non-empty and generate report"
    outputs: {
        report: "Validation report summarizing file check results"
    }
  }

  parameter_meta {
    params_files: "Array of parameter files"
    seg_files: "Array of segment files"
    genome_pdfs: "Array of genome wide plot PDFs"
    allgenome_pdfs: "Array of combined plot PDFs"
    correct_pdfs: "Array of correction comparison PDFs"
    rdata_files: "Array of ichorCNA RData files"
    wig_files: "Array of WIG files generated from BED file data"
  }

  input {
    Array[File] params_files
    Array[File] seg_files
    Array[File] genome_pdfs
    Array[File] allgenome_pdfs
    Array[File] correct_pdfs
    Array[File] rdata_files
    Array[File] wig_files
  }

  command <<<
  set -eo pipefail
  echo "=== IchorCNA Validation Report ===" > ichor_validation_report.txt
  echo "" >> ichor_validation_report.txt

  params_files=("~{sep=" " params_files}")
  seg_files=("~{sep=" " seg_files}")
  genome_pdfs=("~{sep=" " genome_pdfs}")
  allgenome_pdfs=("~{sep=" " allgenome_pdfs}")
  correct_pdfs=("~{sep=" " correct_pdfs}")
  rdata_files=("~{sep=" " rdata_files}")
  wig_files=("~{sep=" " wig_files}")

  validation_passed=true

  for i in "${!wig_files[@]}"; do
      params="${params_files[$i]}"
      seg="${seg_files[$i]}"
      genomewide_pdf="${genome_pdfs[$i]}"
      allgenomewide_pdf="${allgenome_pdfs[$i]}"
      correct_pdf="${correct_pdfs[$i]}"
      rdata="${rdata_files[$i]}"
      wig="${wig_files[$i]}"

      echo "--- Sample: $wig ---" >> ichor_validation_report.txt

      # Parameters file
      if [[ -f "$params" && -s "$params" ]]; then
          size=$(stat -c%s "$params")
          lines=$(wc -l < "$params")
          echo " Parameters: PASS (${size} bytes, ${lines} lines)" \
            >> ichor_validation_report.txt
      else
          echo " Parameters: FAIL - MISSING OR EMPTY" \
            >> ichor_validation_report.txt
          validation_passed=false
      fi

      # Segments file
      if [[ -f "$seg" && -s "$seg" ]]; then
          size=$(stat -c%s "$seg")
          lines=$(wc -l < "$seg")
          echo " Segments: PASS (${size} bytes, ${lines} lines)" \
            >> ichor_validation_report.txt
      else
          echo " Segments: FAIL - MISSING OR EMPTY" >> ichor_validation_report.txt
          validation_passed=false
      fi

      # Genome-wide PDF
      if [[ -f "$genomewide_pdf" && -s "$genomewide_pdf" ]]; then
          size=$(stat -c%s "$genomewide_pdf")
          echo " Genome-wide PDF: PASS (${size} bytes)" \
            >> ichor_validation_report.txt
      else
          echo " Genome-wide PDF: FAIL - MISSING OR EMPTY" \
            >> ichor_validation_report.txt
          validation_passed=false
      fi

      # All genome-wide solutions PDF
      if [[ -f "$allgenomewide_pdf" && -s "$allgenomewide_pdf" ]]; then
          size=$(stat -c%s "$allgenomewide_pdf")
          echo " All Solutions PDF: PASS (${size} bytes)" \
            >> ichor_validation_report.txt
      else
          echo " All Solutions PDF: FAIL - MISSING OR EMPTY" \
            >> ichor_validation_report.txt
          validation_passed=false
      fi

      # Correction PDF
      if [[ -f "$correct_pdf" && -s "$correct_pdf" ]]; then
          size=$(stat -c%s "$correct_pdf")
          echo " Correction PDF: PASS (${size} bytes)" \
            >> ichor_validation_report.txt
      else
          echo " Correction PDF: FAIL - MISSING OR EMPTY" \
            >> ichor_validation_report.txt
          validation_passed=false
      fi

      # RData file
      if [[ -f "$rdata" && -s "$rdata" ]]; then
          size=$(stat -c%s "$rdata")
          echo " RData: PASS (${size} bytes)" >> ichor_validation_report.txt
      else
          echo " RData: FAIL - MISSING OR EMPTY" >> ichor_validation_report.txt
          validation_passed=false
      fi

      # WIG file
      if [[ -f "$wig" && -s "$wig" ]]; then
          size=$(stat -c%s "$wig")
          lines=$(wc -l < "$wig")
          echo " WIG: PASS (${size} bytes, ${lines} lines)" \
            >> ichor_validation_report.txt
      else
          echo " WIG: FAIL - MISSING OR EMPTY" >> ichor_validation_report.txt
          validation_passed=false
      fi

      echo "" >> ichor_validation_report.txt
  done

  echo "=== Validation Summary ===" >> ichor_validation_report.txt
  echo "Total samples processed: ${#wig_files[@]}" \
    >> ichor_validation_report.txt

  if [[ "$validation_passed" == "true" ]]; then
      echo "Overall Status: PASSED" >> ichor_validation_report.txt
  else
      echo "Overall Status: FAILED" >> ichor_validation_report.txt
      exit 1
  fi

  # Also output to stdout for immediate feedback
  cat ichor_validation_report.txt
  >>>

  output {
    File report = "ichor_validation_report.txt"
  }

  runtime {
    docker: "getwilds/ichorcna:0.2.0"
    memory: "4GB"
    cpu: 1
  }
}
