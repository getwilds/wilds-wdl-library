## WILDS WDL module for tumor fraction estimation in cfDNA using ichorCNA.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-bedtools/ww-bedtools.wdl" as ww_bedtools
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
    bed_file: "BED file containing genomic intervals of interest"
    reference_fasta: "Reference genome FASTA file used for analysis"
    reference_index: "Index file for the reference genome"
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
    File? bed_file
    File? reference_fasta
    File? reference_index
    File? wig_gc
    File? wig_map
    File? panel_of_norm_rds
    File? centromeres
    String sex = "male"
    String genome = "hg38"
    String genome_style = "UCSC"
    Int memory_gb = 8
    Int cpus = 2
    Array[String] chrs_list = ["chr1"]
    String chrs_vec = "c(1)"
  }

  String tmp_dir = "/tmp"

  # Determine which genome files to use
  if (!defined(reference_fasta) || !defined(reference_index) || !defined(bed_file)) {
    call ww_testdata.download_ref_data { }
  }
  File genome_fasta = select_first([reference_fasta, download_ref_data.fasta])
  File genome_fasta_index = select_first([reference_index, download_ref_data.fasta_index])
  File final_bed_file = select_first([bed_file, download_ref_data.bed])

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
    call ww_bedtools.makewindows { input:
        bed_file = final_bed_file,
        aligned_bam = sample.bam,
        bam_index = sample.bam_index,
        sample_name = sample.name,
        reference_fasta = genome_fasta,
        reference_index = genome_fasta_index,
        list_chr = chrs_list,
        tmp_dir = tmp_dir,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }

    call generate_tumor_wig { input:
        counts_bed = makewindows.counts_bed,
        sample_name = sample.name,
        chrs = chrs_vec,
        genome_style = genome_style,
        window_size = 500000,
        memory_gb = memory_gb,
        cpus = cpus
    }

    call ichorcna_call { input:
      tumor_wig = generate_tumor_wig.wig,
      wig_gc = final_wig_gc,
      wig_map = final_wig_map,
      panel_of_norm_rds = final_panel_of_norm,
      centromeres = final_centromeres,
      sample_name = sample.name,
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
      wig_files = generate_tumor_wig.wig
  }

  output {
    Array[File] params = ichorcna_call.params
    Array[File] seg = ichorcna_call.seg
    Array[File] genomewide_pdf = ichorcna_call.genomewide_pdf
    Array[File] allgenomewide_pdf = ichorcna_call.allgenomewide_pdf
    Array[File] correct_pdf = ichorcna_call.correct_pdf
    Array[File] rdata = ichorcna_call.rdata
    Array[File] wig = generate_tumor_wig.wig
    File validation_report = validate_outputs.report
  }
}

task generate_tumor_wig {
  meta {
    description: "Generate tumor WIG file from BED read count data"
    outputs: {
        wig: "WIG file created from binned read count data within input BED files"
    }
  }

  parameter_meta {
    counts_bed: "Tarball of per-chromosome BED files of read counts"
    sample_name: "Sample ID for naming output file"
    chrs: "Chromosomes to include in WIG file (R vector format)"
    genome_style: "Chromosome naming convention: NCBI or UCSC"
    window_size: "Window size in base pairs for WIG format"
    memory_gb: "Memory allocated for the task in GB"
    cpus: "Number of CPU cores allocated for the task"
  }

  input {
    File counts_bed
    String sample_name
    String chrs = "c(1:22, \"X\", \"Y\")"
    String genome_style = "UCSC"
    Int window_size = 500000
    Int memory_gb = 4
    Int cpus = 2
  }

  command <<<
    set -eo pipefail
    
    # Extract and sort BED files
    echo "Extracting and sorting BED files..."
    tar -xzf "~{counts_bed}" -O | sort -k 1,1V -k 2,2n > sorted_counts.bed
    echo "Sorted BED file contains $(wc -l < sorted_counts.bed) lines"

    # Parse chromosomes from the chrs parameter (e.g., "c(1:22, \"X\", \"Y\")")
    echo "Parsing expected chromosomes from: ~{chrs}"
    echo "~{chrs}" | \
    sed 's/c(//g' | sed 's/)//g' | \
    sed 's/1:22/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22/g' | \
    tr ',' '\n' | \
    sed 's/[" ]//g' | \
    grep -v '^$' > expected_chromosomes.txt

    echo "Expected chromosomes ($(wc -l < expected_chromosomes.txt) total):"
    cat expected_chromosomes.txt

    # Function to normalize chromosome names based on genome_style
    normalize_chr_name() {
        local chr=$1
        if [[ "~{genome_style}" == "UCSC" ]]; then
            if [[ "$chr" =~ ^[0-9]+$ ]] || [[ "$chr" == "X" ]] || [[ "$chr" == "Y" ]]; then
                echo "chr$chr"
            else
                echo "$chr"
            fi
        else
            # NCBI style - remove chr prefix if present
            echo "$chr" | sed 's/^chr//g'
        fi
    }

    # Create WIG file with proper headers for all expected chromosomes
    echo "Creating WIG file: ~{sample_name}.tumor.wig"
    > "~{sample_name}.tumor.wig"

    while IFS= read -r chr; do
        normalized_chr=$(normalize_chr_name "$chr")
        echo "Processing chromosome: $chr -> $normalized_chr"
        
        # Write WIG header for this chromosome
        echo "fixedStep chrom=$normalized_chr start=1 step=~{window_size} span=~{window_size}" >> "~{sample_name}.tumor.wig"
        
        # Extract data for this chromosome from sorted BED file
        # Column 4 contains the read counts
        chr_data=$(awk -v chr="$normalized_chr" '$1 == chr {print $4}' sorted_counts.bed)
        
        if [[ -n "$chr_data" ]]; then
            data_points=$(echo "$chr_data" | wc -l)
            echo "  Found $data_points data points for $normalized_chr"
            echo "$chr_data" >> "~{sample_name}.tumor.wig"
        else
            echo "  No data found for $normalized_chr - adding single zero value"
            echo "0" >> "~{sample_name}.tumor.wig"
        fi
        
    done < expected_chromosomes.txt

    # Clean up temporary files
    rm -f sorted_counts.bed expected_chromosomes.txt
  >>>

  output {
    File wig = "${sample_name}.tumor.wig"
  }

  runtime {
    docker: "getwilds/ichorcna:0.2.0"
    cpu: cpus
    memory: "~{memory_gb} GB"
  }
}

task ichorcna_call {
  meta {
    description: "Run ichorCNA tumor fraction estimation using pre-generated WIG file"
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
    tumor_wig: "WIG file containing tumor read count data"
    wig_gc: "GC-content WIG file"
    wig_map: "Mappability score WIG file"
    panel_of_norm_rds: "RDS file of median corrected depth from panel of normals"
    centromeres: "Text file containing Centromere locations"
    sample_name: "Sample ID"
    sex: "User-specified: male or female"
    genome: "Genome build (e.g. hg38)"
    genome_style: "Chromosome naming convention (use UCSC if desired output is to have 'chr' string): NCBI or UCSC"
    chrs: "Chromosomes to analyze (default: chr 1-22, X, and Y)"
    ploidy: "Tumor ploidy states to consider"
    normal_states: "Normal contamination states to consider"
    max_cn: "Maximum copy number to consider"
    memory_gb: "Memory allocated for the task in GB"
    cpus: "Number of CPU cores allocated for the task"
  }

  input {
    File tumor_wig
    File wig_gc
    File wig_map
    File panel_of_norm_rds
    File centromeres
    String sample_name
    String sex
    String genome
    String genome_style
    String chrs = "c(1:22, \"X\", \"Y\")"
    String ploidy = "c(2)"
    String normal_states = "c(0.1,0.5,.85)"
    Int max_cn = 3
    Int memory_gb = 16
    Int cpus = 6
  }

  command <<<
    set -eo pipefail

    # Run ichorCNA
    echo "Starting ichorCNA analysis..."
    Rscript /usr/local/bin/ichorCNA/scripts/runIchorCNA.R \
    --id "~{sample_name}" \
    --WIG "~{tumor_wig}" \
    --ploidy "~{ploidy}" \
    --normal "~{normal_states}" \
    --maxCN ~{max_cn} \
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
    --libdir /usr/local/bin/ichorCNA/

    # Filter segments based on sex
    echo "Filtering segments..."
    # Keep only biologically relevant segments by keeping all Y segments if
    # male and only non-Y segments if female
    awk -v G="~{sex}" '$2!~/Y/ || G=="male"' "~{sample_name}.seg.txt" \
    > "~{sample_name}.ichor.segs.txt"

    # Move PDF files to current directory
    echo "Organizing output files..."
    mv "~{sample_name}"/*.pdf .
  >>>

  output {
    File params = "${sample_name}.params.txt"
    File seg = "${sample_name}.ichor.segs.txt"
    File genomewide_pdf = "${sample_name}_genomeWide.pdf"
    File allgenomewide_pdf = "${sample_name}_genomeWide_all_sols.pdf"
    File correct_pdf = "${sample_name}_genomeWideCorrection.pdf"
    File rdata = "${sample_name}.RData"
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
