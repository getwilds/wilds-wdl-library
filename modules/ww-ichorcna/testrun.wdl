version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-testdata/ww-testdata.wdl" as ww_testdata
import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-ichorcna/ww-ichorcna.wdl" as ww_ichorcna

struct IchorSample {
    String name
    File bam
    File bam_index
}

workflow ichorcna_example {
  # Download ichorCNA reference data
  call ww_testdata.download_ichor_data { }

  # Download test sample data
  call ww_testdata.download_bam_data { }

  # Create test samples array
  Array[IchorSample] final_samples = [
    {
      "name": "demo_sample",
      "bam": download_bam_data.bam,
      "bam_index": download_bam_data.bai
    }
  ]

  scatter (sample in final_samples) {
    call ww_ichorcna.readcounter_wig { input:
      bam_file = sample.bam,
      bam_index = sample.bam_index,
      sample_name = sample.name,
      window_size = 500000,
      chromosomes = ["chr1"],
      memory_gb = 8,
      cpus = 2
    }

    call ww_ichorcna.ichorcna_call { input:
      wig_tumor = readcounter_wig.wig_file,
      wig_gc = download_ichor_data.wig_gc,
      wig_map = download_ichor_data.wig_map,
      panel_of_norm_rds = download_ichor_data.panel_of_norm_rds,
      centromeres = download_ichor_data.centromeres,
      name = sample.name,
      sex = "male",
      chrs = "c(1)",
      genome = "hg38",
      genome_style = "UCSC",
      cpus = 2,
      memory_gb = 8
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

task validate_outputs {
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
