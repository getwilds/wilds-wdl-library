## WILDS WDL module for structural variant annotation using AnnotSV.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

import "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/limit-test-runs/modules/ww-testdata/ww-testdata.wdl" as ww_testdata

workflow annotsv_example {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "WDL workflow for structural variant annotation via AnnotSV"
    url: "https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-annotsv"
    outputs: {
        annotated_tsv: "Tab-delimited annotation file with comprehensive SV annotations",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    vcfs: "Optional list of VCF files with structural variants. If not provided, test data will be used."
    genome_build: "Reference genome build (GRCh37, GRCh38)"
    annotation_mode: "Annotation mode: 'full' (one line per annotation) or 'split' (one line per SV)"
    include_ci: "Include confidence intervals in breakpoint coordinates"
    exclude_benign: "Exclude likely benign variants from output"
    sv_min_size: "Minimum SV size to consider for annotation (default: 50)"
    overlap_threshold: "Minimum overlap with genomic features (default: 70)"
    cpus: "Number of CPU cores allocated for each task"
    memory_gb: "Memory allocated for each task in GB"
  }

  input {
    Array[File]? vcfs
    String genome_build = "GRCh38" # GRCh38 or GRCh37
    String annotation_mode = "full"
    Boolean include_ci = true
    Boolean exclude_benign = false
    Int sv_min_size = 50
    Int overlap_threshold = 70
    Int cpus = 4
    Int memory_gb = 8
  }

  # If no VCFs provided, download test data
  if (!defined(vcfs)) {
    call ww_testdata.download_annotsv_vcf { }
  }

  # Determine which VCFs to use: provided ones or test data
  Array[File] vcfs_to_process = select_first([vcfs, [download_annotsv_vcf.test_vcf]])

  scatter (vcf in vcfs_to_process) {
    call annotsv_annotate { input:
        raw_vcf = vcf,
        genome_build = genome_build,
        sv_min_size = sv_min_size,
        annotation_mode = annotation_mode,
        include_ci = include_ci,
        overlap_threshold = overlap_threshold,
        exclude_benign = exclude_benign,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  call validate_outputs { input:
      annotated_tsv_files = annotsv_annotate.annotated_tsv
  }

  output {
    Array[File] annotated_tsv = annotsv_annotate.annotated_tsv
    File validation_report = validate_outputs.report
  }
}

task annotsv_annotate {
  meta {
    description: "Annotate structural variants using AnnotSV with comprehensive genomic and clinical annotations"
    outputs: {
        annotated_tsv: "Tab-delimited file with detailed annotations per SV"
    }
  }

  parameter_meta {
    raw_vcf: "Input VCF file containing structural variants to annotate"
    genome_build: "Reference genome build (GRCh37, GRCh38)"
    tx_source: "Transcript annotation source (ENSEMBL, RefSeq)"
    annotation_mode: "Annotation mode: 'full' (comprehensive) or 'split' (one line per SV)"
    include_ci: "Include confidence intervals in breakpoint coordinates"
    exclude_benign: "Filter out likely benign variants from output"
    sv_min_size: "Minimum SV size in bp to consider for annotation"
    overlap_threshold: "Minimum percentage overlap with genomic features"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File raw_vcf
    String genome_build = "GRCh38" # GRCh38 or GRCh37
    String tx_source = "RefSeq" # ENSEMBL or RefSeq
    String annotation_mode = "full" # full or split
    Boolean include_ci = true
    Boolean exclude_benign = false
    Int sv_min_size = 50
    Int overlap_threshold = 70
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  String sample_name = basename(basename(raw_vcf, ".gz"), ".vcf")
  String output_prefix = "~{sample_name}.annotsv"
  String mode_flag = if annotation_mode == "split" then "-annotationMode split" else "-annotationMode full"
  String ci_flag = if include_ci then "-includeCI 1" else "-includeCI 0"
  String benign_flag = if exclude_benign then "-benignAF 0.01" else ""
  String tx_flag = if tx_source == "ENSEMBL" then "-tx ENSEMBL" else "-tx RefSeq"
  # Sticking with default promoter size for now
  # Can't write the necessary bedfiles for custom values in the Apptainer container

  command <<<
    set -eo pipefail
    
    # Create output directory
    mkdir -p annotsv_output

    # Run AnnotSV
    AnnotSV \
      -SVinputFile "~{raw_vcf}" \
      -genomeBuild "~{genome_build}" \
      -outputFile "annotsv_output/~{output_prefix}" \
      -SVminSize ~{sv_min_size} \
      ~{mode_flag} \
      ~{ci_flag} \
      -overlap ~{overlap_threshold} \
      -promoterSize 500 \
      ~{benign_flag} \
      ~{tx_flag}
    
    # Check if there were any SVs to annotate
    if grep -q "no SV to annotate" stdout; then
      echo "No SV's found in the input VCF. Creating empty output file."
      touch "~{output_prefix}.tsv"
    else
      mv "annotsv_output/~{output_prefix}.tsv" "~{output_prefix}.tsv"
    fi
  >>>

  output {
    File annotated_tsv = "~{output_prefix}.tsv"
  }

  runtime {
    docker: "getwilds/annotsv:3.4.4"
    memory: "~{memory_gb}GB"
    cpu: cpu_cores
  }
}

task validate_outputs {
  meta {
    description: "Validate AnnotSV outputs and generate comprehensive statistics"
    outputs: {
        report: "Validation summary with structural variant annotation statistics"
    }
  }

  parameter_meta {
    annotated_tsv_files: "Array of annotated TSV files to validate"
  }

  input {
    Array[File] annotated_tsv_files
  }

  command <<<
    set -eo pipefail
    
    echo "AnnotSV Validation Report" > validation_report.txt
    echo "=========================" >> validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt
    
    # Validate TSV files
    echo "TSV File Validation:" >> validation_report.txt
    echo "-------------------" >> validation_report.txt
    
    TSV_COUNT=0
    TOTAL_ANNOTATIONS=0
    for tsv_file in ~{sep=" " annotated_tsv_files}; do
      TSV_COUNT=$((TSV_COUNT + 1))
      echo "TSV $TSV_COUNT: $(basename $tsv_file)" >> validation_report.txt
      
      # Count annotations in the TSV file
      ANNOTATION_COUNT=$(wc -l < "$tsv_file")
      TOTAL_ANNOTATIONS=$((TOTAL_ANNOTATIONS + ANNOTATION_COUNT))

      # Check if file exists and is not empty
      if [ -f "$tsv_file" ] && [ -s "$tsv_file" ]; then
        echo "  File exists and is not empty" >> validation_report.txt
      else
        echo "  File missing or empty" >> validation_report.txt
      fi
    done
    
    echo "" >> validation_report.txt
    echo "Summary Statistics:" >> validation_report.txt
    echo "------------------" >> validation_report.txt
    echo "Total TSV files processed: $TSV_COUNT" >> validation_report.txt
    echo "Total annotations generated: $TOTAL_ANNOTATIONS" >> validation_report.txt
    
    echo "" >> validation_report.txt
    echo "Validation completed successfully." >> validation_report.txt
  >>>

  output {
    File report = "validation_report.txt"
  }

  runtime {
    docker: "getwilds/annotsv:3.4.4"
    memory: "4GB"
    cpu: 1
  }
}
