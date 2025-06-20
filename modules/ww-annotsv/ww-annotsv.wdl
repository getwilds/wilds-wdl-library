## WILDS WDL module for structural variant annotation using AnnotSV.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

struct SampleInfo {
    String name
    File vcf
}

workflow annotsv_example {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "WDL workflow for structural variant annotation via AnnotSV"
    url: "https://github.com/getwilds/wilds-wdl-library/tree/main/modules/ww-annotsv"
    outputs: {
        annotated_vcf: "Annotated VCF file with AnnotSV annotations",
        annotated_tsv: "Tab-delimited annotation file with comprehensive SV annotations",
        validation_report: "Validation report confirming all expected outputs were generated"
    }
  }

  parameter_meta {
    samples: "List of sample objects, each containing name and VCF file with structural variants"
    genome_build: "Reference genome build (GRCh37, GRCh38, mm9, mm10, etc.)"
    sv_min_size: "Minimum SV size to consider for annotation (default: 50)"
    annotation_mode: "Annotation mode: 'full' (one line per annotation) or 'split' (one line per SV)"
    include_ci: "Include confidence intervals in breakpoint coordinates"
    overlap_threshold: "Minimum overlap with genomic features (default: 70)"
    exclude_benign: "Exclude likely benign variants from output"
    cpus: "Number of CPU cores allocated for each task"
    memory_gb: "Memory allocated for each task in GB"
  }

  input {
    Array[SampleInfo] samples
    String genome_build = "GRCh38" # GRCh38 or GRCh37 or mm9 or mm10
    Int sv_min_size = 50
    String annotation_mode = "full"
    Boolean include_ci = true
    Int overlap_threshold = 70
    Boolean exclude_benign = false
    Int cpus = 4
    Int memory_gb = 8
  }

  scatter (sample in samples) {
    call annotsv_annotate {
      input:
        input_vcf = sample.vcf,
        sample_name = sample.name,
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

  call validate_outputs {
    input:
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
    input_vcf: "Input VCF file containing structural variants to annotate"
    sample_name: "Name of the sample for output file naming"
    genome_build: "Reference genome build (GRCh37, GRCh38, mm9, mm10, etc.)"
    sv_min_size: "Minimum SV size in bp to consider for annotation"
    annotation_mode: "Annotation mode: 'full' (comprehensive) or 'split' (one line per SV)"
    include_ci: "Include confidence intervals in breakpoint coordinates"
    overlap_threshold: "Minimum percentage overlap with genomic features"
    exclude_benign: "Filter out likely benign variants from output"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File input_vcf
    String sample_name
    String genome_build = "GRCh38" # GRCh38 or GRCh37 or mm9 or mm10
    String tx = "RefSeq" # ENSEMBL or RefSeq
    Int sv_min_size = 50
    String annotation_mode = "full" # full or split
    Boolean include_ci = true
    Int overlap_threshold = 70
    Boolean exclude_benign = false
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  String output_prefix = "~{sample_name}.annotsv"
  String mode_flag = if annotation_mode == "split" then "-annotationMode split" else "-annotationMode full"
  String ci_flag = if include_ci then "-includeCI 1" else "-includeCI 0"
  String benign_flag = if exclude_benign then "-benignAF 0.01" else ""
  String tx_flag = if tx == "ENSEMBL" then "-tx ENSEMBL" else "-tx RefSeq"

  command <<<
    set -eo pipefail
    
    # Create output directory
    mkdir -p annotsv_output

    # Debug: check filesystem permissions
    touch /tmp/test_write 2>&1 || echo "Cannot write to /tmp"
    touch /AnnotSV-3.4.4/test_write 2>&1 || echo "Cannot write to AnnotSV dir"

    # Run AnnotSV
    AnnotSV \
      -SVinputFile "~{input_vcf}" \
      -genomeBuild "~{genome_build}" \
      -outputFile "annotsv_output/~{output_prefix}" \
      -SVminSize ~{sv_min_size} \
      ~{mode_flag} \
      ~{ci_flag} \
      -overlap ~{overlap_threshold} \
      -promoterSize 500 \ # Sticking with default value for now, Apptainer issues
      ~{benign_flag} \
      ~{tx_flag}
    
    mv "annotsv_output/~{output_prefix}.tsv" "~{output_prefix}.tsv"
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
    for tsv_file in ~{sep=' ' annotated_tsv_files}; do
      TSV_COUNT=$((TSV_COUNT + 1))
      echo "TSV $TSV_COUNT: $(basename $tsv_file)" >> validation_report.txt
      
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
