## WILDS WDL module for structural variant annotation using AnnotSV.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task annotsv_annotate {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
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
