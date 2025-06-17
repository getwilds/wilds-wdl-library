## WILDS WDL module for structural variant annotation using AnnotSV.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

struct SampleInfo {
    String name
    File vcf
    File? vcf_index
}

struct RefGenome {
    String name
    File? fasta
    File? fasta_index
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
    genome_build: "Reference genome build (hg19, hg38, mm9, mm10, etc.)"
    sv_min_size: "Minimum SV size to consider for annotation (default: 50)"
    annotation_mode: "Annotation mode: 'full' (one line per annotation) or 'split' (one line per SV)"
    include_ci: "Include confidence intervals in breakpoint coordinates"
    overlap_threshold: "Minimum overlap with genomic features (default: 70)"
    promoter_size: "Size of promoter regions in bp (default: 1000)"
    exclude_benign: "Exclude likely benign variants from output"
    custom_bed_files: "Optional array of custom BED files for additional annotations"
    cpus: "Number of CPU cores allocated for each task"
    memory_gb: "Memory allocated for each task in GB"
  }

  input {
    Array[SampleInfo] samples
    String genome_build = "hg38"
    Int sv_min_size = 50
    String annotation_mode = "full"
    Boolean include_ci = true
    Float overlap_threshold = 70.0
    Int promoter_size = 1000
    Boolean exclude_benign = false
    Array[File]? custom_bed_files
    Int cpus = 4
    Int memory_gb = 8
  }

  scatter (sample in samples) {
    call annotsv_annotate {
      input:
        input_vcf = sample.vcf,
        input_vcf_index = sample.vcf_index,
        sample_name = sample.name,
        genome_build = genome_build,
        sv_min_size = sv_min_size,
        annotation_mode = annotation_mode,
        include_ci = include_ci,
        overlap_threshold = overlap_threshold,
        promoter_size = promoter_size,
        exclude_benign = exclude_benign,
        custom_bed_files = custom_bed_files,
        cpu_cores = cpus,
        memory_gb = memory_gb
    }
  }

  call validate_outputs {
    input:
      annotated_vcf_files = annotsv_annotate.annotated_vcf,
      annotated_tsv_files = annotsv_annotate.annotated_tsv
  }

  output {
    Array[File] annotated_vcf = annotsv_annotate.annotated_vcf
    Array[File] annotated_tsv = annotsv_annotate.annotated_tsv
    File validation_report = validate_outputs.report
  }
}

task annotsv_annotate {
  meta {
    description: "Annotate structural variants using AnnotSV with comprehensive genomic and clinical annotations"
    outputs: {
        annotated_vcf: "VCF file with AnnotSV annotations added",
        annotated_tsv: "Tab-delimited file with detailed annotations per SV"
    }
  }

  parameter_meta {
    input_vcf: "Input VCF file containing structural variants to annotate"
    input_vcf_index: "Index file for the input VCF (optional but recommended)"
    sample_name: "Name of the sample for output file naming"
    genome_build: "Reference genome build (hg19, hg38, mm9, mm10, etc.)"
    sv_min_size: "Minimum SV size in bp to consider for annotation"
    annotation_mode: "Annotation mode: 'full' (comprehensive) or 'split' (one line per SV)"
    include_ci: "Include confidence intervals in breakpoint coordinates"
    overlap_threshold: "Minimum percentage overlap with genomic features"
    promoter_size: "Size of promoter regions to consider in base pairs"
    exclude_benign: "Filter out likely benign variants from output"
    custom_bed_files: "Array of custom BED files for additional region-based annotations"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File input_vcf
    File? input_vcf_index
    String sample_name
    String genome_build = "hg38"
    Int sv_min_size = 50
    String annotation_mode = "full"
    Boolean include_ci = true
    Float overlap_threshold = 70.0
    Int promoter_size = 1000
    Boolean exclude_benign = false
    Array[File]? custom_bed_files
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  String output_prefix = "~{sample_name}.annotsv"
  String mode_flag = if annotation_mode == "split" then "-annotationMode split" else "-annotationMode full"
  String ci_flag = if include_ci then "-includeCI yes" else "-includeCI no"
  String benign_flag = if exclude_benign then "-benignAF 0.01" else ""

  command <<<
    set -eo pipefail
    
    # Create output directory
    mkdir -p annotsv_output
    
    # Handle custom BED files if provided
    CUSTOM_BED_ARGS=""
    if [ ! -z "~{sep=' ' custom_bed_files}" ]; then
      mkdir -p custom_annotations
      count=1
      for bedfile in ~{sep=' ' custom_bed_files}; do
        cp "$bedfile" "custom_annotations/custom_${count}.bed"
        CUSTOM_BED_ARGS="$CUSTOM_BED_ARGS -bedFile custom_annotations/custom_${count}.bed"
        count=$((count + 1))
      done
    fi
    
    # Run AnnotSV
    AnnotSV \
      -SVinputFile "~{input_vcf}" \
      -genomeBuild "~{genome_build}" \
      -outputFile "annotsv_output/~{output_prefix}" \
      -SVminSize ~{sv_min_size} \
      ~{mode_flag} \
      ~{ci_flag} \
      -overlap ~{overlap_threshold} \
      -promoterSize ~{promoter_size} \
      ~{benign_flag} \
      $CUSTOM_BED_ARGS \
      -tx ENSEMBL \
      -typeOfAnnotation both
    
    # AnnotSV outputs both annotated VCF and TSV files
    # Move outputs to expected locations
    if [ -f "annotsv_output/~{output_prefix}.annotated.vcf" ]; then
      mv "annotsv_output/~{output_prefix}.annotated.vcf" "~{output_prefix}.annotated.vcf"
    fi
    
    if [ -f "annotsv_output/~{output_prefix}.annotated.tsv" ]; then
      mv "annotsv_output/~{output_prefix}.annotated.tsv" "~{output_prefix}.annotated.tsv"
    fi
    
    # Compress and index VCF if it exists
    if [ -f "~{output_prefix}.annotated.vcf" ]; then
      bgzip "~{output_prefix}.annotated.vcf"
      tabix -p vcf "~{output_prefix}.annotated.vcf.gz"
    fi
    
    # Ensure TSV output exists (AnnotSV always creates TSV)
    if [ ! -f "~{output_prefix}.annotated.tsv" ]; then
      echo "Error: AnnotSV did not produce expected TSV output"
      exit 1
    fi
    
    # Create summary statistics
    echo "AnnotSV annotation completed for sample: ~{sample_name}" > annotation_summary.txt
    echo "Genome build: ~{genome_build}" >> annotation_summary.txt
    echo "Minimum SV size: ~{sv_min_size} bp" >> annotation_summary.txt
    echo "Annotation mode: ~{annotation_mode}" >> annotation_summary.txt
    
    # Count total variants annotated
    TOTAL_VARIANTS=$(tail -n +2 "~{output_prefix}.annotated.tsv" | wc -l)
    echo "Total variants annotated: $TOTAL_VARIANTS" >> annotation_summary.txt
    
    # Count variants by type if possible
    if [ -f "~{output_prefix}.annotated.tsv" ]; then
      echo "Annotation summary:" >> annotation_summary.txt
      tail -n +2 "~{output_prefix}.annotated.tsv" | cut -f4 | sort | uniq -c >> annotation_summary.txt
    fi
  >>>

  output {
    File annotated_vcf = "~{output_prefix}.annotated.vcf.gz"
    File annotated_tsv = "~{output_prefix}.annotated.tsv"
    File summary = "annotation_summary.txt"
  }

  runtime {
    docker: "getwilds/annotsv:3.4.6"
    memory: "~{memory_gb}GB"
    cpu: cpu_cores
    disks: "local-disk 50 SSD"
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
    annotated_vcf_files: "Array of annotated VCF files to validate"
    annotated_tsv_files: "Array of annotated TSV files to validate"
  }

  input {
    Array[File] annotated_vcf_files
    Array[File] annotated_tsv_files
  }

  command <<<
    set -eo pipefail
    
    echo "AnnotSV Validation Report" > validation_report.txt
    echo "=========================" >> validation_report.txt
    echo "Generated on: $(date)" >> validation_report.txt
    echo "" >> validation_report.txt
    
    # Validate VCF files
    echo "VCF File Validation:" >> validation_report.txt
    echo "-------------------" >> validation_report.txt
    
    VCF_COUNT=0
    for vcf_file in ~{sep=' ' annotated_vcf_files}; do
      VCF_COUNT=$((VCF_COUNT + 1))
      echo "VCF $VCF_COUNT: $(basename $vcf_file)" >> validation_report.txt
      
      # Check if file exists and is not empty
      if [ -f "$vcf_file" ] && [ -s "$vcf_file" ]; then
        echo "  ✓ File exists and is not empty" >> validation_report.txt
        
        # Count variants in VCF
        if command -v zcat >/dev/null 2>&1; then
          VARIANT_COUNT=$(zcat "$vcf_file" | grep -v "^#" | wc -l)
          echo "  ✓ Contains $VARIANT_COUNT variants" >> validation_report.txt
        fi
      else
        echo "  ✗ File missing or empty" >> validation_report.txt
      fi
    done
    
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
        echo "  ✓ File exists and is not empty" >> validation_report.txt
        
        # Count annotations
        ANNOTATION_COUNT=$(tail -n +2 "$tsv_file" | wc -l)
        TOTAL_ANNOTATIONS=$((TOTAL_ANNOTATIONS + ANNOTATION_COUNT))
        echo "  ✓ Contains $ANNOTATION_COUNT annotations" >> validation_report.txt
        
        # Check for key AnnotSV columns
        HEADER=$(head -1 "$tsv_file")
        if echo "$HEADER" | grep -q "AnnotSV_ID"; then
          echo "  ✓ Contains AnnotSV_ID column" >> validation_report.txt
        fi
        if echo "$HEADER" | grep -q "SV_chrom"; then
          echo "  ✓ Contains SV_chrom column" >> validation_report.txt
        fi
        if echo "$HEADER" | grep -q "SV_type"; then
          echo "  ✓ Contains SV_type column" >> validation_report.txt
        fi
        if echo "$HEADER" | grep -q "AnnotSV_ranking_score"; then
          echo "  ✓ Contains AnnotSV ranking score" >> validation_report.txt
        fi
      else
        echo "  ✗ File missing or empty" >> validation_report.txt
      fi
    done
    
    echo "" >> validation_report.txt
    echo "Summary Statistics:" >> validation_report.txt
    echo "------------------" >> validation_report.txt
    echo "Total VCF files processed: $VCF_COUNT" >> validation_report.txt
    echo "Total TSV files processed: $TSV_COUNT" >> validation_report.txt
    echo "Total annotations generated: $TOTAL_ANNOTATIONS" >> validation_report.txt
    
    # Check if all files were processed successfully
    if [ $VCF_COUNT -eq ${#annotated_vcf_files[@]} ] && [ $TSV_COUNT -eq ${#annotated_tsv_files[@]} ]; then
      echo "✓ All expected output files were generated successfully" >> validation_report.txt
    else
      echo "✗ Some expected output files were missing" >> validation_report.txt
    fi
    
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
    disks: "local-disk 10 SSD"
  }
}
