## WILDS WDL for RNA-seq quality control analysis using RSeQC.
## Designed to be a modular component within the WILDS ecosystem that can be used
## independently or integrated with other WILDS workflows.

version 1.0

task run_rseqc {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Run comprehensive RSeQC quality control metrics on aligned RNA-seq data"
    outputs: {
        read_distribution: "Distribution of reads across genomic features (CDS, UTR, introns, intergenic)",
        gene_body_coverage: "Gene body coverage plot showing 5' to 3' bias",
        infer_experiment: "Results of strand specificity inference",
        bam_stat: "Basic alignment statistics from the BAM file",
        junction_annotation: "Splice junction annotation and classification",
        rseqc_summary: "Summary report of all RSeQC metrics"
    }
  }

  parameter_meta {
    bam_file: "Aligned reads in BAM format"
    bam_index: "Index file for the aligned BAM file"
    ref_bed: "Reference genome annotation in BED format (12-column)"
    sample_name: "Sample name for output files"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bam_file
    File bam_index
    File ref_bed
    String sample_name
    Int cpu_cores = 2
    Int memory_gb = 4
  }

  command <<<
    set -eo pipefail

    echo "Running RSeQC quality control analyses..."

    # Read distribution across genomic features
    echo "1. Analyzing read distribution..."
    read_distribution.py -i "~{bam_file}" -r "~{ref_bed}" > "~{sample_name}.read_distribution.txt"

    # Gene body coverage (5' to 3' bias)
    echo "2. Analyzing gene body coverage..."
    geneBody_coverage.py -i "~{bam_file}" -r "~{ref_bed}" -o "~{sample_name}"

    # Infer experiment (strand specificity)
    echo "3. Inferring strand specificity..."
    infer_experiment.py -i "~{bam_file}" -r "~{ref_bed}" > "~{sample_name}.infer_experiment.txt"

    # Basic BAM statistics
    echo "4. Computing BAM statistics..."
    bam_stat.py -i "~{bam_file}" > "~{sample_name}.bam_stat.txt"

    # Junction annotation
    echo "5. Annotating splice junctions..."
    junction_annotation.py -i "~{bam_file}" -r "~{ref_bed}" -o "~{sample_name}" 2>&1 | tee "~{sample_name}.junction_annotation.log"

    # Create summary report
    echo "6. Creating summary report..."
    cat > "~{sample_name}.rseqc_summary.txt" << EOF
RSeQC Quality Control Summary
Sample: ~{sample_name}
Date: $(date)

================================================================================
ANALYSES PERFORMED:
================================================================================
1. Read Distribution - Distribution of reads across genomic features
2. Gene Body Coverage - 5' to 3' coverage bias assessment
3. Infer Experiment - Strand specificity determination
4. BAM Statistics - Basic alignment statistics
5. Junction Annotation - Splice junction classification

================================================================================
OUTPUT FILES:
================================================================================
- ~{sample_name}.read_distribution.txt
- ~{sample_name}.geneBodyCoverage.txt
- ~{sample_name}.geneBodyCoverage.curves.pdf
- ~{sample_name}.infer_experiment.txt
- ~{sample_name}.bam_stat.txt
- ~{sample_name}.junction.bed
- ~{sample_name}.junction_annotation.log
- ~{sample_name}.splice_events.pdf
- ~{sample_name}.splice_junction.pdf

For detailed interpretation, please refer to RSeQC documentation:
http://rseqc.sourceforge.net/
EOF

    echo "RSeQC analysis complete!"
  >>>

  output {
    File read_distribution = "~{sample_name}.read_distribution.txt"
    File gene_body_coverage = "~{sample_name}.geneBodyCoverage.txt"
    File gene_body_coverage_plot = "~{sample_name}.geneBodyCoverage.curves.pdf"
    File infer_experiment = "~{sample_name}.infer_experiment.txt"
    File bam_stat = "~{sample_name}.bam_stat.txt"
    File junction_bed = "~{sample_name}.junction.bed"
    File junction_log = "~{sample_name}.junction_annotation.log"
    File splice_events_plot = "~{sample_name}.splice_events.pdf"
    File splice_junction_plot = "~{sample_name}.splice_junction.pdf"
    File rseqc_summary = "~{sample_name}.rseqc_summary.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/rseqc:5.0.1--py39hf95cd2a_1"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task gtf_to_bed {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Convert GTF annotation file to BED12 format for RSeQC compatibility"
    outputs: {
        bed_file: "BED12-formatted annotation file ready for use with RSeQC"
    }
  }

  parameter_meta {
    gtf_file: "GTF annotation file to convert"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File gtf_file
    Int cpu_cores = 1
    Int memory_gb = 2
  }

  command <<<
    set -eo pipefail

    echo "Converting GTF to BED12 format..."
    # Convert GTF to BED12 format
    # Extract transcript features and convert to BED12 format
    awk 'BEGIN{OFS="\t"}
    $3=="transcript" {
      # Extract chromosome, start, end, strand
      chr=$1
      start=$4-1  # Convert to 0-based
      end=$5
      strand=$7

      # Extract gene_id and transcript_id from attributes using simpler string functions
      # Look for gene_id "value" pattern
      gene_id = ""
      trans_id = ""

      # Split the attributes field and extract IDs
      n = split($0, parts, /gene_id "/)
      if (n > 1) {
        split(parts[2], gene_parts, /"/)
        gene_id = gene_parts[1]
      }

      n = split($0, parts, /transcript_id "/)
      if (n > 1) {
        split(parts[2], trans_parts, /"/)
        trans_id = trans_parts[1]
      }

      # Print BED12 format (simplified - using transcript as single exon)
      # chr, start, end, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
      if (trans_id != "") {
        print chr, start, end, trans_id, 0, strand, start, end, "0,0,0", 1, end-start, 0
      }
    }' "~{gtf_file}" > annotation.bed
  >>>

  output {
    File bed_file = "annotation.bed"
  }

  runtime {
    docker: "ubuntu:20.04"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
