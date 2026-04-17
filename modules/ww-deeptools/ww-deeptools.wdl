## WILDS WDL module for deepTools - tools for exploring deep sequencing data.
## deepTools provides a suite of tools for analyzing and visualizing high-throughput
## sequencing data, including coverage track generation, signal comparison, and
## heatmap/profile visualization of genomic regions.

version 1.0

#### TASK DEFINITIONS ####

task bam_coverage {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Generates a normalized coverage track (bigWig or bedGraph) from a BAM file using deepTools bamCoverage"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deeptools/ww-deeptools.wdl"
    outputs: {
        coverage_file: "Normalized coverage track file in bigWig or bedGraph format"
    }
    topic: "genomics,transcriptomics,epigenomics,ribosome_profiling"
    species: "any"
    operation: "quantification"
    in_sample_req: "bam:nucleic_acid_sequence_alignment:bam,bai:data_index:bai"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "coverage_file:annotation_track:bigwig|bedgraph"
    out_ref: "none"
  }

  parameter_meta {
    bam: "Input BAM file (must be sorted and indexed)"
    bai: "BAM index file (.bai)"
    sample_name: "Name identifier for the sample, used as output file prefix"
    output_format: "Output file format: bigwig or bedgraph"
    bin_size: "Length in bases of the output bins (default: 50)"
    normalize_using: "Normalization method: RPKM, CPM, BPM, RPGC, or None"
    effective_genome_size: "Effective genome size for RPGC normalization (e.g., 2913022398 for hg38)"
    extend_reads: "Extend each read to the given fragment length (useful for single-end data)"
    ignore_duplicates: "If true, reads flagged as duplicates are ignored"
    min_mapping_quality: "Minimum mapping quality for reads to be considered"
    extra_args: "Additional arguments to pass to bamCoverage"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File bam
    File bai
    String sample_name
    String output_format = "bigwig"
    Int bin_size = 50
    String normalize_using = "RPKM"
    Int? effective_genome_size
    Int? extend_reads
    Boolean ignore_duplicates = false
    Int min_mapping_quality = 0
    String extra_args = ""
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  String ext = if output_format == "bigwig" then ".bw" else ".bedgraph"

  command <<<
    set -eo pipefail

    bamCoverage \
      --bam ~{bam} \
      --outFileName "~{sample_name}~{ext}" \
      --outFileFormat ~{output_format} \
      --binSize ~{bin_size} \
      --normalizeUsing ~{normalize_using} \
      ~{if defined(effective_genome_size) then "--effectiveGenomeSize " + effective_genome_size else ""} \
      ~{if defined(extend_reads) then "--extendReads " + extend_reads else ""} \
      ~{if ignore_duplicates then "--ignoreDuplicates" else ""} \
      --minMappingQuality ~{min_mapping_quality} \
      --numberOfProcessors ~{cpu_cores} \
      ~{extra_args}
  >>>

  output {
    File coverage_file = "~{sample_name}~{ext}"
  }

  runtime {
    docker: "getwilds/deeptools:3.5.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task bam_compare {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Compares two BAM files and generates a coverage comparison track using deepTools bamCompare"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deeptools/ww-deeptools.wdl"
    outputs: {
        comparison_file: "Coverage comparison track in bigWig or bedGraph format"
    }
    topic: "genomics,transcriptomics,epigenomics,ribosome_profiling"
    species: "any"
    operation: "statistical_calculation"
    in_sample_req: "treatment_bam:nucleic_acid_sequence_alignment:bam,treatment_bai:data_index:bai,control_bam:nucleic_acid_sequence_alignment:bam,control_bai:data_index:bai"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "comparison_file:annotation_track:bigwig|bedgraph"
    out_ref: "none"
  }

  parameter_meta {
    treatment_bam: "Treatment/ChIP BAM file (sorted and indexed)"
    treatment_bai: "Treatment BAM index file (.bai)"
    control_bam: "Control/input BAM file (sorted and indexed)"
    control_bai: "Control BAM index file (.bai)"
    sample_name: "Name identifier for the comparison output"
    output_format: "Output file format: bigwig or bedgraph"
    bin_size: "Length in bases of the output bins (default: 50)"
    operation: "Operation to perform: log2, ratio, subtract, add, mean, reciprocal_ratio, first, second"
    normalize_using: "Normalization method: RPKM, CPM, BPM, RPGC, or None"
    scale_factors_method: "Method to compute scale factors: None or readCount. Must be None when normalizeUsing is not None"
    extra_args: "Additional arguments to pass to bamCompare"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File treatment_bam
    File treatment_bai
    File control_bam
    File control_bai
    String sample_name
    String output_format = "bigwig"
    Int bin_size = 50
    String operation = "log2"
    String normalize_using = "RPKM"
    String scale_factors_method = "None"
    String extra_args = ""
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  String ext = if output_format == "bigwig" then ".bw" else ".bedgraph"

  command <<<
    set -eo pipefail

    bamCompare \
      --bamfile1 ~{treatment_bam} \
      --bamfile2 ~{control_bam} \
      --outFileName "~{sample_name}_compare~{ext}" \
      --outFileFormat ~{output_format} \
      --binSize ~{bin_size} \
      --operation ~{operation} \
      --normalizeUsing ~{normalize_using} \
      --scaleFactorsMethod ~{scale_factors_method} \
      --numberOfProcessors ~{cpu_cores} \
      ~{extra_args}
  >>>

  output {
    File comparison_file = "~{sample_name}_compare~{ext}"
  }

  runtime {
    docker: "getwilds/deeptools:3.5.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task compute_matrix {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Computes a matrix of signal values from bigWig files over genomic regions using deepTools computeMatrix"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deeptools/ww-deeptools.wdl"
    outputs: {
        matrix_gz: "Compressed matrix file for use with plotHeatmap and plotProfile",
        matrix_tab: "Tab-separated matrix values file"
    }
    topic: "genomics,transcriptomics,epigenomics,ribosome_profiling"
    species: "any"
    operation: "statistical_calculation"
    in_sample_req: "bigwig_files:annotation_track:bigwig,regions_file:annotation_track:bed|gtf"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "matrix_gz:sequence_report:matrix,matrix_tab:sequence_report:tsv"
    out_ref: "none"
  }

  parameter_meta {
    bigwig_files: "Array of bigWig score files to compute the matrix from"
    regions_file: "BED or GTF file defining the genomic regions of interest"
    sample_name: "Name identifier for the output files"
    mode: "Computation mode: reference-point or scale-regions"
    reference_point: "Reference point for reference-point mode: TSS, TES, or center"
    before_region: "Distance upstream of the reference point or start of region (bp)"
    after_region: "Distance downstream of the reference point or end of region (bp)"
    bin_size: "Bin size for averaging signal values (default: 10)"
    extra_args: "Additional arguments to pass to computeMatrix"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    Array[File] bigwig_files
    File regions_file
    String sample_name
    String mode = "reference-point"
    String reference_point = "TSS"
    Int before_region = 3000
    Int after_region = 3000
    Int bin_size = 10
    String extra_args = ""
    Int cpu_cores = 4
    Int memory_gb = 16
  }

  Boolean is_reference_point = mode == "reference-point"
  String ref_point_arg = if is_reference_point then "--referencePoint " + reference_point else ""

  command <<<
    set -eo pipefail

    computeMatrix ~{mode} \
      --scoreFileName ~{sep=" " bigwig_files} \
      --regionsFileName ~{regions_file} \
      --outFileName "~{sample_name}_matrix.gz" \
      --outFileNameMatrix "~{sample_name}_matrix.tab" \
      ~{ref_point_arg} \
      --beforeRegionStartLength ~{before_region} \
      --afterRegionStartLength ~{after_region} \
      --binSize ~{bin_size} \
      --numberOfProcessors ~{cpu_cores} \
      ~{extra_args}
  >>>

  output {
    File matrix_gz = "~{sample_name}_matrix.gz"
    File matrix_tab = "~{sample_name}_matrix.tab"
  }

  runtime {
    docker: "getwilds/deeptools:3.5.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task plot_heatmap {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Creates a heatmap visualization of signal enrichment over genomic regions using deepTools plotHeatmap"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deeptools/ww-deeptools.wdl"
    outputs: {
        heatmap: "Heatmap image file in the specified format"
    }
    topic: "genomics,transcriptomics,epigenomics,ribosome_profiling"
    species: "any"
    operation: "visualisation"
    in_sample_req: "matrix_gz:sequence_report:matrix"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "heatmap:plot:png|pdf|svg|eps"
    out_ref: "none"
  }

  parameter_meta {
    matrix_gz: "Compressed matrix file from computeMatrix"
    sample_name: "Name identifier for the output file"
    plot_format: "Output image format: png, pdf, svg, or eps"
    color_map: "Matplotlib colormap name (e.g., RdYlBu, viridis, Blues)"
    extra_args: "Additional arguments to pass to plotHeatmap"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File matrix_gz
    String sample_name
    String plot_format = "png"
    String color_map = "RdYlBu"
    String extra_args = ""
    Int cpu_cores = 1
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    plotHeatmap \
      --matrixFile ~{matrix_gz} \
      --outFileName "~{sample_name}_heatmap.~{plot_format}" \
      --colorMap ~{color_map} \
      ~{extra_args}
  >>>

  output {
    File heatmap = "~{sample_name}_heatmap.~{plot_format}"
  }

  runtime {
    docker: "getwilds/deeptools:3.5.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task plot_profile {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Creates a profile plot of average signal enrichment over genomic regions using deepTools plotProfile"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deeptools/ww-deeptools.wdl"
    outputs: {
        profile: "Profile plot image file in the specified format"
    }
    topic: "genomics,transcriptomics,epigenomics,ribosome_profiling"
    species: "any"
    operation: "visualisation"
    in_sample_req: "matrix_gz:sequence_report:matrix"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "profile:plot:png|pdf|svg|eps"
    out_ref: "none"
  }

  parameter_meta {
    matrix_gz: "Compressed matrix file from computeMatrix"
    sample_name: "Name identifier for the output file"
    plot_format: "Output image format: png, pdf, svg, or eps"
    extra_args: "Additional arguments to pass to plotProfile"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File matrix_gz
    String sample_name
    String plot_format = "png"
    String extra_args = ""
    Int cpu_cores = 1
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    plotProfile \
      --matrixFile ~{matrix_gz} \
      --outFileName "~{sample_name}_profile.~{plot_format}" \
      ~{extra_args}
  >>>

  output {
    File profile = "~{sample_name}_profile.~{plot_format}"
  }

  runtime {
    docker: "getwilds/deeptools:3.5.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task multi_bam_summary {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Computes read coverage summary across multiple BAM files using deepTools multiBamSummary"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deeptools/ww-deeptools.wdl"
    outputs: {
        summary_npz: "Compressed numpy array with coverage data for correlation and PCA analysis",
        raw_counts: "Tab-separated file with raw read counts per bin"
    }
    topic: "genomics,transcriptomics,epigenomics,ribosome_profiling"
    species: "any"
    operation: "quantification,sequencing_quality_control"
    in_sample_req: "bam_files:nucleic_acid_sequence_alignment:bam,bai_files:data_index:bai"
    in_sample_opt: "regions_file:annotation_track:bed"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "summary_npz:sequence_report:npz,raw_counts:sequence_report:tsv"
    out_ref: "none"
  }

  parameter_meta {
    bam_files: "Array of sorted and indexed BAM files to summarize"
    bai_files: "Array of BAM index files (.bai) corresponding to the BAM files"
    sample_name: "Name identifier for the output files"
    mode: "Computation mode: bins (genome-wide) or BED-file (specific regions)"
    bin_size: "Bin size in bp for bins mode (default: 10000)"
    regions_file: "Optional BED file for BED-file mode"
    extra_args: "Additional arguments to pass to multiBamSummary"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    Array[File] bam_files
    Array[File] bai_files
    String sample_name
    String mode = "bins"
    Int bin_size = 10000
    File? regions_file
    String extra_args = ""
    Int cpu_cores = 4
    Int memory_gb = 16
  }

  Boolean is_bins_mode = mode == "bins"

  command <<<
    set -eo pipefail

    multiBamSummary ~{mode} \
      --bamfiles ~{sep=" " bam_files} \
      --outFileName "~{sample_name}_summary.npz" \
      --outRawCounts "~{sample_name}_raw_counts.tab" \
      ~{if is_bins_mode then "--binSize " + bin_size else ""} \
      ~{if defined(regions_file) then "--BED " + regions_file else ""} \
      --numberOfProcessors ~{cpu_cores} \
      ~{extra_args}
  >>>

  output {
    File summary_npz = "~{sample_name}_summary.npz"
    File raw_counts = "~{sample_name}_raw_counts.tab"
  }

  runtime {
    docker: "getwilds/deeptools:3.5.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task plot_correlation {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Generates a correlation heatmap or scatterplot from multiBamSummary or multiBigwigSummary output using deepTools plotCorrelation"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deeptools/ww-deeptools.wdl"
    outputs: {
        correlation_plot: "Correlation plot image file",
        correlation_matrix: "Tab-separated correlation matrix values"
    }
    topic: "genomics,transcriptomics,epigenomics,ribosome_profiling"
    species: "any"
    operation: "visualisation,statistical_calculation"
    in_sample_req: "summary_npz:sequence_report:npz"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "correlation_plot:plot:png|pdf|svg|eps,correlation_matrix:sequence_report:tsv"
    out_ref: "none"
  }

  parameter_meta {
    summary_npz: "Compressed numpy array from multiBamSummary or multiBigwigSummary"
    sample_name: "Name identifier for the output files"
    correlation_method: "Correlation method: spearman or pearson"
    plot_type: "Plot type: heatmap or scatterplot"
    plot_format: "Output image format: png, pdf, svg, or eps"
    extra_args: "Additional arguments to pass to plotCorrelation"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File summary_npz
    String sample_name
    String correlation_method = "spearman"
    String plot_type = "heatmap"
    String plot_format = "png"
    String extra_args = ""
    Int cpu_cores = 1
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    plotCorrelation \
      --corData ~{summary_npz} \
      --plotFile "~{sample_name}_correlation.~{plot_format}" \
      --outFileCorMatrix "~{sample_name}_correlation_matrix.tab" \
      --corMethod ~{correlation_method} \
      --whatToPlot ~{plot_type} \
      ~{extra_args}
  >>>

  output {
    File correlation_plot = "~{sample_name}_correlation.~{plot_format}"
    File correlation_matrix = "~{sample_name}_correlation_matrix.tab"
  }

  runtime {
    docker: "getwilds/deeptools:3.5.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task plot_pca {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Generates a PCA plot from multiBamSummary or multiBigwigSummary output using deepTools plotPCA"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deeptools/ww-deeptools.wdl"
    outputs: {
        pca_plot: "PCA plot image file",
        pca_data: "Tab-separated PCA data with eigenvalues and coordinates"
    }
    topic: "genomics,transcriptomics,epigenomics,ribosome_profiling"
    species: "any"
    operation: "visualisation,statistical_calculation"
    in_sample_req: "summary_npz:sequence_report:npz"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "pca_plot:plot:png|pdf|svg|eps,pca_data:sequence_report:tsv"
    out_ref: "none"
  }

  parameter_meta {
    summary_npz: "Compressed numpy array from multiBamSummary or multiBigwigSummary"
    sample_name: "Name identifier for the output files"
    plot_format: "Output image format: png, pdf, svg, or eps"
    extra_args: "Additional arguments to pass to plotPCA"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    File summary_npz
    String sample_name
    String plot_format = "png"
    String extra_args = ""
    Int cpu_cores = 1
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    plotPCA \
      --corData ~{summary_npz} \
      --plotFile "~{sample_name}_pca.~{plot_format}" \
      --outFileNameData "~{sample_name}_pca_data.tab" \
      ~{extra_args}
  >>>

  output {
    File pca_plot = "~{sample_name}_pca.~{plot_format}"
    File pca_data = "~{sample_name}_pca_data.tab"
  }

  runtime {
    docker: "getwilds/deeptools:3.5.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task plot_fingerprint {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "Generates a fingerprint plot to assess ChIP enrichment quality using deepTools plotFingerprint"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-deeptools/ww-deeptools.wdl"
    outputs: {
        fingerprint_plot: "Fingerprint plot image file",
        quality_metrics: "Tab-separated file with quality metrics including JSD and CHANCE divergence"
    }
    topic: "epigenomics,data_quality_management"
    species: "any"
    operation: "quality_control"
    in_sample_req: "bam_files:nucleic_acid_sequence_alignment:bam,bai_files:data_index:bai"
    in_sample_opt: "none"
    in_ref_req: "none"
    in_ref_opt: "none"
    out_sample: "fingerprint_plot:plot:png|pdf|svg|eps,quality_metrics:sequence_report:tsv"
    out_ref: "none"
  }

  parameter_meta {
    bam_files: "Array of sorted and indexed BAM files to analyze"
    bai_files: "Array of BAM index files (.bai) corresponding to the BAM files"
    sample_name: "Name identifier for the output files"
    labels: "Optional space-separated labels for each BAM file in the plot"
    plot_format: "Output image format: png, pdf, svg, or eps"
    extra_args: "Additional arguments to pass to plotFingerprint"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
  }

  input {
    Array[File] bam_files
    Array[File] bai_files
    String sample_name
    String? labels
    String plot_format = "png"
    String extra_args = ""
    Int cpu_cores = 4
    Int memory_gb = 8
  }

  command <<<
    set -eo pipefail

    plotFingerprint \
      --bamfiles ~{sep=" " bam_files} \
      --plotFile "~{sample_name}_fingerprint.~{plot_format}" \
      --outQualityMetrics "~{sample_name}_quality_metrics.tab" \
      ~{if defined(labels) then "--labels " + labels else ""} \
      --numberOfProcessors ~{cpu_cores} \
      ~{extra_args}
  >>>

  output {
    File fingerprint_plot = "~{sample_name}_fingerprint.~{plot_format}"
    File quality_metrics = "~{sample_name}_quality_metrics.tab"
  }

  runtime {
    docker: "getwilds/deeptools:3.5.6"
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
