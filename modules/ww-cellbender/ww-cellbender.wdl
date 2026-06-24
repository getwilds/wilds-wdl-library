## WILDS WDL Module: ww-cellbender
## Description: Module for removing ambient RNA and barcode swapping from single-cell RNA-seq data using CellBender.
## CellBender's remove-background command uses a deep generative model to separate
## real cell signal from technical artifacts in UMI-based scRNA-seq, snRNA-seq, and CITE-seq data.

version 1.0

#### TASK DEFINITIONS ####

task remove_background {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Remove ambient RNA contamination and barcode swapping artifacts from a 10x Genomics raw feature-barcode matrix using CellBender remove-background"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-cellbender/ww-cellbender.wdl"
    outputs: {
        output_h5: "Full cleaned count matrix in H5 format (all input barcodes retained)",
        filtered_h5: "Filtered cleaned count matrix containing only barcodes with >50% cell probability",
        report_html: "Interactive HTML report with training diagnostics and quality plots",
        summary_pdf: "PDF summary of the inference procedure",
        log_file: "Run log with diagnostics and parameter summary",
        cell_barcodes_csv: "Plain-text list of cell barcodes exceeding 50% cell probability threshold",
        metrics_csv: "Run metrics CSV (one row per FPR value)",
        checkpoint_tar: "Trained model checkpoint tarball (can be used to resume or inspect training)"
    }
    topic: "transcriptomics,single_cell"
    species: "human,mouse,eukaryote"
    operation: "data_handling"
    input_sample_required: "input_h5:gene_expression_matrix:h5"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "output_h5:gene_expression_matrix:h5,filtered_h5:gene_expression_matrix:h5,report_html:quality_control_report:html,summary_pdf:plot:pdf,cell_barcodes_csv:gene_report:csv,metrics_csv:gene_report:csv"
    output_reference: "none"
  }

  parameter_meta {
    input_h5: "Raw feature-barcode matrix HDF5 file from Cell Ranger (raw_feature_bc_matrix.h5)"
    sample_name: "Sample name used as output file prefix"
    expected_cells: "Optional estimated number of real cells in the sample; helps CellBender set priors"
    total_droplets_included: "Optional total number of droplets to include in the analysis (raw barcodes ranked by UMI count); defaults to CellBender auto-selection"
    fpr: "False positive rate target(s) for the output matrix, space-separated if multiple (e.g. '0.01 0.05')"
    epochs: "Number of training epochs"
    learning_rate: "Base learning rate for the variational inference optimizer"
    model: "CellBender model variant: naive, simple, ambient, swapping, or full"
    low_count_threshold: "Droplets with total UMI count below this value are excluded from analysis"
    exclude_feature_types: "Optional space-separated list of feature types to exclude (e.g. 'Antibody Capture')"
    gpu_enabled: "Enable GPU acceleration (adds --cuda flag and requests 1 GPU in runtime); set to false for CPU-only execution"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
    docker_image: "Docker image to use for this task"
  }

  input {
    File input_h5
    String sample_name
    Int? expected_cells
    Int? total_droplets_included
    String fpr = "0.01"
    Int epochs = 150
    Float learning_rate = 0.0001
    String model = "full"
    Int low_count_threshold = 5
    String? exclude_feature_types
    Boolean gpu_enabled = true
    Int cpu_cores = 4
    Int memory_gb = 32
    String docker_image = "getwilds/cellbender:0.3.2"
  }

  command <<<
    set -eo pipefail

    cellbender remove-background \
      --input "~{input_h5}" \
      --output "~{sample_name}_out.h5" \
      --model ~{model} \
      --epochs ~{epochs} \
      --learning-rate ~{learning_rate} \
      --fpr ~{fpr} \
      --low-count-threshold ~{low_count_threshold} \
      ~{"--expected-cells " + expected_cells} \
      ~{"--total-droplets-included " + total_droplets_included} \
      ~{"--exclude-feature-types " + exclude_feature_types} \
      ~{if gpu_enabled then "--cuda" else ""} \
      --cpu-threads ~{cpu_cores}
  >>>

  output {
    File output_h5             = "~{sample_name}_out.h5"
    File filtered_h5           = "~{sample_name}_out_filtered.h5"
    File report_html           = "~{sample_name}_out_report.html"
    File summary_pdf           = "~{sample_name}_out.pdf"
    File log_file              = "~{sample_name}_out.log"
    File cell_barcodes_csv     = "~{sample_name}_out_cell_barcodes.csv"
    File metrics_csv           = "~{sample_name}_out_metrics.csv"
    File checkpoint_tar        = "ckpt.tar.gz"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
    gpus: if gpu_enabled then "1" else "0"
  }
}
