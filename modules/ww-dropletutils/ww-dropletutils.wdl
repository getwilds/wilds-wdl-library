## WILDS WDL Module: ww-dropletutils
## Description: Module for loading 10x Genomics feature-barcode matrices and calling
## real cells from empty droplets using Bioconductor DropletUtils. Produces a
## SingleCellExperiment object (RDS) suitable as input to downstream modules such as
## ww-scran, ww-scater, and ww-singler.

version 1.0

#### TASK DEFINITIONS ####

task read10x_counts {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Load a 10x Genomics feature-barcode matrix (H5 format) into a SingleCellExperiment object using DropletUtils::read10xCounts"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-dropletutils/ww-dropletutils.wdl"
    outputs: {
        sce_rds: "SingleCellExperiment object (RDS) containing the loaded counts, gene, and barcode metadata"
    }
    topic: "transcriptomics,single_cell"
    species: "human,mouse,eukaryote"
    operation: "data_loading"
    input_sample_required: "h5_matrix:gene_expression_matrix:h5"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "sce_rds:gene_expression_matrix:rds"
    output_reference: "none"
  }

  parameter_meta {
    h5_matrix: "10x Genomics feature-barcode matrix in HDF5 format (filtered or raw)"
    sample_name: "Sample name used as output file prefix and SCE column metadata"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
    docker_image: "Docker image to use for this task"
  }

  input {
    File h5_matrix
    String sample_name
    Int cpu_cores = 2
    Int memory_gb = 8
    String docker_image = "getwilds/dropletutils:1.32.0"
  }

  command <<<
    set -euo pipefail

    Rscript -e '
      library(DropletUtils)

      sce <- read10xCounts("~{h5_matrix}", col.names = TRUE)
      sce$Sample <- "~{sample_name}"
      counts(sce) <- as(counts(sce), "CsparseMatrix")

      saveRDS(sce, file = "~{sample_name}.sce.rds")
    '
  >>>

  output {
    File sce_rds = "~{sample_name}.sce.rds"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}

task empty_drops_filter {
  meta {
    author: "WILDS Team"
    email: "wilds@fredhutch.org"
    description: "Call real cells from a raw (unfiltered) 10x feature-barcode matrix using DropletUtils::emptyDrops, filter the SingleCellExperiment to called cells, and produce a barcode-rank QC plot"
    url: "https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/main/modules/ww-dropletutils/ww-dropletutils.wdl"
    outputs: {
        filtered_sce_rds: "SingleCellExperiment object (RDS) filtered to barcodes called as cells by emptyDrops",
        empty_drops_csv: "Per-barcode emptyDrops statistics (Total, LogProb, PValue, FDR) for all tested barcodes",
        barcode_rank_pdf: "Barcode-rank plot (log-log UMI count vs. rank) with knee and inflection points marked"
    }
    topic: "transcriptomics,single_cell"
    species: "human,mouse,eukaryote"
    operation: "data_handling"
    input_sample_required: "raw_h5_matrix:gene_expression_matrix:h5"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "filtered_sce_rds:gene_expression_matrix:rds,empty_drops_csv:gene_report:csv,barcode_rank_pdf:plot:pdf"
    output_reference: "none"
  }

  parameter_meta {
    raw_h5_matrix: "Raw (unfiltered) 10x Genomics feature-barcode matrix in HDF5 format, containing all detected barcodes"
    sample_name: "Sample name used as output file prefix and SCE column metadata"
    fdr_threshold: "False discovery rate threshold below which a barcode is called a real cell"
    lower_umi_threshold: "UMI count threshold below which barcodes are assumed to be empty droplets and used to estimate the ambient RNA profile"
    random_seed: "Random seed for emptyDrops Monte Carlo p-value computation, for reproducibility"
    cpu_cores: "Number of CPU cores allocated for the task"
    memory_gb: "Memory allocated for the task in GB"
    docker_image: "Docker image to use for this task"
  }

  input {
    File raw_h5_matrix
    String sample_name
    Float fdr_threshold = 0.01
    Int lower_umi_threshold = 100
    Int random_seed = 100
    Int cpu_cores = 2
    Int memory_gb = 8
    String docker_image = "getwilds/dropletutils:1.32.0"
  }

  command <<<
    set -euo pipefail

    Rscript -e '
      library(DropletUtils)

      set.seed(~{random_seed})

      sce <- read10xCounts("~{raw_h5_matrix}", col.names = TRUE)
      sce$Sample <- "~{sample_name}"

      ranks <- barcodeRanks(counts(sce))
      pdf("~{sample_name}.barcode_rank.pdf")
      plot(ranks$rank, ranks$total, log = "xy", xlab = "Rank", ylab = "Total UMI count",
           main = "~{sample_name} Barcode Rank Plot")
      abline(h = metadata(ranks)$knee, col = "dodgerblue", lty = 2)
      abline(h = metadata(ranks)$inflection, col = "forestgreen", lty = 2)
      legend("bottomleft", legend = c("Knee", "Inflection"),
             col = c("dodgerblue", "forestgreen"), lty = 2)
      dev.off()

      empty <- emptyDrops(counts(sce), lower = ~{lower_umi_threshold})
      write.csv(as.data.frame(empty), file = "~{sample_name}.empty_drops.csv", row.names = TRUE)

      is_cell <- which(empty$FDR <= ~{fdr_threshold})
      filtered_sce <- sce[, is_cell]
      counts(filtered_sce) <- as(counts(filtered_sce), "CsparseMatrix")

      saveRDS(filtered_sce, file = "~{sample_name}.filtered_sce.rds")
    '
  >>>

  output {
    File filtered_sce_rds = "~{sample_name}.filtered_sce.rds"
    File empty_drops_csv = "~{sample_name}.empty_drops.csv"
    File barcode_rank_pdf = "~{sample_name}.barcode_rank.pdf"
  }

  runtime {
    docker: docker_image
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
  }
}
