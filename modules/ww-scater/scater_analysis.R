#!/usr/bin/env Rscript

library(scater)
library(SingleCellExperiment)
library(ggplot2)

# Parse simple --key=value arguments (no optparse dependency; not present
# in the getwilds/scater image)
parse_args <- function(args, defaults) {
  opt <- defaults
  for (arg in args) {
    kv <- sub("^--", "", arg)
    key <- sub("=.*$", "", kv)
    value <- sub("^[^=]*=", "", kv)
    opt[[key]] <- value
  }
  opt
}

defaults <- list(
  sce_rds = NULL,
  sample_name = NULL,
  n_pcs = "50",
  n_hvgs = "2000",
  random_seed = "100",
  output_prefix = NULL
)

opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
if (is.null(opt$output_prefix)) opt$output_prefix <- opt$sample_name
opt$n_pcs <- as.integer(opt$n_pcs)
opt$n_hvgs <- as.integer(opt$n_hvgs)
opt$random_seed <- as.integer(opt$random_seed)

set.seed(opt$random_seed)


##################
## 1. Load data ##
##################

message("Loading SingleCellExperiment from: ", opt$sce_rds)
sce <- readRDS(opt$sce_rds)
message("Cells loaded: ", ncol(sce))

if (!"logcounts" %in% assayNames(sce)) {
  stop("Input SingleCellExperiment has no logcounts assay; ",
       "normalize (e.g. with ww-scran) before running scater")
}


########################################################
## 2. Highly variable genes for dimensionality reduction ##
########################################################

# Rank genes by variance of log-expression to pick the top HVGs
# (self-contained: does not assume upstream HVG selection, e.g. scran's)
gene_var <- apply(logcounts(sce), 1, var)
n_hvgs <- min(opt$n_hvgs, length(gene_var))
top_hvgs <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_hvgs)]
message("HVGs used for dimensionality reduction: ", length(top_hvgs))


############################################
## 3. PCA, UMAP, and QC/reduced-dim plots ##
############################################

sce <- runPCA(sce, subset_row = top_hvgs, ncomponents = opt$n_pcs)
sce <- runUMAP(sce, dimred = "PCA")

pca_plot <- plotReducedDim(sce, dimred = "PCA") +
  ggtitle(paste("PCA:", opt$sample_name))
ggsave(paste0(opt$output_prefix, "_pca.png"),
       plot = pca_plot, dpi = 300, width = 6, height = 5, device = "png")

umap_plot <- plotReducedDim(sce, dimred = "UMAP") +
  ggtitle(paste("UMAP:", opt$sample_name))
ggsave(paste0(opt$output_prefix, "_umap.png"),
       plot = umap_plot, dpi = 300, width = 6, height = 5, device = "png")

qc_metrics <- perCellQCMetrics(sce)
colData(sce) <- cbind(colData(sce), qc_metrics)

qc_plot <- plotColData(sce, x = "sum", y = "detected",
                        colour_by = "total") +
  ggtitle(paste("Per-Cell QC:", opt$sample_name)) +
  xlab("Total UMI Counts") +
  ylab("Detected Genes")
ggsave(paste0(opt$output_prefix, "_qc.png"),
       plot = qc_plot, dpi = 300, width = 7, height = 6, device = "png")

saveRDS(sce, file = paste0(opt$output_prefix, ".rds"))
message("Done. Output prefix: ", opt$output_prefix)
