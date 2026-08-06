#!/usr/bin/env Rscript

library(scran)
library(scater)
library(SingleCellExperiment)
library(ggplot2)

# Parse simple --key=value arguments (no optparse dependency; not present
# in the getwilds/scran image)
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
  mito_pattern = "^MT-",
  nmads = "3.0",
  min_mean = "0.1",
  n_hvgs = "2000",
  output_prefix = NULL
)

opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
if (is.null(opt$output_prefix)) opt$output_prefix <- opt$sample_name
opt$nmads <- as.double(opt$nmads)
opt$min_mean <- as.double(opt$min_mean)
opt$n_hvgs <- as.integer(opt$n_hvgs)

set.seed(4)


##################
## 1. Load data ##
##################

message("Loading SingleCellExperiment from: ", opt$sce_rds)
counts_sce <- readRDS(opt$sce_rds)
counts_sce <- counts_sce[!duplicated(rownames(counts_sce)), ]
message("Cells loaded: ", ncol(counts_sce))


##################
## 2. QC filter ##
##################

is_mito <- grepl(opt$mito_pattern, rownames(counts_sce))
message("Mitochondrial genes detected: ", sum(is_mito))

qc_stats <- perCellQCMetrics(counts_sce, subsets = list(Mito = is_mito))
qc_filter <- quickPerCellQC(qc_stats,
                            sub.fields = "subsets_Mito_percent",
                            nmads = opt$nmads)

colData(counts_sce) <- cbind(colData(counts_sce), qc_stats)
counts_sce$discard <- qc_filter$discard

# QC scatter plot of library size vs. detected features, colored by discard status
qc_plot <- plotColData(counts_sce,
                       x = "sum",
                       y = "detected",
                       colour_by = "discard") +
  ggtitle(paste("Per-Cell QC for Sample:", opt$sample_name)) +
  xlab("Total UMI Counts") +
  ylab("Detected Genes")

ggsave(paste0(opt$output_prefix, "_qc.png"),
       plot = qc_plot,
       dpi = 300,
       width = 7,
       height = 6,
       device = "png")

sce <- counts_sce[, !qc_filter$discard]
message("Cells after QC: ", ncol(sce))


#################################################
## 3. Normalize with scran deconvolution       ##
#################################################

clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters, min.mean = opt$min_mean)
sce <- logNormCounts(sce)

sf_plot <- ggplot(data.frame(size_factor = sizeFactors(sce)),
                  aes(x = size_factor)) +
  geom_histogram(bins = 40) +
  ggtitle(paste("Size Factor Distribution for Sample:", opt$sample_name)) +
  xlab("Size Factor") +
  ylab("Number of Cells")

ggsave(paste0(opt$output_prefix, "_size_factors.png"),
       plot = sf_plot,
       dpi = 300,
       width = 6,
       height = 5,
       device = "png")


#####################################################
## 4. Model variance and pick highly variable genes ##
#####################################################

dec <- modelGeneVar(sce, min.mean = opt$min_mean)
top_hvgs <- getTopHVGs(dec, n = opt$n_hvgs)

hvg_df <- as.data.frame(dec)
hvg_df$gene <- rownames(hvg_df)
write.csv(hvg_df[order(-hvg_df$bio), ],
          paste0(opt$output_prefix, "_hvgs.csv"),
          row.names = FALSE)

mean_var_plot <- ggplot(hvg_df, aes(x = mean, y = total)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(y = tech), color = "dodgerblue", linewidth = 1) +
  ggtitle(paste("Mean-Variance Trend for Sample:", opt$sample_name)) +
  xlab("Mean Log-Expression") +
  ylab("Variance")

ggsave(paste0(opt$output_prefix, "_mean_variance.png"),
       plot = mean_var_plot,
       dpi = 300,
       width = 6,
       height = 5,
       device = "png")


############################
## 5. PCA and save object ##
############################

sce <- runPCA(sce, subset_row = top_hvgs)

saveRDS(sce, file = paste0(opt$output_prefix, ".rds"))
message("Done. Output prefix: ", opt$output_prefix)
